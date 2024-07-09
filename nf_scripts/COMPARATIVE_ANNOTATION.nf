#!/usr/bin/env nextflow
// this script was optimized for K-BDS hosted by KISTI in South Korea. 
nextflow.enable.dsl=2
params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.outdir_Base = "$HOME" // modify this if you want output files stored any other directory. 
params.inputDir = "${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT/bins_quality_passedFinal"
params.outdir = "${params.outdir_Base}/results/metagenome/COMPARATIVE_ANNOTATION"
params.metadata ="${params.inputDir}/../final_report/quality_taxonomy_combined_final.csv"
params.metacol =""
params.module_completeness = 0.5
//params.kofam_db = 
params.cpus = 64
params.pan_identity=0.8
params.pan_coverage=0.8


if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir."
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : bin directory."
    }
}

process run_skani {
    publishDir "${params.outdir}/ani" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(genome_dir)
    path metadata_file 
    output:
    path "*"

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    skani triangle *.fa -t 32  --full-matrix | tail -n +2 | sed -e 's/.fa//g' > skani_fullmatrix
    conda activate R432_environment
    Rscript ${params.scripts_baseDir}/skani_visualization.R \
    -i skani_fullmatrix -m ${metadata_file}  -mc ${params.metacol} -out heatmap_skani.pdf -out_html skani_interactive

    """
}




process run_prokka {
    publishDir "${params.outdir}/prokka" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    tuple val(id), path(genome)
    output:
    path "*"

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    prokka --prefix ${id} --noanno --cpus ${params.cpus} ${genome} --centre  MGSSB --compliant 
    """
}

process run_panaroo {
    publishDir "${params.outdir}" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(prokka_dir)
    output:
    path "*"

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    panaroo -i */*.gff -o panaroo_result --clean-mode moderate --remove-invalid-genes \
    -c 0.9 -f 0.5 --merge_paralogs --threads ${params.cpus}
    """
}

process run_ppanggolin {
    publishDir "${params.outdir}" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(prokka_dir)
    output:
    path "*"

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    find "$(pwd)" -type f -name "*.gff" -exec sh -c 'echo -e "$(basename "{}" .gff)\t$(realpath "{}")"' \; > gff_paths.tsv
    ppanggolin 
    ppanggolin all --anno gff_paths.tsv -o ppanggolin_result  --cpu ${params.cpus}\
    --coverage ${params.pan_identity} --identity ${params.pan_coverage} \ 
     --kingdom bacteria --rarefaction
    #panaroo -i */*.gff -o ppanggorin_result --clean-mode moderate --remove-invalid-genes \
    #-c 0.9 -f 0.5 --merge_paralogs --threads ${params.cpus}
    """
}



process run_kofamscan {
    publishDir "${params.outdir}/kofamscan" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(panaroo_dir)
    path metadata_file 

    output:
    path("ko_matrix.csv")
    path("KEGG_module_visualization_shiny")
    path("KEGG_module_completeness.csv")
    path("heatmap_KEGG.pdf")
    //path("pangene_kofamscan_KO.tsv")

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    ppanggolin -p ppanggolin_result/pangenome.h5 -o pan_faa  --genes all 

    transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    conda activate kofamscan
    time /scratch/tools/microbiome_analysis/program/kofam_scan-1.3.0/exec_annotation \
    -o kofamscan_result --cpu ${params.cpus} \
    -c ${params.db_baseDir}/kofam/27Nov2023/config-template.yml \
    -f mapper pan_genome_reference.faa

    #sed -i  's/_[0-9]*\t/\t/g' kofamscan_result
    #sed -i 's/_[0-9]*\$//g' kofamscan_result
    conda activate COMPARATIVE_ANNOTATION
    
    python ${params.scripts_baseDir}/kofam_to_geneID_KO.py \
    -i panaroo_result/gene_presence_absence.csv -o output_gene_pa_ko.csv
    
    python ${params.scripts_baseDir}/KO_Genome_long2matrix.py \
    -i output_gene_pa_ko.csv  -o ko_matrix.csv

    conda activate R432_environment
    Rscript ${params.scripts_baseDir}/KO_module_visualization.R \
    -i  ko_matrix.csv -m ${metadata_file}  -n ${params.metacol} -mc ${params.module_completeness}\
    -out_table KEGG_module_completeness.csv -out_html KEGG_module_visualization_shiny
    """
}

process run_VFDB {
    publishDir "${params.outdir}/VFDB" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(panaroo_dir)
    //path("${params.metadata}")
    path metadata_file 


    output:
    path("pangene_vfdb_result.txt")
    path("gene_PA_VFDB_added.csv")
    path("heatmap_VFDB.pdf")
    path("VFDB_interactive")

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    conda activate rgi603
    diamond blastp -d ${params.db_baseDir}/VFDB/VFDB_setB_prot_Aug2023.dmnd \
    -q pan_genome_reference.faa -o pangene_vfdb_result.txt -f 6 --max-target-seqs 1 --id 50 \
    --subject-cover 80 -e 1e-10 

    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt -g panaroo_result/gene_presence_absence.Rtab \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_PA_VFDB_added.csv

    conda activate R432_environment
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i  gene_PA_VFDB_added.csv -m ${metadata_file}  -mc ${params.metacol} -out_html VFDB_interactive
    """
}


process run_rgi_CARD {
    publishDir "${params.outdir}/CARD" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(panaroo_dir)
    //path("${params.metadata}")
    path metadata_file 

    output:
    path("pangene_rgi_CARD_result.txt")
    path("gene_PA_CARD_added.csv")
    path("heatmap_CARD.pdf")
    path("CARD_interactive")

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    sed -i "s/\\*\$//g" pan_genome_reference.faa
    conda activate rgi603
    rgi load -i /scratch/tools/microbiome_analysis/database/CARD/card.json --local 
    rgi main -i pan_genome_reference.faa -o pangene_rgi_CARD_result --clean  \
    -t protein -n ${params.cpus} --include_nudge --local 

    python ${params.scripts_baseDir}/add_rgi_togenePA.py -i pangene_rgi_CARD_result.txt\
    -o gene_PA_CARD_added.csv -r ${params.db_baseDir}/CARD/aro_index.tsv \
    -gpa panaroo_result/gene_presence_absence.Rtab

    conda activate R432_environment
    Rscript ${params.scripts_baseDir}/CARD_visualization.R \
    -i  gene_PA_CARD_added.csv -m ${metadata_file}  -mc ${params.metacol} -out heatmap_CARD.pdf -out_html CARD_interactive
    """
}

/*
process run_kofamscan 
*/
workflow{
    bins_ch = Channel.fromFilePairs("${params.inputDir}/*.fa",size:1)
    
    prokka_result = run_prokka(bins_ch)
    prokka_collect = prokka_result.collect()
    panaroo_result = run_panaroo(prokka_collect)
    metadata_ch = Channel.fromPath("${params.metadata}")
    run_kofamscan(panaroo_result,metadata_ch)

    //metadata_ch = file(params.metadata)
    allbins=Channel.fromPath("${params.inputDir}/*.fa")
    run_skani(allbins.collect(),metadata_ch)

    run_VFDB(panaroo_result,metadata_ch)
    run_rgi_CARD(panaroo_result,metadata_ch)
}
