#!/usr/bin/env nextflow
// this script was optimized apptainer in any compurational environment. 
nextflow.enable.dsl=2
//params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.db_baseDir = "/opt/database"

params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.outdir_Base = "${launchDir}" // modify this if you want output files stored any other directory. 
params.inputDir = "${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT/bins_quality_passedFinal"
params.outdir = "${params.outdir_Base}/results/metagenome/COMPARATIVE_ANNOTATION"
//params.metadata ="${params.inputDir}/../final_report/quality_taxonomy_combined_final.csv"
params.metadata ="${launchDir}/selected_metadata.csv"
// metadata loc should be modified. 
params.metacol =""
params.samplecol = 1
params.module_completeness = 0.5
//params.kofam_db = 
params.cpus = 40
params.pan_identity=0.8
params.pan_coverage=0.8
params.kingdom="bacteria"


if (params.metacol == null || params.metacol == 0) {
    log.error 'Metadata column parameter --metacol cannot be empty.\n' + 
               'You can specify by --metacol {column number} ' + 
               'You can check {column number} by using  $head -n1 selected_metadata.csv | tr "," "\n" | nl'
    System.exit(1) // Exits the script with an error status
}

if (!new File(params.inputDir).exists()) {
    error """Input directory does not exist: ${params.inputDir}.\n
    Please specify a valid directory with --inputDir.\n
    
    You would select your genomes with genome metadata file and script \n   : select_genomes_toCOMPARATIVE_ANNOTATION.sh. \n
    
    Your combined_metadata.csv   could be generated using python script.
    Suppose second column of ### metageome_metadata.csv ### is WMS accession ID extracted from paired end metagenome file name 
    python combine_metadata_WMS_genome.py -g  quality_taxonomy_combined.csv -m metageome_metadata.csv -o combined_metadata.csv -a 2

    sh select_genomes_toCOMPARATIVE_ANNOTATION.sh \"s__Bacteroides uniformis\"

    You need to also provide metadata to program if you do not want to use select_genomes_toCOMPARATIVE_ANNOTATION.sh
"""
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : bin directory."
    }
}

//######################



log.info """\
    Pangenome analysis, species functional annotation 
    =================================================
    PPPanGGOLiN option          : Identity = ${params.pan_identity}, 
                                Coverage = ${params.pan_coverage}
                                modify it if you want by specifying --pan_identity [0~1] --pan_coverage [0~1] in command line 
    KEGG module present threshold : ${params.module_completeness}
                                modify threshold if you want to change it : --module_completeness [0~1]
    
    """
    .stripIndent(true)



process run_skani {
    publishDir "${params.outdir}/ani" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(genome_dir)
    path metadata_file 
    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    skani triangle *.fa -t ${task.cpus}  --full-matrix | tail -n +2 | sed -e 's/.fa//g' > skani_fullmatrix
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/skani_visualization.R \
    -i skani_fullmatrix -m ${metadata_file}  -mc ${params.metacol} -out heatmap_skani.pdf -out_html skani_interactive

    """
}




process run_prokka {
    publishDir "${params.outdir}/prokka" , mode : 'copy'
    container = 'docker://staphb/prokka:latest'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    //cpus params.cpus

    input:
    tuple val(id), path(genome)
    output:
    path "*"

    script:
    """
    #eval "\$(micromamba shell hook --shell bash)"
    #micromamba  activate COMPARATIVE_ANNOTATION
    #export TMPDIR='/tmp/\${uuidgen()}'
    #export JAVA_OPTS="-XX:-UsePerfData -Djava.io.tmpdir=\$TMPDIR"

    prokka --prefix ${id} --noanno --cpus 2 ${genome} --compliant 
    #--centre  MGSSB 
    """
}
//prokka --prefix ${id} --noanno --cpus ${task.cpus} ${launchDir}/${params.inputDir}/${genome} --centre  MGSSB --compliant 
process run_panaroo {
    publishDir "${params.outdir}" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(prokka_dir)
    output:
    path "*"

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    panaroo -i */*.gff -o panaroo_result --clean-mode moderate --remove-invalid-genes \
    -c 0.9 -f 0.5 --merge_paralogs --threads ${task.cpus}
    """
}


process run_ppanggolin {
    publishDir "${params.outdir}" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(prokka_dir)
    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    ls */*.gff | while read gff_file; do 
    printf "%s\t%s\n" "\$(basename "\$gff_file" .gff)" "\$(realpath "\$gff_file")"
    done > gff_paths.tsv
    ppanggolin 
    ppanggolin all --anno gff_paths.tsv -o ppanggolin_result --cpu ${task.cpus}\
    --coverage ${params.pan_identity} --identity ${params.pan_coverage} \
    --kingdom bacteria --rarefaction
    cd ppanggolin_result 
    ppanggolin fasta -p pangenome.h5 --prot_families all -o pan_genes  --genes all 
    python ${params.scripts_baseDir}/make_gene_count_table_ppanggolin205.py ./
    # genome_id_linking.tsv, gene_count_matrix.tsv is generated 
    """
}
//#panaroo -i */*.gff -o ppanggorin_result --clean-mode moderate --remove-invalid-genes \
//    #-c 0.9 -f 0.5 --merge_paralogs --threads ${params.cpus}
//find "\$(pwd)" -type f -name "*.gff" -exec sh -c 'echo -e "\$(basename {} .gff)\t\$(realpath "{}")"' \; > gff_paths.tsv    

process run_genePA_cluster {
    publishDir "${params.outdir}/genePA_cluster" , mode : 'copy'
    cpus params.cpus
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    input:
    path(ppanggolin_dir)
    path metadata 
    output:
    path("pcoa_plot_interactive.html")    

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate scoary-2
    python ${params.scripts_baseDir}/genome_PA_stat.py --cpus ${task.cpus}

    """

}

process run_kofamscan {
    publishDir "${params.outdir}/kofamscan" , mode : 'copy'
    cpus params.cpus
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(ppanggolin_dir)
    path metadata_file 

    output:
    path("ko_matrix.csv")
    path("KEGG_module_visualization_shiny")
    path("KEGG_module_completeness.csv")
    path("heatmap_KEGG.pdf")
    //path("pangene_kofamscan_KO.tsv")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    #ppanggolin fasta -p ppanggolin_result/pangenome.h5 -o pan_faa  --genes all 
    #transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    #sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    #micromamba activate kofamscan
    #time /scratch/tools/microbiome_analysis/program/kofam_scan-1.3.0/exec_annotation
    
    #sed -i "s|/scratch/tools/microbiome_analysis/database/kofam/27Nov2023/profiles|/opt/database/kofam/27Nov2023/profiles|g" ${params.db_baseDir}/kofam/27Nov2023/config-template.yml
    #sed -i "s|/scratch/tools/microbiome_analysis/database/kofam/27Nov2023/ko_list|/opt/database/kofam/27Nov2023/ko_list|g" ${params.db_baseDir}/kofam/27Nov2023/config-template.yml
    echo running kofamscan
    time exec_annotation \
    -o kofamscan_result --cpu ${task.cpus}  --e-value 0.00001 \
    -c ${params.db_baseDir}/kofam/27Nov2023/config-template.yml \
    -f mapper ppanggolin_result/pan_genes/all_protein_families.faa
    
    head kofamscan_result
    #sed -i  's/_[0-9]*\t/\t/g' kofamscan_result # this is for panaroo
    #sed -i 's/_[0-9]*\$//g' kofamscan_result # this is for panaroo
    #micromamba activate COMPARATIVE_ANNOTATION
    
    # this is for panaroo
    #python ${params.scripts_baseDir}/kofam_to_geneID_KO.py \
    #-i panaroo_result/gene_presence_absence.csv -o output_gene_pa_ko.csv

    # for ppanggolin205
    python ${params.scripts_baseDir}/kofam_to_geneID_KO_ppanggolin205.py \
    -i ppanggolin_result/matrix.csv -o output_gene_pa_ko.csv

    #python ${params.scripts_baseDir}/KO_Genome_long2matrix.py 
    python ${params.scripts_baseDir}/KO_Genome_long2matrix_ppanggolin205.py \
    -i output_gene_pa_ko.csv  -o ko_matrix.csv

    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/KO_module_visualization_apptainer.R \
    -i  ko_matrix.csv -m ${metadata_file}  -n ${params.metacol} -mc ${params.module_completeness}\
    -out_table KEGG_module_completeness.csv -out_html KEGG_module_visualization_shiny
    """
}

process run_VFDB {
    publishDir "${params.outdir}/VFDB" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 


    output:
    path("pangene_vfdb_result.txt")
    path("gene_PA_VFDB_added.csv")
    path("heatmap_VFDB.pdf")
    path("VFDB_interactive")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    #transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    #sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    #ppanggolin_result/pan_genes/all_protein_families.faa
    micromamba activate rgi603
    diamond blastp -d ${params.db_baseDir}/VFDB/VFDB_setB_prot_Aug2023.dmnd \
    -q ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_vfdb_result.txt -f 6 --max-target-seqs 1 --id 50 \
    --subject-cover 80 -e 1e-10 

    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt -g ppanggolin_result/gene_presence_absence.Rtab \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_PA_VFDB_added.csv

    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i  gene_PA_VFDB_added.csv -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive

    """
}


process run_rgi_CARD {
    publishDir "${params.outdir}/CARD" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 

    output:
    path("pangene_rgi_CARD_result.txt")
    path("gene_PA_CARD_added.csv")
    path("heatmap_CARD.pdf")
    path("CARD_interactive")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    #transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    #sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    #sed -i "s/\\*\$//g" pan_genome_reference.faa
    micromamba activate rgi603
    rgi load -i ${params.db_baseDir}/CARD/card.json --local 
    rgi main -i ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_rgi_CARD_result --clean  \
    -t protein -n ${task.cpus} --include_nudge --local 

    python ${params.scripts_baseDir}/add_rgi_togenePA.py -i pangene_rgi_CARD_result.txt\
    -o gene_PA_CARD_added.csv -r ${params.db_baseDir}/CARD/aro_index.tsv \
    -gpa ppanggolin_result/gene_presence_absence.Rtab

    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/CARD_visualization.R \
    -i  gene_PA_CARD_added.csv -m ${metadata_file}  -mc ${params.metacol} -out heatmap_CARD.pdf -out_html CARD_interactive
    """
}

process run_dbCAN {
    publishDir "${params.outdir}/dbCAN" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 

    output:
    path("dbcan_gene_count_matrix.csv")
    path("dbcan_hmmerfamily_count_matrix.csv")
    path("heatmap_dbCAN.pdf")
    path("dbCAN_interactive")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate run_dbcan4
    export OMP_NUM_THREADS=${task.cpus}

    time run_dbcan ppanggolin_result/pan_genes/all_protein_families.faa  protein --dia_cpu ${task.cpus}  \
    --hmm_cpu ${task.cpus}  --out_dir db_can_out  --db_dir ${params.db_baseDir}/dbCAN

    cut -f1,3 db_can_out/overview.txt | grep -Pv "\\t-\$" |  awk -F'\\t' '{split(\$2,a,"+"); for (i in a) {print \$1 "\\t" a[i];}}' | \
    awk -F"\\t" '{split(\$2,a,"_"); \$2=a[1]; print \$1 "\\t" \$2}'  | sed -e "s/(.*//" > dbcan_family.tsv

    python ${params.scripts_baseDir}/dbcan_result_parse.py 
    micromamba activate R432_environment
    Rscript  ${params.scripts_baseDir}/dbcan_visualization.R -i dbcan_hmmerfamily_count_matrix.csv  \
    -m ${metadata_file}  -mc ${params.metacol} -out heatmap_dbCAN.pdf -out_html dbCAN_interactive
    """
}

process run_eggNOG {
    publishDir "${params.outdir}/eggNOG" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 
    output:
    path("eggnog_mmseqs.emapper*")
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate eggnog-mapper2112
    emapper.py -m mmseqs --data_dir ${params.db_baseDir}/eggNOG5 --output eggnog_mmseqs  --cpu ${task.cpus} \
    -i  ppanggolin_result/pan_genes/all_protein_families.faa --itype proteins
    """
}

process run_scoary2 {
    publishDir "${params.outdir}/scoary2" , mode : 'copy'
    cpus 8

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 
    output:
    path("scoary_out")
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate scoary-2
    python ${params.scripts_baseDir}/make_scoary_binary_trait.py -i ${metadata_file} -m ${params.metacol} -s ${params.samplecol}
    # genome_traits.tsv is generated

    # Set MPLCONFIGDIR to a writable directory
    export MPLCONFIGDIR=\$PWD/mpl_config
    mkdir -p \$MPLCONFIGDIR

    # Set NUMBA_CACHE_DIR to a writable directory
    export NUMBA_CACHE_DIR=\$PWD/numba_cache
    mkdir -p \$NUMBA_CACHE_DIR
    #export CONFINT_DB=\$PWD/confint_cache
    #mkdir -p \$CONFINT_DB
    export CONFINT_DB=$PWD/CONFINT_DB

    scoary2 --genes ppanggolin_result/gene_count_matrix.tsv  \
    --gene-data-type "gene-count:\\t" --traits genome_traits.tsv --trait-data-type "binary:\\t" --n-permut 1000 \
    --n-cpus ${task.cpus} --outdir scoary_out
    """
}


/*
process run_kofamscan 
*/
workflow{
    bins_ch = Channel.fromFilePairs("${params.inputDir}/*.fa",size:1)
    
    prokka_result = run_prokka(bins_ch)
    prokka_collect = prokka_result.collect()
    ppanggolin_result = run_ppanggolin(prokka_collect)
    metadata_ch = Channel.fromPath("${params.metadata}")
    run_genePA_cluster(ppanggolin_result,metadata_ch)
    run_kofamscan(ppanggolin_result,metadata_ch)
    run_VFDB(ppanggolin_result,metadata_ch)
    run_rgi_CARD(ppanggolin_result,metadata_ch)
    run_dbCAN(ppanggolin_result,metadata_ch)
    run_eggNOG(ppanggolin_result,metadata_ch)
    run_scoary2(ppanggolin_result,metadata_ch)
    //metadata_ch = file(params.metadata)
    allbins=Channel.fromPath("${params.inputDir}/*.fa")
    run_skani(allbins.collect(),metadata_ch)

    //run_VFDB(panaroo_result,metadata_ch)
    //run_rgi_CARD(panaroo_result,metadata_ch)
}
