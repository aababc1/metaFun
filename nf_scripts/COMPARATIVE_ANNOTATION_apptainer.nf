#!/usr/bin/env nextflow
// this script was optimized apptainer in any compurational environment. 
nextflow.enable.dsl=2
//params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.db_baseDir = "/opt/database"

params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.outdir_Base = "${launchDir}" // modify this if you want output files stored any other directory. 
params.inputDir = "${params.outdir_Base}/results/metagenome/ASSEMBLY_BINNING/final_bins"
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

params.kofamscan_eval= 0.00001 
params.VFDB_identity=50
params.VFDB_coverage= 80
params.VFDB_e_value = 1e-10

params.CAZyme_hmm_eval=1e-15
params.CAZyme_hmm_cov=0.35

params.run_drep = true
params.drep_ani = 0.995
params.drep_cov = 0.3
params.drep_ignore_quality = true
params.drep_algorithm = "skani"


def current_time = new Date().format('yyyyMMddHHmmss')

//params.outdir = "${params.outdir_Base}/results/metagenome/COMPARATIVE_ANNOTATION"
params.run_name = current_time
params.outdir = "${params.outdir_Base}/results/metagenome/COMPARATIVE_ANNOTATION/${params.run_name}"
params.visualization_results = "${params.outdir}/visualization_results" 
params.annotation_results =  "${params.outdir}/annotation_results" 
//params.gene_matrix = null // gene_PA or gene_count

def logo_file = file("${workflow.projectDir}/../scripts/logo.txt")
//def logo_file = file('logo.txt')
if (logo_file.exists()) {
    custom_logo = logo_file.text
    println(custom_logo)
} else {
    custom_logo = "Where is the face of metaFun?? :l"
    println(custom_logo)
}
println("\nVisit documentation to find out description of this workflow: \n")
println("https://metafun-doc.readthedocs.io/en/latest/workflows/COMPARATIVE_ANNOTATION.html\n")

def checkProcessOutput(String outputPath) {
    return file(outputPath).exists()
}


// def determineRunMode() {
//     def annotationExists = file("${params.annotation_results}").exists()
//     def metadataExists = params.metadata != null
//     def metacolExists = params.metacol != null

//     if (annotationExists && metadataExists && metacolExists) {
//         return 'visualize'
//     } else if (metadataExists && metacolExists) {
//         return 'full'
//     } else {
//         return 'annotate'
//     }
// }
//params.run_mode = determineRunMode() //'full', 'annotate', 'visualize'

// color declaration 
def ANSI_RESET = "\u001B[0m"
def ANSI_BLUE = "\033[38;2;78;149;217m"
def ANSI_YELLOW = "\u001B[33m"
def ANSI_GREEN = "\u001B[32m"
def ANSI_CYAN = "\033[38;2;0;255;255m"
def ANSI_RED = "\u001B[31m"
def ANSI_ORANGE = "\033[38;2;255;149;0m"


params.help = null
if (params.help) {
    log.info """
Usage:
    nextflow run nf_scripts/${ANSI_BLUE}COMPARATIVE_ANNOTATION_apptainer.nf${ANSI_RESET}  --metadata "\${your_metadata}" --metacol "\${your_metacol}"
${ANSI_YELLOW}COMPARATIVE_ANNOTATION : Comparative genomic analysis workflow${ANSI_RESET}
=========================================================
 ${ANSI_GREEN}PPanGGOLiN${ANSI_RESET}        : Identity = ${params.pan_identity}, Coverage = ${params.pan_coverage}
 ${ANSI_GREEN}KEGG${ANSI_RESET}               : Module presence threshold for presence = ${params.module_completeness}, E-value = ${params.kofamscan_eval}
 ${ANSI_GREEN}VFDB${ANSI_RESET}               : Identity = ${params.VFDB_identity}, Coverage = ${params.VFDB_coverage}, E-value = ${params.VFDB_e_value}
 ${ANSI_GREEN}CARD${ANSI_RESET}               : Using RGI
 ${ANSI_GREEN}dbCAN${ANSI_RESET}              : Using run_dbcan4
 ${ANSI_GREEN}eggNOG${ANSI_RESET}             : Using eggnog-mapper


CPUs used          : ${ANSI_BLUE}${params.cpus}${ANSI_RESET}
Metadata           : ${ANSI_BLUE}${params.metadata}${ANSI_RESET}
Metadata column    : ${ANSI_BLUE}${params.metacol}${ANSI_RESET}
Input directory    : ${ANSI_BLUE}${params.inputDir}${ANSI_RESET}
Output directory   : ${ANSI_BLUE}${params.outdir}${ANSI_RESET}

Additional options:
--pan_identity     : PPanGGOLiN identity threshold (default: ${ANSI_BLUE}${params.pan_identity}${ANSI_RESET})
--pan_coverage     : PPanGGOLiN coverage threshold (default: ${ANSI_BLUE}${params.pan_coverage}${ANSI_RESET})
--module_completeness : KEGG module completeness threshold (default: ${ANSI_BLUE}${params.module_completeness}${ANSI_RESET})
--kofamscan_eval : KEGG KO e-value threshold (default: ${ANSI_BLUE}${params.kofamscan_eval}${ANSI_RESET})
--kingdom          : Kingdom for PPPanGGOLiN (default: ${ANSI_BLUE}${params.kingdom}${ANSI_RESET}). Choose bacteria or archaea
--VFDB_identity : Virulence factor identity threshold (default: ${ANSI_BLUE}${params.VFDB_identity}${ANSI_RESET} )
--VFDB_coverage : Virulence factor identity threshold (default: ${ANSI_BLUE}${params.VFDB_coverage}${ANSI_RESET})
--VFDB_e_value : Virulence factor e_value threshold (default: ${ANSI_BLUE}${params.VFDB_e_value}${ANSI_RESET})
--CAZyme_hmm_eval  : CAZyme hmm e_value threshold (default: ${ANSI_BLUE}${params.CAZyme_hmm_eval}${ANSI_RESET})
--CAZyme_hmm_cov  : CAZyme hmm coverage threshold (default: ${ANSI_BLUE}${params.CAZyme_hmm_cov}${ANSI_RESET})

"""
    exit 0
}

if (params.metacol == null || params.metacol == "" || params.metacol.toInteger() == 0) {
    log.warn "metacol is not specified. Visualization step will be skipped. Only other processes will be executed."
    params.run_visualization = false
} else {
    params.run_visualization = true
}


// if (params.metacol == null || params.metacol == 0) {
//     log.error 'Metadata column parameter --metacol cannot be empty.\n' + 
//                'Metadata column should be selected. '
//                'You can specify by --metacol {column number} ' + 
//                'You can check {column number} by using  $head -n1 selected_metadata.csv | tr "," "\n" | nl'
//     System.exit(1) // Exits the script with an error status
// }

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
// help log of comparative annotation . 



log.info """\

===========================================
${ANSI_BLUE}COMPARATIVE_ANNOTATION${ANSI_RESET} - Functional Comparison of MAGs
===========================================
(Apptainer Container with Nextflow)

Key Features:
  ✓ KEGG Orthology annotation (KofamScan)
  ✓ Orthology assignment and functional annotation (eggNOG)
  ✓ Virulence factor annotation (VFDB)
  ✓ Antimicrobial resistance gene annotation (CARD/RGI)
  ✓ Carbohydrate-Active enzyme annotation (dbCAN)
  ✓ Pangenome analysis (PPanGGOLiN)
  ✓ Genome similarity and subspecies dereplication (skani/PCoA/dRep)


${ANSI_RED}Required Parameters:${ANSI_RESET}
------------------
${ANSI_RED}--inputDir${ANSI_RESET}         : Directory containing bins from ASSEMBLY_BINNING or  your own genomes
                     Typically output from ASSEMBLY_BINNING (final_bins)
${ANSI_RED}--metadata${ANSI_RESET}         : CSV file with additional metadata for samples
                     Typically output from GENOME_SELECTOR (genome_selector_result.csv)
${ANSI_RED}--metacol${ANSI_RESET}          : Column number in metadata file to use for grouping


Optional Parameters:
------------------
${ANSI_YELLOW}--module_completeness${ANSI_RESET} : KEGG module completeness threshold (default: ${params.module_completeness})
${ANSI_YELLOW}--pan_identity${ANSI_RESET}     : Identity threshold for PPanGGOLiN (default: ${params.pan_identity})
${ANSI_YELLOW}--pan_coverage${ANSI_RESET}     : Coverage threshold for PPanGGOLiN (default: ${params.pan_coverage})
${ANSI_YELLOW}--kofamscan_eval${ANSI_RESET}   : E-value threshold for KofamScan (default: ${params.kofamscan_eval})
${ANSI_YELLOW}--run_drep${ANSI_RESET}         : Run dereplication (default: ${params.run_drep})
${ANSI_YELLOW}--drep_ani${ANSI_RESET}         : ANI threshold for dereplication (default: ${params.drep_ani})
${ANSI_YELLOW}--drep_algorithm${ANSI_RESET}   : Algorithm for ANI calculation (default: ${params.drep_algorithm})
${ANSI_YELLOW}--kingdom${ANSI_RESET}          : Kingdom for annotation (default: ${params.kingdom})
${ANSI_YELLOW}--cpus${ANSI_RESET}             : Number of CPUs to use (default: ${params.cpus})


${ANSI_CYAN}Database Information:${ANSI_RESET}
------------------
${ANSI_CYAN}KEGG${ANSI_RESET}               : KofamScan profiles and KEGG Orthology database
${ANSI_CYAN}eggNOG${ANSI_RESET}             : eggNOG v5.0 orthologous groups and annotations
${ANSI_CYAN}VFDB${ANSI_RESET}               : Virulence Factor Database (${params.VFDB_identity}% identity, ${params.VFDB_coverage}% coverage)
${ANSI_CYAN}CARD${ANSI_RESET}               : Comprehensive Antibiotic Resistance Database
${ANSI_CYAN}dbCAN${ANSI_RESET}              : CAZyme annotation (${params.CAZyme_hmm_eval} e-value, ${params.CAZyme_hmm_cov} coverage)
${ANSI_CYAN}PPanGGOLiN${ANSI_RESET}         : Pangenome partitioning for ${params.kingdom}


${ANSI_CYAN}Output Directories:${ANSI_RESET}
------------------
${ANSI_CYAN}- annotation_results/${ANSI_RESET}   : Individual genome annotations
${ANSI_CYAN}- visualization_results/${ANSI_RESET} : Comparative analysis and visualizations

${params.run_drep ? ANSI_YELLOW + 'dRep will be run for genome dereplication with ' + params.drep_algorithm + ' algorithm' + ANSI_RESET : ANSI_YELLOW + 'dRep will be skipped' + ANSI_RESET}

\033[38;2;180;180;180mNote:\033[0m Results can be visualized with: 
${ANSI_ORANGE}INTERACTIVE_COMPARATIVE${ANSI_RESET}
"""
.stripIndent(true)

    // Pangenome analysis, species functional annotation 
    // =================================================
    // PPPanGGOLiN option          : Identity = ${params.pan_identity}, 
    //                             Coverage = ${params.pan_coverage}
    //                             modify it if you want by specifying --pan_identity [0~1] --pan_coverage [0~1] in command line 
    // KEGG module present threshold : ${params.module_completeness}
    //                             modify threshold if you want to change it : --module_completeness [0~1]
    
    // """
    // .stripIndent(true)


process prepare_genomes {
    publishDir "${params.outdir}/selected_genomes", mode: 'copy'

    input:
    path metadata_file
    path input_dir

    output:
    path "selected_genomes${params.run_name}"

    script:
    """
    mkdir -p selected_genomes${params.run_name}
    awk -F',' 'NR>1 {print \$${params.samplecol}}' ${metadata_file} | while read genome; do
        if [ -f "${input_dir}/\${genome}.fa" ]; then
            cp "${input_dir}/\${genome}.fa" selected_genomes${params.run_name}/
        else
            echo "Warning: \${genome}.fa not found in ${input_dir}" >&2
        fi
    done
    """
}



// create metadata header information file.
process create_metadata_summary {
    publishDir "${params.outdir}/annotation_results", mode: 'copy'
    //publishDir "${launchDir}", mode: 'copy'
    input:
    path metadata_file
    output:
    path "metadata_column_COMPARATIVE_ANNOTATION_summary.tsv"

    // when:
    // params.metacol != "" && params.metacol.toInteger() > 0

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv("${metadata_file}", sep=',', header=0)

    summary = pd.DataFrame({
        'Column_Index': range(len(df.columns)),
        'Column_Name': df.columns
    })

    summary.to_csv('metadata_column_COMPARATIVE_ANNOTATION_summary.tsv', sep='\t', index=False)
    """
}

process run_skani_annotation {
    publishDir "${params.annotation_results}/ani", mode: 'copy'
    cpus params.cpus

    input:
    path(genome_dir)

    output:
    path "skani_fullmatrix"
    path "skani_ANI_dist.tsv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    skani triangle *.fa -t ${task.cpus} --full-matrix | tail -n +2 | sed -e 's/.fa//g' > skani_fullmatrix

    skani triangle *.fa -t ${task.cpus} --full-matrix --distance > test
    #skani triangle *.fa -t ${task.cpus} --full-matrix -l rl --distance > test
    cut -f1 test | tail -n+2 | tr '\n' '\t' | sed -e 's/^/\t/g' | sed -e 's/\t\$/\\n/' > skani_ANI_dist.tsv
    tail -n+2 test >> skani_ANI_dist.tsv
    sed -i 's/\\.fa//g' skani_ANI_dist.tsv
    rm test


    """
}
process run_skani_visualization {
    publishDir "${params.visualization_results}/ani/column_${params.metacol}", mode: 'copy'
    when: params.run_visualization

    input:
    path(skani_fullmatrix)
    path(metadata_file)

    output:
    path "heatmap_skani.pdf"
    path "skani_interactive"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/skani_visualization.R \
    -i ${skani_fullmatrix} -m ${metadata_file} -mc ${params.metacol} -out heatmap_skani.pdf -out_html skani_interactive
    """
}

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
    //container = 'docker://staphb/prokka:latest'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    //cpus params.cpus

    input:
    tuple val(id), path(genome)
    output:
    path "*"

    script:
    """
    #eval "\$(micromamba shell hook --shell bash)"
    #micromamba  activate prokka
    #export TMPDIR='/tmp/\${uuidgen()}'
    #export JAVA_OPTS="-XX:-UsePerfData -Djava.io.tmpdir=\$TMPDIR"
    JAVA_TOOL_OPTIONS="-Xlog:disable -Xlog:all=warning:stderr:uptime,level,tags"

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
    --kingdom  ${params.kingdom} --rarefaction
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
    python ${params.scripts_baseDir}/genome_PA_stat.py --cpus ${task.cpus} --metadata ${params.metadata}

    """

}


process run_kofamscan_annotation {
    publishDir "${params.annotation_results}/kofamscan", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("ko_matrix.csv")
    path("KO_definition_GeneID_countgenomes.csv")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    
    echo running kofamscan
    time exec_annotation \
    -o kofamscan_result --cpu ${task.cpus} --e-value ${params.kofamscan_eval}\
    -c ${params.db_baseDir}/kofam/27Nov2023/config-template.yml \
    -f mapper ppanggolin_result/pan_genes/all_protein_families.faa
    
    python ${params.scripts_baseDir}/kofam_to_geneID_KO_ppanggolin205.py \
    -i ppanggolin_result/matrix.csv -o output_gene_pa_ko.csv
#KO_definition_GeneID_countgenomes.py generate KO_definition_GeneID_countgenomes.csv 
    python ${params.scripts_baseDir}/KO_definition_GeneID_countgenomes.py -i  output_gene_pa_ko.csv \
    -o KO_definition_GeneID_countgenomes.csv -k ${params.db_baseDir}/kofam/27Nov2023/ko_list

    python ${params.scripts_baseDir}/KO_Genome_long2matrix_ppanggolin205.py \
    -i output_gene_pa_ko.csv -o ko_matrix.csv
    """
}

process run_kofamscan_visualization {
    publishDir "${params.visualization_results}/kofamscan", mode: 'copy'
    publishDir "${params.visualization_results}/kofamscan/column_${params.metacol}", mode: 'copy'
    when: params.run_visualization
    
    input:
    path(ko_matrix)
    path(metadata_file)
    
    output:
    path("KEGG_module_visualization_shiny")
    path("KEGG_module_completeness.csv")
    path("heatmap_KEGG.pdf")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/KO_module_visualization_apptainer.R \
    -i ${ko_matrix} -m ${metadata_file} -n ${params.metacol} -mc ${params.module_completeness} \
    -out_table KEGG_module_completeness.csv -out_html KEGG_module_visualization_shiny
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

// process run_VFDB_annotation_tmp {
//     publishDir "${params.annotation_results}/VFDB", mode: 'copy'
//     cpus params.cpus

//     input:
//     path(ppanggolin_dir)

//     output:
//     path("pangene_vfdb_result.txt"), emit: vfdb_result
//     path("gene_PA_VFDB_added.csv"), emit: gene_pa    
    

//     script:
//     """
//     eval "\$(micromamba shell hook --shell bash)"
//     micromamba activate COMPARATIVE_ANNOTATION
//     micromamba activate rgi603
//     diamond blastp -d ${params.db_baseDir}/VFDB/VFDB_setB_prot_Aug2023.dmnd \
//     -q ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_vfdb_result.txt -f 6 --max-target-seqs 1 --id ${params.VFDB_identity} \
//     --subject-cover ${params.VFDB_coverage} -e ${params.VFDB_e_value} 

//     python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
//     -d pangene_vfdb_result.txt -g ppanggolin_result/gene_presence_absence.Rtab \
//     -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
//     -o gene_PA_VFDB_added.csv
//     """
// }
// process run_VFDB_visualization_tmp {
//     publishDir "${params.visualization_results}/VFDB", mode: 'copy'

//     input:
//     path(gene_PA_VFDB_added)
//     path(metadata_file)

//     output:
//     path("heatmap_VFDB.pdf")
//     path("VFDB_interactive")

//     script:
//     """
//     eval "\$(micromamba shell hook --shell bash)"
//     micromamba activate R432_environment
//     Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
//     -i ${gene_PA_VFDB_added} -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive
//     """
// }

//process run_VFDB_count_annotation {
process run_VFDB_annotation {
    publishDir "${params.annotation_results}/VFDB", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("pangene_vfdb_result.txt"), emit: vfdb_result
    path("gene_PA_VFDB_added.csv"), emit: vfdb_gene_pa    
    path("gene_count_VFDB_added.csv"), emit: vfdb_gene_count

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    micromamba activate rgi603
    diamond blastp -d ${params.db_baseDir}/VFDB/VFDB_setB_prot_Aug2023.dmnd \
    -q ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_vfdb_result.txt -f 6 --max-target-seqs 1 --id ${params.VFDB_identity} \
    --subject-cover ${params.VFDB_coverage} -e ${params.VFDB_e_value} 


    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt -g ppanggolin_result/gene_presence_absence.Rtab \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_PA_VFDB_added.csv


    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt \
    -g ppanggolin_result/gene_count_matrix.tsv \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_count_VFDB_added.csv     
    """
}
//process run_VFDB_count_visualization {

process run_VFDB_visualization {
    publishDir "${params.visualization_results}/VFDB", mode: 'copy'
    //publishDir "${params.visualization_results}/VFDB/${params.gene_matrix}_column_${params.metacol}", mode: 'copy'
    when: params.run_visualization

    input:
    path(gene_PA_VFDB_added)
    path(gene_count_VFDB_added)
    path(metadata_file)

    output:
    //path("heatmap_VFDB.pdf")
    //path("VFDB_interactive")
    path("heatmap_VFDB_gene_PA_metadata${params.metacol}th_col.pdf")
    path("VFDB_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_VFDB_gene_count_metadata${params.metacol}th_col.pdf")
    path("VFDB_interactive_gene_count_metadata${params.metacol}th_col")

    script:
    //def input_file = params.gene_matrix == 'gene_PA' ? gene_PA_VFDB_added : gene_count_VFDB_added
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
        # Run visualization for gene_PA
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i ${gene_PA_VFDB_added} -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_VFDB.pdf heatmap_VFDB_gene_PA_metadata${params.metacol}th_col.pdf

    # Run visualization for gene_count
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i ${gene_count_VFDB_added} -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_VFDB.pdf heatmap_VFDB_gene_count_metadata${params.metacol}th_col.pdf
    """
}
 // Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
  //     -i ${gene_count_VFDB_added} -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive

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

process run_rgi_CARD_annotation {
    publishDir "${params.annotation_results}/CARD", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("pangene_rgi_CARD_result.txt"), emit: rgi_result
    path("gene_PA_CARD_added.csv"), emit: card_gene_pa
    path("gene_count_CARD_added.csv"), emit: card_gene_count
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    micromamba activate rgi603
    rgi load -i ${params.db_baseDir}/CARD/card.json --local 
    rgi main -i ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_rgi_CARD_result --clean  \
    -t protein -n ${task.cpus} --include_nudge --local 

    python ${params.scripts_baseDir}/add_rgi_togenePA.py -i pangene_rgi_CARD_result.txt \
    -o gene_PA_CARD_added.csv -r ${params.db_baseDir}/CARD/aro_index.tsv \
    -gpa ppanggolin_result/gene_presence_absence.Rtab

    python ${params.scripts_baseDir}/add_rgi_togenePA.py -i pangene_rgi_CARD_result.txt \
    -o gene_count_CARD_added.csv -r ${params.db_baseDir}/CARD/aro_index.tsv \
    -gpa ppanggolin_result/gene_count_matrix.tsv    
    """
}

process run_rgi_CARD_visualization {
    publishDir "${params.visualization_results}/CARD", mode: 'copy'
    when: params.run_visualization
    input:
    path(gene_PA_CARD_added)
    path(gene_count_CARD_added)
    path(metadata_file)

    output:
    path("heatmap_CARD_gene_PA_metadata${params.metacol}th_col.pdf")
    path("CARD_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_CARD_gene_count_metadata${params.metacol}th_col.pdf")
    path("CARD_interactive_gene_count_metadata${params.metacol}th_col")


    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    
    # Run visualization for gene_PA
    Rscript ${params.scripts_baseDir}/CARD_visualization.R \
    -i ${gene_PA_CARD_added} -m ${metadata_file} -mc ${params.metacol} -out_html CARD_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_CARD.pdf heatmap_CARD_gene_PA_metadata${params.metacol}th_col.pdf

    # Run visualization for gene_count
    Rscript ${params.scripts_baseDir}/CARD_visualization.R \
    -i ${gene_count_CARD_added} -m ${metadata_file} -mc ${params.metacol} -out_html CARD_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_CARD.pdf heatmap_CARD_gene_count_metadata${params.metacol}th_col.pdf
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

process run_defensefinder_annotation {
    publishDir "${params.annotation_results}/defensefinder", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("all_protein_families_defense_finder_genes.tsv"), emit: defensefinder_result
    path("gene_PA_defense.csv"), emit: defense_gene_pa
    path("gene_count_defense.csv"), emit: defense_gene_count

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate defensefinder

    defense-finder run -w ${task.cpus} --db-type unordered --models-dir ${params.db_baseDir}/defensefinder/models --antidefensefinder -o defensefinder_out --db-type unordered ppanggolin_result/pan_genes/all_protein_families.faa

    mv defensefinder_out/all_protein_families_defense_finder_genes.tsv .
    micromamba activate COMPARATIVE_ANNOTATION
    python ${params.scripts_baseDir}/gene_matrix_add_defense.py -d all_protein_families_defense_finder_genes.tsv -g ppanggolin_result/gene_presence_absence.Rtab -o gene_PA_defense.csv

    python ${params.scripts_baseDir}/gene_matrix_add_defense.py -d all_protein_families_defense_finder_genes.tsv -g ppanggolin_result/gene_count_matrix.tsv -o gene_count_defense.csv
    """
}

process run_defensefinder_visualization {
    publishDir "${params.visualization_results}/defensefinder", mode: 'copy'
    when: params.run_visualization

    input:
    path(gene_PA_defense)
    path(gene_count_defense)
    path(metadata_file)

    output:
    path("heatmap_defensefinder_gene_PA_metadata${params.metacol}th_col.pdf")
    path("defensefinder_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_defensefinder_gene_count_metadata${params.metacol}th_col.pdf")
    path("defensefinder_interactive_gene_count_metadata${params.metacol}th_col")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    
    # Run visualization for gene_PA
    Rscript ${params.scripts_baseDir}/Defense_visualization.R \
    -i ${gene_PA_defense} -m ${metadata_file} -mc ${params.metacol} -out_html defensefinder_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_defensefinder.pdf heatmap_defensefinder_gene_PA_metadata${params.metacol}th_col.pdf

    # Run visualization for gene_count
    Rscript ${params.scripts_baseDir}/Defense_visualization.R \
    -i ${gene_count_defense} -m ${metadata_file} -mc ${params.metacol} -out_html defensefinder_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_defensefinder.pdf heatmap_defensefinder_gene_count_metadata${params.metacol}th_col.pdf
    """
}









process run_dbCAN_annotation_bak {
    publishDir "${params.annotation_results}/dbCAN", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("dbcan_hmmerfamily_count_matrix.csv")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate run_dbcan4
    export OMP_NUM_THREADS=${task.cpus}

    time run_dbcan ppanggolin_result/pan_genes/all_protein_families.faa protein --dia_cpu ${task.cpus} \
    --hmm_cpu ${task.cpus} --out_dir db_can_out --db_dir ${params.db_baseDir}/dbCAN

    cut -f1,3 db_can_out/overview.txt | grep -Pv "\\t-\$" | awk -F'\\t' '{split(\$2,a,"+"); for (i in a) {print \$1 "\\t" a[i];}}' | \
    awk -F"\\t" '{split(\$2,a,"_"); \$2=a[1]; print \$1 "\\t" \$2}' | sed -e "s/(.*//" > dbcan_family.tsv

    python ${params.scripts_baseDir}/dbcan_result_parse.py
    """
}


process run_dbCAN_annotation {
    publishDir "${params.annotation_results}/dbCAN", mode: 'copy'
    cpus params.cpus  // 

    input:
    path(ppanggolin_dir)

    output:
    path("db_can_out"), emit: db_can_out
    path("dbcan_raw_gene_count_data.csv"), emit: raw_gene_count_data
    path("dbcan_geneID_HMMER_count_gene_count.csv"), emit: geneID_HMMER_count_gene_count
    path("dbcan_HMMER_count_gene_count_matrix.csv"), emit: HMMER_count_gene_count_matrix
    path("dbcan_raw_gene_PA_data.csv"), emit: raw_gene_PA_data
    path("dbcan_geneID_HMMER_count_gene_PA.csv"), emit: geneID_HMMER_count_gene_PA
    path("dbcan_HMMER_count_gene_PA_matrix.csv"), emit: HMMER_count_gene_PA_matrix

    script:
    def parallel_jobs = (params.cpus as int).intdiv(4)  // divide maximum cpu by 4 
    """
    # Activate environment for file splitting
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC

    # Split the input file
    seqkit split2 ppanggolin_result/pan_genes/all_protein_families.faa -p ${parallel_jobs} -O ./

    # Switch environment for dbCAN analysis
    micromamba activate run_dbcan4

    # Run dbCAN on each part in parallel
    ls all_protein_families.part*.faa | /opt/conda/envs/COMPARATIVE_ANNOTATION/bin/parallel -j ${parallel_jobs} run_dbcan {} protein --dia_cpu 4 --hmm_cpu 4 --out_dir {}_dir --db_dir ${params.db_baseDir}/dbCAN \
    --hmm_eval  ${params.CAZyme_hmm_eval} --hmm_cov ${params.CAZyme_hmm_cov}

    # Merge the results
    mkdir db_can_out

    for file in overview.txt dbcan-sub.hmm.out diamond.out hmmer.out; do
        # Get the header from the first file
        head -n 1 \$(ls all_protein_families.part*_dir/\$file | head -n 1) > db_can_out/\$file
        # Append all data rows (excluding headers) from all files
        for f in all_protein_families.part*_dir/\$file; do
            tail -n +2 \$f >> db_can_out/\$file
        done
    done
    
    cat all_protein_families.part*_dir/uniInput > db_can_out/uniInput
    # cat all_protein_families.part*_dir/overview.txt > db_can_out/overview.txt
    # cat all_protein_families.part*_dir/dbcan-sub.hmm.out > db_can_out/dbcan-sub.hmm.out
    # cat all_protein_families.part*_dir/diamond.out > db_can_out/diamond.out
    # cat all_protein_families.part*_dir/hmmer.out > db_can_out/hmmer.out


    # Process the merged results
    cut -f1,3 db_can_out/overview.txt | grep -Pv "\\t-\$" | awk -F'\\t' '{split(\$2,a,"+"); for (i in a) {print \$1 "\\t" a[i];}}' | \
    awk -F"\\t" '{split(\$2,a,"_"); \$2=a[1]; print \$1 "\\t" \$2}' | sed -e "s/(.*//" > dbcan_family.tsv

    python ${params.scripts_baseDir}/dbcan_result_parse.py

    # Clean up intermediate files
    rm -rf all_protein_families.part*.faa all_protein_families.part*_dir  
    """
}



process run_dbCAN_visualization {
   // publishDir "${params.visualization_results}/dbCAN/${params.gene_matrix}_column_${params.metacol}", mode: 'copy'
    publishDir "${params.visualization_results}/dbCAN", mode: 'copy'
    when: params.run_visualization
    
    input:
    path(dbcan_HMMER_count_gene_count_matrix)
    path(dbcan_HMMER_count_gene_PA_matrix)
    path(metadata_file)

    output:
    path("heatmap_dbCAN_gene_PA_metadata${params.metacol}th_col.pdf")
    path("dbCAN_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_dbCAN_gene_count_metadata${params.metacol}th_col.pdf")
    path("dbCAN_interactive_gene_count_metadata${params.metacol}th_col")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/dbcan_visualization.R -i ${dbcan_HMMER_count_gene_PA_matrix} \
    -m ${metadata_file} -mc ${params.metacol} -out_html dbCAN_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_dbCAN.pdf heatmap_dbCAN_gene_PA_metadata${params.metacol}th_col.pdf

    Rscript ${params.scripts_baseDir}/dbcan_visualization.R -i ${dbcan_HMMER_count_gene_count_matrix} \
    -m ${metadata_file} -mc ${params.metacol} -out_html dbCAN_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_dbCAN.pdf heatmap_dbCAN_gene_count_metadata${params.metacol}th_col.pdf
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
    publishDir "${params.annotation_results}/eggNOG" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(ppanggolin_dir)

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
    publishDir "${params.visualization_results}/scoary2" , mode : 'copy'
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


process run_drep_dereplication {
    publishDir "${params.outdir}/drep", mode: 'copy'
    cpus params.cpus

    input:
    path(genome_dir)

    output:
    path "drep_output"
    path "dereplicated_genomes"
    path "subspecies_clusters.tsv"


    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    #micromamba activate drep350

      #-g ${genome_dir}/*.fa \

    # Run dRep dereplication
    dRep dereplicate drep_output \
      -g *.fa \
      -p ${task.cpus} \
      --ignoreGenomeQuality \
      --S_algorithm ${params.drep_algorithm} \
      -sa ${params.drep_ani} \
      -nc ${params.drep_cov}
      
    # Copy dereplicated genomes to a separate directory for easier access
    mkdir -p dereplicated_genomes
    cp drep_output/dereplicated_genomes/*.fa dereplicated_genomes/

    echo -e "genome\tsubspecies_cluster\tprimary_cluster" > subspecies_clusters.tsv
    awk -F, 'NR>1 {print \$1"\\t"\$2"\\t"\$6}' drep_output/data_tables/Cdb.csv >> subspecies_clusters.tsv

    """
}



/*
process run_kofamscan 
*/
/*
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
*/

process create_shiny_dashboard {
    publishDir "${params.outdir}/shiny_dashboard", mode: 'copy'

    input:
    path visualization_results_dir

    output:
    path "shiny_dashboard.html"

    script:
    """
    #!/usr/bin/env python3
    import os

    html_content = "<html><body>"
    html_content += "<h1>Shiny App Links</h1>"

    for root, dirs, files in os.walk("${visualization_results_dir}"):
        for file in files:
            if file == "htShiny.sh":
                app_path = os.path.join(root, file)
                app_name = os.path.basename(os.path.dirname(app_path))
                relative_path = os.path.relpath(app_path, "${params.outdir}")
                html_content += f'<p><a href="{relative_path}" target="_blank">{app_name} Shiny App</a></p>'

    html_content += "</body></html>"

    with open("shiny_dashboard.html", "w") as f:
        f.write(html_content)
    """
}
process run_multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus 1

    input:
    path visualization_results
    path shiny_dashboard

    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    export MPLCONFIGDIR=~/.config/matplotlib
    mkdir -p ~/.config/matplotlib


    multiqc ${visualization_results}  -o ${params.outdir}/multiqc
    cp ${shiny_dashboard} ${params.outdir}/multiqc/
    echo '<p>See <a href="shiny_dashboard.html">Shiny Apps Dashboard</a></p>' >> ${params.outdir}/multiqc/multiqc_report.html
    """
}

process make_sequence_db {
    publishDir "${params.outdir}/sequence_db", mode: 'copy'
    cpus params.cpus

    input:
    path metadata_file
    path prokka_dir
    path ppanggolin_dir

    output:
    path "sequences.h5"
    path "logs"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    #micromamba activate COMPARATIVE_ANNOTATION

    mkdir -p logs
    mkdir -p prokka
    mkdir -p ppanggolin_need 
        find -L ./ -name "*.faa" | grep -v ppanggolin_result | xargs -I {} cp {} prokka/
    find -L ./ -name "*.ffn" | grep -v ppanggolin_result | xargs -I {} cp {} prokka/
    
    cp ${ppanggolin_dir}/gene_families.tsv ppanggolin_need/
    cp ${ppanggolin_dir}/genome_id_linking.tsv ppanggolin_need/

    
    Rscript /scratch/tools/microbiome_analysis/comparative_annotation/make_sequence_db.R \
      --input ${metadata_file} \
      --processors ${task.cpus} \
      --output ./ \
      --data ./
    
 #rm -rf prokka
 #rm -rf ppanggolin_need
    """
}


workflow {
    log.info "Starting analysis : your output directory is at  ${params.outdir}"
    //log.info "Determined run mode: ${params.run_mode}"

    metadata_ch = Channel.fromPath(params.metadata)
    input_dir_ch = Channel.fromPath(params.inputDir)

    selected_genomes = prepare_genomes(metadata_ch, input_dir_ch)

    def runProkka = !file("${params.outdir}/prokka").exists()
    def runPpanggolin = !file("${params.outdir}/ppanggolin_result").exists()
    def runKofamscan = !file("${params.annotation_results}/kofamscan/ko_matrix.csv").exists()
    def runVFDB = !file("${params.annotation_results}/VFDB/gene_count_VFDB_added.csv").exists()
    def runCARD = !file("${params.annotation_results}/CARD/gene_count_CARD_added.csv").exists()
    def runDbCAN = !file("${params.annotation_results}/dbCAN/dbcan_hmmerfamily_count_matrix.csv").exists()
    def runSkani = !file("${params.annotation_results}/ani/skani_fullmatrix").exists()
    def runEggNOG = !file("${params.annotation_results}/eggNOG/eggnog_mmseqs.emapper.annotations").exists()
    def runMetadataSummary = !file("${params.annotation_results}/metadata_column_COMPARATIVE_ANNOTATION_summary.tsv").exists()
    def runGenePACluster = !file("${params.outdir}/genePA_cluster/pcoa_plot_interactive.html").exists()
    def run_sequence_db = !file("${params.outdir}/sequence_db/sequences.h5").exists()
    //def runDefenseFinder = !file("${params.outdir}/defensefinder/all_protein_families_defense_finder_genes.tsv").exists()

    def runDrep = params.run_drep && !file("${params.outdir}/drep/subspecies_clusters.tsv").exists()

    // if (params.run_mode == 'full') {
    //     // Full mode: Run both annotation and visualization
    //     bins_ch = selected_genomes.flatMap { dir ->
    //         file("${dir}/*.fa").collect { fa_file ->
    //             tuple(fa_file.baseName, fa_file)
    //         }

    bins_ch = selected_genomes.map { dir ->
        // 디렉토리 내 모든 .fa 파일 찾기
        def fa_files = file("${dir}/*.fa")
        if (fa_files.size() == 0) {
            log.warn "No .fa files found in ${dir}"
        }
        return fa_files
    }.flatten().map { fa_file ->
        // 각 파일을 tuple로 변환
        return tuple(fa_file.baseName, fa_file)
    }
    
    // 확인을 위해 파일 목록 로깅
    bins_ch.view { id, file -> "Processing genome: ${id} (${file})" }
    prokka_result = run_prokka(bins_ch)
    prokka_collect = prokka_result.collect()

    // if (runProkka) {
    //     bins_ch = selected_genomes.flatMap { dir ->
    //         file("${dir}/*.fa").collect { fa_file ->
    //             tuple(fa_file.baseName, fa_file)
    //         }
    //     }
    //     prokka_result = run_prokka(bins_ch)
    //     prokka_collect = prokka_result.collect()
    // } else {
    //     prokka_collect = Channel.fromPath("${params.outdir}/prokka")
    // }

    if (runPpanggolin) {
        ppanggolin_result = run_ppanggolin(prokka_collect)
    } else {
        ppanggolin_result = Channel.fromPath("${params.outdir}/ppanggolin_result")
    }
    if (runKofamscan) {
        kofamscan_result = run_kofamscan_annotation(ppanggolin_result)
    } else {
        kofamscan_result = Channel.fromPath("${params.annotation_results}/kofamscan/ko_matrix.csv")
    }
    if (runVFDB) {
        vfdb_result = run_VFDB_annotation(ppanggolin_result)
        vfdb_gene_pa = vfdb_result.vfdb_gene_pa
        vfdb_gene_count = vfdb_result.vfdb_gene_count
    } else {
        vfdb_gene_pa = Channel.fromPath("${params.annotation_results}/VFDB/gene_PA_VFDB_added.csv")
        vfdb_gene_count = Channel.fromPath("${params.annotation_results}/VFDB/gene_count_VFDB_added.csv")
    }

    if (runCARD) {
        card_result = run_rgi_CARD_annotation(ppanggolin_result)
        card_gene_pa = card_result.card_gene_pa
        card_gene_count = card_result.card_gene_count
    } else {
        card_gene_pa = Channel.fromPath("${params.annotation_results}/CARD/gene_PA_CARD_added.csv")
        card_gene_count = Channel.fromPath("${params.annotation_results}/CARD/gene_count_CARD_added.csv")
    }
    // if (runDefenseFinder) {
    //     defensefinder_result = run_defensefinder_annotation(ppanggolin_result)
    //     defensefinder_gene_pa = defensefinder_result.defense_gene_pa
    //     defensefinder_gene_count = defensefinder_result.defense_gene_count
    // } else {
    //     defensefinder_gene_pa = Channel.fromPath("${params.annotation_results}/defensefinder/gene_PA_defense.csv")
    //     defensefinder_gene_count = Channel.fromPath("${params.annotation_results}/defensefinder/gene_count_defense.csv")
    // }

    if (runDbCAN) {
        dbcan_result = run_dbCAN_annotation(ppanggolin_result)
        dbcan_count_matrix = dbcan_result.HMMER_count_gene_count_matrix
        dbcan_pa_matrix = dbcan_result.HMMER_count_gene_PA_matrix
    } else {
        dbcan_count_matrix = Channel.fromPath("${params.annotation_results}/dbCAN/dbcan_HMMER_count_gene_count_matrix.csv")
        dbcan_pa_matrix = Channel.fromPath("${params.annotation_results}/dbCAN/dbcan_HMMER_count_gene_PA_matrix.csv")
    }
    if (runSkani) {
        all_genomes = selected_genomes.flatMap { dir -> file("${dir}/*.fa") }.collect()
        skani_result = run_skani_annotation(all_genomes)
        skani_fullmatrix = skani_result[0]
        skani_dist_result = skani_result[1]
    } else {
        skani_fullmatrix = Channel.fromPath("${params.annotation_results}/ani/skani_fullmatrix")
        skani_dist_result = Channel.fromPath("${params.annotation_results}/ani/skani_ANI_dist.tsv")
    }

    if (runEggNOG) {
        eggnog_result = run_eggNOG(ppanggolin_result)
    } else {
        log.info "eggNOG already exists."

//        eggnog_result = Channel.fromPath("${params.annotation_results}/eggNOG/eggnog_mmseqs.emapper.annotations")
    }
    if (runMetadataSummary) {
        create_metadata_summary(metadata_ch)
    }
    if (run_sequence_db) {
        sequence_db = make_sequence_db(metadata_ch, prokka_collect, ppanggolin_result)
    } else {
        sequence_db = Channel.fromPath("${params.outdir}/sequence_db/sequences.h5")
    }
    
    def genomeCount = new File("selected_genomes${params.run_name}").exists() ? new File("selected_genomes${params.run_name}").listFiles().size() : 0

    if (runGenePACluster) {
        if (genomeCount >= 5) {
            run_genePA_cluster(ppanggolin_result, metadata_ch)
        } else {
            log.info "Genome count is less than 5. GenePA Cluster will not be run."
        }

    } else {
        log.info "GenePA Cluster already exists."
    }
     if (runDrep) {
         all_genomes = selected_genomes.flatMap { dir -> file("${dir}/*.fa") }.collect()
         drep_result = run_drep_dereplication(all_genomes)
     }


        // prokka_result = run_prokka(bins_ch)
        // prokka_collect = prokka_result.collect()
//        ppanggolin_result = run_ppanggolin(prokka_collect)
        //kofamscan_result = run_kofamscan_annotation(ppanggolin_result)
        //vfdb_result = run_VFDB_annotation(ppanggolin_result)
        //card_result = run_rgi_CARD_annotation(ppanggolin_result)
        //dbcan_result = run_dbCAN_annotation(ppanggolin_result)
        //all_genomes = selected_genomes.flatMap { dir ->
        //    file("${dir}/*.fa")
        //}.collect()
        //skani_result = run_skani_annotation(all_genomes)

        // Visualization part
        metadata_ch = Channel.fromPath("${params.metadata}")
       // genePA will be updated to provide PA and  count         
        //genePA_result = run_genePA_cluster(ppanggolin_result, metadata_ch)
        //kofamscan_vis_result = run_kofamscan_visualization(kofamscan_result, metadata_ch)


if (params.run_visualization) {

        kofamscan_vis_result = run_kofamscan_visualization(
    kofamscan_result[0],
    metadata_ch          
)
        //vfdb_vis_result = run_VFDB_visualization(vfdb_result.gene_pa, metadata_ch)
        vfdb_vis_result = run_VFDB_visualization(vfdb_gene_pa, vfdb_gene_count, metadata_ch)
        //card_vis_result = run_rgi_CARD_visualization(card_result.gene_pa, metadata_ch)
        card_vis_result = run_rgi_CARD_visualization(card_gene_pa, card_gene_count, metadata_ch)     
        dbcan_vis_result = run_dbCAN_visualization(dbcan_count_matrix, dbcan_pa_matrix, metadata_ch)
        //skani_vis_result = run_skani_visualization(skani_result[0], metadata_ch)
        skani_vis_result = run_skani_visualization(skani_fullmatrix, metadata_ch)
        //defensefinder_vis_result = run_defensefinder_visualization(defensefinder_gene_pa, defensefinder_gene_count, metadata_ch)

        scoary2_result = run_scoary2(ppanggolin_result, metadata_ch)

        vis_results = Channel.empty()
            .mix(kofamscan_vis_result)
            .mix(vfdb_vis_result)
            .mix(card_vis_result)
            .mix(dbcan_vis_result)
            .mix(skani_vis_result)
            .mix(scoary2_result)
            //.mix(defensefinder_vis_result) 
            .collect()
} else {
    log.info "metacol is not specified. Visualization step will be skipped."
}
        // Shiny dashboard creation and MultiQC execution after all visualizations are done
        //shiny_dashboard = create_shiny_dashboard(vis_results)
        //run_multiqc(vis_results, shiny_dashboard)
    

}
