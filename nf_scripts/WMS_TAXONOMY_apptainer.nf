#!/usr/bin/env nextflow

nextflow.enable.dsl=2 

params.db_baseDir = "/opt/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.inputDir = "${launchDir}/results/metagenome/RAWREAD_QC/read_filtered"
params.outdir =  "${launchDir}/results/metagenome/WMS_TAXONOMY"
params.metadata = ""
params.sampleIDcolumn = 1
params.analysiscolumn = 0
params.confidence_filter = 0.1  // default confidence value for kraken2
params.relab_filter = 0.0001  // optimal filtration for kraken2 with bracken 0.01% 1is maximum 
params.cpus = 15
params.kraken_method = 'default'
params.profiler = 'sylph'
params.sylph_abundance_type = 'relative_abundance'  

if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir."
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : paired-end read files."
    }
}

if (params.metadata.isEmpty()) {
    error "Metadata file is not specified. Please specify.  Please specify a valid metadata in CSV(comma , separated format) file with --metadata."
}
if (params.sampleIDcolumn <= 0) {
    error "Sample ID column is not specified. Please specify a sample ID column number(integer) with --sampleIDcolumn."
}
if (params.analysiscolumn <= 0) {
    println """
    Analysis ID column is not specified. You would want to analysis it separately with script. \
    do nextflow run WMS_TAXA_stat.nf --columnID {column number}
    Or specify a sample ID column number(integer) with --analysiscolumn you want to analyze."
    """
}

def ANSI_RESET = "\u001B[0m"
def ANSI_BLUE = "\033[38;2;8;70;250m"  // WMS_TAXONOMY 색상
def ANSI_YELLOW = "\u001B[33m"
def ANSI_GREEN = "\u001B[32m"
def ANSI_CYAN = "\033[38;2;0;255;255m"
def ANSI_RED = "\u001B[31m"
def ANSI_ORANGE = "\033[38;2;255;149;0m"

log.info """\
 Usage:
    nextflow run nf_scripts/${ANSI_BLUE}WMS_TAXONOMY_apptainer.nf${ANSI_RESET} --inputDir "\${filtered_reads_dir}" --metadata "\${metadata_file}" --sampleIDcolumn "\${column_number}" --profiler "\${taxonomy_tool}"

${ANSI_BLUE}WMS_TAXONOMY : Taxonomic profiling of whole metagenome sequencing data${ANSI_RESET}
=========================================================
 ${ANSI_BLUE}Kraken2${ANSI_RESET}           : Fast taxonomic classification using k-mer matching
 ${ANSI_BLUE}Bracken${ANSI_RESET}           : Abundance estimation from Kraken2 results
 ${ANSI_BLUE}Sylph${ANSI_RESET}             : Fast and accurate taxonomic profiling

${ANSI_RED}Required Parameters:${ANSI_RESET}
------------------
${ANSI_RED}--inputDir${ANSI_RESET}         : Directory containing quality filtered paired-end reads
${ANSI_RED}--metadata${ANSI_RESET}         : CSV file with sample metadata (comma-separated format)
${ANSI_RED}--sampleIDcolumn${ANSI_RESET}   : Column number in metadata file containing sample IDs
${ANSI_RED}--profiler${ANSI_RESET}         : Taxonomic profiler to use (kraken2 or sylph)

${ANSI_YELLOW}Optional Parameters:${ANSI_RESET}
------------------
${ANSI_YELLOW}--analysiscolumn${ANSI_RESET}   : Column number for statistical analysis grouping
${ANSI_YELLOW}--kraken_method${ANSI_RESET}    : Kraken2 execution method (default: ${params.kraken_method})
${ANSI_YELLOW}--confidence_filter${ANSI_RESET} : Confidence threshold for Kraken2 (default: ${params.confidence_filter})
${ANSI_YELLOW}--relab_filter${ANSI_RESET}     : Relative abundance filter for Bracken (default: ${params.relab_filter})
${ANSI_YELLOW}--sylph_abundance_type${ANSI_RESET} : Abundance type for Sylph (default: ${params.sylph_abundance_type})
${ANSI_YELLOW}--cpus${ANSI_RESET}             : Number of CPUs to use (default: ${params.cpus})
${ANSI_YELLOW}--outdir${ANSI_RESET}           : Output directory (default: ${params.outdir})

${ANSI_CYAN}Output Files:${ANSI_RESET}
------------------
${ANSI_CYAN}- kraken2/${ANSI_RESET}          : Kraken2 classification results
${ANSI_CYAN}- bracken/${ANSI_RESET}          : Bracken abundance estimates
${ANSI_CYAN}- sylph/${ANSI_RESET}            : Sylph profiling results
${ANSI_CYAN}- phyloseq/${ANSI_RESET}         : Phyloseq objects for R analysis
${ANSI_CYAN}- stats_analysis/${ANSI_RESET}    : Statistical analysis results (if --analysiscolumn specified)

Final output is a Phyloseq object that can be used with the ${ANSI_ORANGE}INTERACTIVE_TAXONOMY${ANSI_RESET} module for interactive visualization.

"""
.stripIndent(true)


process kraken2_run{
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    tag "kraken2  ${accession}"
    cpus params.cpus
    
    input:
    tuple val(accession), path(reads)

    output:
    val(accession)
    path("${accession}_kraken2.report")

    script:
    def db_path = params.kraken_method == 'memory-mapping' ? "/dev/shm/kraken2_GTDB/" : "${params.db_baseDir}/kraken2_GTDBr220/"
    def memory_mapping = params.kraken_method == 'memory-mapping' ? "--memory-mapping" : ""

    """
    kraken2 --db ${db_path} --output ${accession}_kraken2.out --report ${accession}_kraken2.report \
    --confidence ${params.confidence_filter} --gzip-compressed --paired ${memory_mapping} \
    --threads ${task.cpus} ${reads[0]} ${reads[1]}
    """
}
    /*
    if(params.kraken_method == 'default')
    """
    kraken2 --db ${params.db_baseDir}/kraken2_GTDBr220/ --output ${accession}_kraken2.out --report ${accession}_kraken2.report \
    --confidence 0.25 --gzip-compressed --paired  \
    --threads ${task.cpus}  ${reads[0]} ${reads[1]} 

    """
    else if (params.kraken_method == 'memory-mapping')
    """
    kraken2 --db /dev/shm/kraken2_GTDB/ --output ${accession}_kraken2.out --report ${accession}_kraken2.report \
    --confidence 0.25 --gzip-compressed --paired  --memory-mapping \
    --threads ${task.cpus}  ${reads[0]} ${reads[1]} 
    """
    
}
*/

process bracken_run{
    publishDir "${params.outdir}/bracken", mode: 'copy'
    tag "bracken  ${accession}"
    containerOptions {"--bind ${launchDir}:${launchDir}"}

    input:
    val(accession)
    path(kraken_output)
    output:
    path "*"
    //path("${accession}_bracken.out")
    //awk -v sample=$sample '{if ($1 == sample "fastp_hg38_1") print $(NF - 2)}' ${params.multiqcFile}

    script:
    """
    #ls /data1/leehg/OMD_pipeline/year3/WMS_TAXONOMY/read_filtered
    ls ${launchDir}/results/metagenome/RAWREAD_QC/multiqc/multiqc_data/multiqc_general_stats.txt 
    
    len=\$(awk -v accession=${accession} -F"\\t" '{if (\$1 ~ accession "_fastp.*_1") print \$(NF -2)}' \
    ${launchDir}/results/metagenome/RAWREAD_QC/multiqc/multiqc_data/multiqc_general_stats.txt )
    read_length=0
    if [ "\$len" -le 125 ]; then
        read_length=100
    elif [ "\$len" -le 175 ]; then
        read_length=150
    elif [ "\$len" -le 225 ]; then
        read_length=200
    elif [ "\$len" -le 275 ]; then
        read_length=250    
    else
        read_length=300
    fi
    
    /scratch/tools/microbiome_analysis/program/Bracken/bracken -d ${params.db_baseDir}/kraken2_GTDBr220/ \
    -i ${kraken_output} -o ${accession}_bracken.out -r \${read_length} -l S 
    # relative abundance filtering 
    awk '\$NF >= ${params.relab_filter}' ${accession}_bracken.out > ${accession}_bracken.out.tmp
    mv ${accession}_bracken.out.tmp ${accession}_bracken.out

    """

}

process phyloseq_creation {
    publishDir "${params.outdir}/phyloseq", mode: 'copy'
    // below line is tmp . 

    containerOptions {"--bind /data2/leehg/OMD3/Generated_pipeline/scratch/tools/microbiome_analysis/scripts:/scratch/tools/microbiome_analysis/scripts"}

    input: 
    path(bracken_files)
    path(metadata)

    output:
    path("phyloseq_object.RDS")

    //3 custome script : 2 for bracken relab data manage, 1 for phyloseq object creation 
    //1 database information file : taxID to species of kraken database 
    //   BrackentaxID_GTDBtaxa_completed.csv

    script: 
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    python ${params.scripts_baseDir}/combine_bracken_outputs.py --files *bracken.out -o combined_bracken_out
    python ${params.scripts_baseDir}/extract_frac_combined_bracken.py combined_bracken_out bracken_species_relab.tsv
    sed -i 's/_bracken.out_frac//g' bracken_species_relab.tsv
    
    Rscript ${params.scripts_baseDir}/phyloseq_creation.R  bracken_species_relab.tsv  \
    ${params.scripts_baseDir}/BrackentaxID_GTDBtaxa_completed.csv ${params.metadata} ${params.sampleIDcolumn}\
    ./

    """
}


process statistical_analysis {
    publishDir "${params.outdir}/stats_analysis", mode: 'copy'
    tag "stats_analysis"
    when: params.analysiscolumn >= 0

    input:
    path(phyloseq_object)

    output:
    path "*" // Adjust output as necessary

    script:
    def adjusted_analysiscolumn = params.analysiscolumn -1 

    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment 
    Rscript ${params.scripts_baseDir}/phyloseq_WMS_taxa_analysis.R ${phyloseq_object} \
    WMS_taxonomy_analysis ${adjusted_analysiscolumn}
    """
}

process sylph_sketch_all {
    publishDir "${params.outdir}/sylph", mode: 'copy'
    tag "sylph_sketch_all"
    cpus params.cpus
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.paired.sylsp")
    
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate sylph061    
    sylph sketch -t ${task.cpus} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -S ${sample_id} \
        -c 200
    
    """
}

process sylph_process_all {
    publishDir "${params.outdir}/sylph", mode: 'copy'
    cpus params.cpus

    input:
    path('*')  
    
    output:
    path("merged_sylph_species.tsv")
    path("all.profile-sylph.tsv")
    
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate sylph061    

    # 1. Profiling using database 
    sylph profile *.sylsp ${params.db_baseDir}/sylph/sylph_c200_gtdbr220.syldb \
        -t ${task.cpus} \
        -o all.profile-sylph.tsv

    # 2. Convert to sylphmpa format
    python ${params.scripts_baseDir}/sylph-utils/sylph_to_taxprof.py \
        -s all.profile-sylph.tsv \
        -m ${params.scripts_baseDir}/gtdbr220_taxonomy_list_TaxID.tsv

    # 3. Merge sylphmpa files

    python ${params.scripts_baseDir}/sylph-utils/merge_sylph_taxprof.py \
        -o merged_sylphmpa.tsv \
        --column ${params.sylph_abundance_type} \
        *.sylphmpa

    # 4. Convert to TaxID
    python ${params.scripts_baseDir}/convert_species_TAXID_merged_syhph_mpa.py \
        -m ${params.scripts_baseDir}/gtdbr220_taxonomy_list_TaxID.tsv \
        -i merged_sylphmpa.tsv \
        -o merged_sylph_species.tsv
    """
}


process phyloseq_creation_sylph {
    publishDir "${params.outdir}/phyloseq", mode: 'copy'
    
    input:
    path(merged_species)
    path(metadata)
    
    output:
    path("phyloseq_object_sylph.RDS")
    
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    
    Rscript ${params.scripts_baseDir}/phyloseq_creation.R \
        ${merged_species} \
        ${params.scripts_baseDir}/BrackentaxID_GTDBtaxa_completed.csv \
        ${metadata} \
        ${params.sampleIDcolumn} \
        ./
        mv phyloseq_object.RDS phyloseq_object_sylph.RDS
    """
}




workflow{
    reads_ch = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fq}{.gz,}", checkIfExists: true)
                .map { accession, reads -> 
                    // Split the accession name to remove the '_fastp' part
                    //def sam_accession = accession.split('_')[0]
                    def sam_accession = accession.replaceAll(/_fastp.*$/, '')
                    tuple(sam_accession, reads) 
                }                
    metadata_ch = Channel.fromPath("${params.metadata}")

    reads_ch.view()


    if (params.profiler == 'kraken2') {
        kraken2_result = kraken2_run(reads_ch)
        bracken_result = bracken_run(kraken2_result)
        phyloseq_input = bracken_result.collect()
        phyloseq_ch = phyloseq_creation(phyloseq_input, metadata_ch)
    }
    else if (params.profiler == 'sylph') {
        sylph_sketches = sylph_sketch_all(reads_ch)
        sylph_results = sylph_process_all(sylph_sketches.collect())
        phyloseq_ch = phyloseq_creation_sylph(sylph_results[0], metadata_ch)
    }
    


    //phyloseq_ch = phyloseq_creation(phyloseq_input,metadata_ch)
    if (params.analysiscolumn > 0) {
        statistical_analysis(phyloseq_ch)
    }
    //reads_ch = Channel.fromFilePairs("", checkIfExists : true )
}