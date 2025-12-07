#!/usr/bin/env nextflow
// this script was optimized apptainer in any compurational environment. 
nextflow.enable.dsl=2

//params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.db_baseDir = "/opt/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"
params.filter = 'human'
params.human_index = "${params.db_baseDir}/host_genome/GRCh38/GRCh38_p12.dna.primary_assembly.fa"
params.custom_index = "${params.db_baseDir}/host_genome/${params.filter}/${params.filter}"
//params.reads = "${baseDir}/raw_reads/*_{1,2}.fastq.gz"
params.outdir = "${launchDir}/results/metagenome/RAWREAD_QC"
params.cpus = 4
params.inputDir = 'raw_reads' 

if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir. \n \
    Input paired-end read files should be gzipped fastq files"
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : paired-end read files.\n \
        If there is problem with execution of assembly, check the fastq files with vdb-validate from sra-tools"
    }
}

def ANSI_RESET = "\u001B[0m"
def ANSI_YELLOW = "\u001B[33m"
def ANSI_GREEN = "\u001B[32m"
def ANSI_BLUE = "\u001B[34m"
def ANSI_RED= "\u001B[31m"

log.info """\

===========================================
    ${ANSI_RED}RAWREAD_QC${ANSI_RESET} - Raw Read Quality Control  
===========================================
(Apptainer Container with Nextflow)

Key Features:
  ✓ FastQC-based quality assessment
  ✓ fastp-based read quality control
  ✓ Host DNA filtering (Human/Custom/None)
  ✓ MultiQC generate HTML reports


\033[31mRequired Parameters:\033[0m 
------------------
\033[31m--inputDir\033[0m         : Directory containing raw read files
                     Files should be paired-end gzipped FASTQ files
                     with naming pattern: *_1.fastq.gz and *_2.fastq.gz

\033[33m--filter\033[0m           : Type of host reads to filter out
                     Available options:
                    \033[31mhuman\033[0m : Filter human reads (default)
                    \033[31mnone\033[0m :  \033[1;31mSkip host filtering\033[0m
                    \033[31mcustom\033[0m : Use custom host genome index
                          1. Prepare index: 
                             \033[1m\$ metafun -module PREPARE_CUSTOM_HOST -i <genome.fasta> -f <custom_name>\033[0m
                             Your custom genome is indexed by Bowtie2 using \033[31mcustom_name\033[0m. 
                          2. Use with: metafun -module RAWREAD_QC  --filter <custom_name>


Optional Parameters:
------------------
\033[33m--custom_index\033[0m     : Path to custom genome index for filtering
                     Used only when --filter=custom
\033[33m--outdir\033[0m           : Output directory 
                     (default: ${launchDir}/results/metagenome/RAWREAD_QC)
\033[33m--cpus\033[0m             : Number of CPUs to use (default: 4)


\033[38;2;0;255;255mOutput Files:\033[0m
-----------
\033[38;2;0;255;255m- read_filtered/\033[0m           : Directory containing filtered reads
\033[38;2;0;255;255m- fastqc_raw/\033[0m              : Raw read quality reports
\033[38;2;0;255;255m- fastqc_filtered/\033[0m         : Filtered read quality reports
\033[38;2;0;255;255m- multiqc/\033[0m                 : Combined quality reports


\033[38;2;180;180;180mNote:\033[0m Filtered reads from this module will be used in: 
\033[38;2;255;149;0mASSEMBLY_BINNING\033[0m 
\033[38;2;0;176;80mBIN_ASSESSMENT\033[0m 
\033[38;2;8;70;250mWMS_TAXONOMY\033[0m 
\033[38;2;112;48;160mWMS_FUNCTION\033[0m

    """
    .stripIndent(true)

workflow {
    reads_ch = Channel.fromFilePairs("${params.inputDir}/*{1,2}.{fastq,fq}{.gz,}", checkIfExists: true)
                .map {accession, reads -> 
                def sam_accession = accession.replaceAll(/[1,2]\.(fastq|fq)(\.gz)?$/, '')
                tuple(sam_accession, reads)
                }


    //read_pairs = Channel.fromFilePairs(params.reads)
    fastqc_raw_out = fastqc_raw(reads_ch)
    if (params.filter == "none") {
    //fastp_out = fastp_only(reads_ch)
    fastp_out = fastp_with_publish(reads_ch)
    fastqc_filtered_out=fastqc_filtered(fastp_out.map { it -> tuple(it[0], it[1]) })
    qc_results = Channel.empty()
                    .mix(fastqc_raw_out)
                    .mix(fastqc_filtered_out)
                    .mix(fastp_out.map{ it[2] })
                    .collect()
    }
    else {
    fastp_out = fastp(reads_ch)
    humanread_filter_out = humanread_filter(fastp_out.map { it -> tuple(it[0], it[1]) })
    fastqc_filtered_out = fastqc_filtered(humanread_filter_out)

    qc_results = Channel.empty()
                        .mix(fastqc_raw_out)
                        .mix(fastqc_filtered_out)
                        .mix(fastp_out.map{ it[2] })
                        .collect()
    }
    //fastp_out = fastp(reads_ch)
    //humanread_filter_out = humanread_filter(fastp_out.map { it -> tuple(it[0], it[1]) })
    //fastqc_filtered_out = fastqc_filtered(humanread_filter_out)

    //qc_results = Channel.empty()
    //                   .mix(fastqc_raw_out)
    //                    .mix(fastqc_filtered_out)
    //                    .mix(fastp_out.map{ it[2] })
    //                    .collect()
    multiQC(qc_results)

    //multiQC([fastqc_raw.out, fastqc_filtered.out])
}

process fastqc_raw {
    tag "fastqc ${sample_id}"
    //conda "$HOME/miniforge3/envs/RAWREAD_QC"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}_raw_fastqc")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    mkdir ${sample_id}_raw_fastqc
    fastqc ${reads} -o ${sample_id}_raw_fastqc -t ${task.cpus}
    """
}

process fastp {
    tag "fastp ${sample_id}"
    //publishDir "${params.outdir}/fastp", mode: 'copy'
    //conda "$HOME/miniforge3/envs/RAWREAD_QC"
    cpus params.cpus
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastp_*.fastq.gz"), path("${sample_id}_fastp.json") //, emit: json

    script:
    //if(params.filter == 'human')
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC    
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \
    --out1 ${sample_id}_fastp_1.fastq.gz --out2 ${sample_id}_fastp_2.fastq.gz \
    -w ${task.cpus} -j ${sample_id}_fastp.json
    """
    //else
    //    error "Invalid alignment mode: ${params.filter}"
    //-j ${sample_id}_fastp.json
}

/*
mode = 'tcoffee'

process align {
    input:
    path sequences

    script:
    if( mode == 'tcoffee' )
        """
        t_coffee -in $sequences > out_file
        """

    else if( mode == 'mafft' )
        """
        mafft --anysymbol --parttree --quiet $sequences > out_file
        """

    else if( mode == 'clustalo' )
        """
        clustalo -i $sequences -o out_file
        """

    else
        error "Invalid alignment mode: ${mode}"
}
*/


process fastp_with_publish {
    tag "fastp ${sample_id}"
    publishDir "${params.outdir}/read_filtered", mode: 'copy'
    cpus params.cpus
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastp_*.fastq.gz"), path("${sample_id}_fastp.json")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC    
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \
    --out1 ${sample_id}_fastp_1.fastq.gz --out2 ${sample_id}_fastp_2.fastq.gz \
    -w ${task.cpus} -j ${sample_id}_fastp.json
    """
}



process humanread_filter {
    tag "read filter based on ${params.filter} : ${sample_id}"
    //conda "$HOME/miniforge3/envs/RAWREAD_QC"
    publishDir "${params.outdir}/read_filtered", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastp_*.fastq.gz")

    script:
    if (params.filter=='human')
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    bowtie2 --very-sensitive -x ${params.human_index} -1 ${reads[0]} -2 ${reads[1]} \
    -p ${task.cpus} --un-conc-gz ${sample_id}_fastp_hg38 -S /dev/null
    mv ${sample_id}_fastp_hg38.1 ${sample_id}_fastp_hg38_1.fastq.gz
    mv ${sample_id}_fastp_hg38.2 ${sample_id}_fastp_hg38_2.fastq.gz
    """
    else if(params.filter != 'human')
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    bowtie2 --very-sensitive -x ${params.custom_index} -1 ${reads[0]} -2 ${reads[1]} \
    -p ${task.cpus} --un-conc-gz ${sample_id}_fastp_custom -S /dev/null
    mv ${sample_id}_fastp_custom.1 ${sample_id}_fastp_custom_1.fastq.gz
    mv ${sample_id}_fastp_custom.2 ${sample_id}_fastp_custom_2.fastq.gz
    """
    else 
    // Optionally handle the default case
    // This block can be left empty if no action is needed
    """
    echo "No matching condition found for filter: $filter"
    """ 
}

process fastqc_filtered {
    tag "${sample_id}"
    //conda "$HOME/miniforge3/envs/RAWREAD_QC"
    publishDir "${params.outdir}/fastqc_filtered", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}_raw_fastqc_filtered")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    mkdir ${sample_id}_raw_fastqc_filtered
    fastqc ${reads} -o ${sample_id}_raw_fastqc_filtered  -t ${task.cpus}
    """
}

process multiQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    //conda "$HOME/miniforge3/envs/RAWREAD_QC"

    input:
    path reports

    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC

    multiqc .
    """
}
