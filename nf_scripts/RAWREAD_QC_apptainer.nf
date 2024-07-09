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

log.info """\
    Quality contorl of raw reads of whole metagenome sequencing data. 
    ===================================
    output_directory : ${params.outdir}
    
    human_index : ${params.human_index}
    
    used_index  : ${params.custom_index}
                
    filter      : ${params.filter}

                If you want to use other index , specify it with --filter \$value in command line 
                You need to index genomes with in specific folder 
                ${filter} is the prefix name of your gneome index. 
                Index files should be in the database/host_genome subdirectory 
                \$ mkdir /database/host_genome/\${filter} 
                \$ bowtie2-build \$your_genome database/\${filter}/\${filter}  

    inputDir    : ${params.inputDir}
    outdir      : ${params.outdir}
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
    fastp_out = fastp_only(reads_ch)
    qc_results = Channel.empty()
                    .mix(fastqc_raw_out)
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
