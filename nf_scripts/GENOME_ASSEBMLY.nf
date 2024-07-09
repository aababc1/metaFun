#!/usr/bin/env nextflow
// this script was optimized for K-BDS hosted by KISTI in South Korea. 
nextflow.enable.dsl=2 
params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.inputReadType="" // should be one of short,long, SAG
params.longReadType="--pacbio-hifi" // should be one of  --pacbio-hifi and --pacbio-raw. default is pacbio-hifi

params.inputDir = ""
params.outdir_Base = "$HOME"
params.outdir =  "${params.outdir_Base}/results/genome/GENOME_ASSEMBLY/${params.inputReadType}"
params.cpus = 32
params.assembly_info = "${params.outdir_Base}/results/genome/GENOME_ASSEMBLY/assemblies_info.csv"

if (!new File(params.inputDir).exists()){
    error "Input directory dose not exist: ${params.inputDir}. Please specify a valid input directory."
}
if (!(['short', 'long', 'SAG'].contains(params.inputReadType))) {
    error "Invalid input read data type: ${params.inputReadType}. Please select one of the following: 'short', 'long', 'SAG'"
}
if (!(['--pacbio-hifi', '--pacbio-raw'].contains(params.longReadType))) {
    error "Invalid input read data type: ${params.longReadType}. Please select one of the following: '--pacbio-hifi', '--pacbio-raw'"
}

log.info """\
    De novo assembly, binning and binning refinement 
    =================================================
    inputReadType              : ${params.inputReadType}
    inputDir                   : ${params.inputDir}
    outdir                     : ${params.outdir}
    """
    .stripIndent(true)

process short_assembly {
    publishDir "${params.outdir}", mode : 'copy'
    publishDir "${params.outdir}/fastp", mode : 'copy', pattern: "*_fastp.json"
    conda "$HOME/miniforge3/envs/GENOME_ASSEMBLY" 
    maxForks 1
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(prefix), path(reads)

    output:
    tuple val(prefix), path("${prefix}.fa"), path("${prefix}_fastp.json")
    path("${prefix}_info.txt"), emit: short_info_files


    script:
    """
    source activate GENOME_ASSEMBLY
    fastp -i ${reads[0]} -I  ${reads[1]} -o ${prefix}_fastp_1.fastq.gz -O ${prefix}_fastp_2.fastq.gz \
    -w ${params.cpus} -j ${prefix}_fastp.json 

    spades.py -1 ${prefix}_fastp_1.fastq.gz -2 ${prefix}_fastp_2.fastq.gz -o ${prefix}_assembly -t ${params.cpus}
    cp ${prefix}_assembly/scaffolds.fasta ./${prefix}.fa
    echo "${prefix},Isolate_short" > ${prefix}_info.txt 

    """    
}

process long_assembly {
    publishDir "${params.outdir}", mode : 'copy'
    conda "$HOME/miniforge3/envs/GENOME_ASSEMBLY" 
    maxForks 1
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(fastq)
    
    output:
    path("*.fa")
    path("*_info.txt"), emit: long_info_files

    script:
    """
    source activate GENOME_ASSEMBLY
    file_name=\$(basename ${fastq})
    if [[ "\$file_name" == *.fastq.gz ]]; then
        file_prefix=\${file_name%.fastq.gz}
    elif [[ "\$file_name" == *.fastq ]]; then
        file_prefix=\${file_name%.fastq}
    fi

    flye ${params.longReadType}  ${fastq} --out-dir \${file_prefix}_assembly -t ${params.cpus}
    cp \${file_prefix}_assembly/assembly.fasta ./\${file_prefix}.fa 
    echo "\${file_prefix},Isolate_long" > \${file_prefix}_info.txt 

    """
}
process SAG_assembly { 
    publishDir "${params.outdir}", mode : 'copy'
    conda "$HOME/miniforge3/envs/GENOME_ASSEMBLY" 
    maxForks 1
    //errorStrategy 'retry'
    //maxRetries 1

    input:
    tuple val(prefix), path(reads)
    output:
    tuple val(prefix), path("${prefix}.fa"), path("${prefix}_fastp.json")
    path("${prefix}_info.txt"), emit: sag_info_files

    script:
    """
    source activate GENOME_ASSEMBLY
    fastp -i ${reads[0]} -I  ${reads[1]} -o ${prefix}_fastp_1.fastq.gz -O ${prefix}_fastp_2.fastq.gz \
    -w ${params.cpus} -j ${prefix}_fastp.json 
    
    spades.py -1 ${prefix}_fastp_1.fastq.gz -2 ${prefix}_fastp_2.fastq.gz -o ${prefix}_assembly \
    --sc --careful --disable-rr --disable-gzip-output -t ${params.cpus}
    cp ${prefix}_assembly/contigs.fasta ./${prefix}.fa
    echo "${prefix},SAG" > ${prefix}_info.txt 

    """    
}

// to modify 

process generate_csv {
    publishDir "${params.outdir}", mode : 'copy'

    input:
    path(info_files)

    output:
    path("assemblies_info.csv")

    script:
    """
    source activate GENOME_ASSEMBLY

    if [[ -f assemblies_info.csv ]]; then
        # Append without header
        for file in ${info_files}; do
            tail -n +2 \$file >> assemblies_info.csv
        done
    else
        # Create new file with header
        echo "accession,inputReadType" > assemblies_info.csv
        cat ${info_files.join(' ')} >> assemblies_info.csv
    fi
    """
}

workflow{
    // read_pair = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fastq.gz}", checkIfExists: true)
    // read_long = Channel.fromPath("${params.inputDir}/*.{fastq,fastq.gz}", checkIfExists: true)

    if (params.inputReadType == "short") {
        read_input = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fastq.gz}", checkIfExists: true)
        short_assembly(read_input).short_info_files.set {read_info}
        //short_assembly.short_info_files.set {read_info}
    } else if (params.inputReadType == "long") {
        read_input = Channel.fromPath("${params.inputDir}/*.{fastq,fastq.gz}", checkIfExists: true)
        long_assembly(read_input).long_info_files.set {read_info}
        //long_assembly.long_info_files.set {read_info}
    } else if (params.inputReadType == "SAG") {
        read_input = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fastq.gz}", checkIfExists: true)
        SAG_assembly(read_input).sag_info_files.set {read_info}
        //SAG_assembly.sag_info_files.set {read_info}
    }
    generate_csv(read_info)
    // if (params.inputReadType == "short"){
    //     //read_pair = fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fastq.gz}", checkIfExists: true)
    //     short_assembly(read_pair)
    // } else if (params.inputReadType == "long"){
    //     //read = Channel.fromPath("${params.inputDir}/*.{fastq,fastq.gz}", checkIfExists: true)
    //     long_assembly(read_long)
    // } else if (params.inputReadType == "SAG"){
    //     //read_pair = fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fastq.gz}", checkIfExists: true)
    //     SAG_assembly(read_pair)
    // }
}