#!/usr/bin/env nextflow

nextflow.enable.dsl=2 

params.db_baseDir = "/opt/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.inputDir = "${launchDir}/results/metagenome/RAWREAD_QC/read_filtered"
params.outdir =  "${launchDir}/results/metagenome/WMS_FUNCTION"

params.metadata = "meta.csv"
params.sampleIDcolumn = 1
params.analysiscolumn = 0

params.metadatacolumn = ""
params.cpus = 36


if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir."
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : paired-end read files."
    }
}

// if (params.metadata.isEmpty()) {
//     error "Metadata file is not specified. Please specify.  Please specify a valid metadata in CSV(comma , separated format) file with --metadata."
// }
if (params.sampleIDcolumn <= 0) {
    error "Sample ID column is not specified. Please specify a sample ID column number(integer) with --sampleIDcolumn."
}

process humann3_run{
    publishDir "${params.outdir}/humann3", mode: 'copy'
    tag "humann3  ${accession}"
    errorStrategy 'retry'
    maxRetries 1
    cpus 6
    //containerOptions {"--bind ${params.db_baseDir}/humann3/metaphlan:/opt/conda/envs/humann3/lib/python3.10/site-packages/metaphlan/metaphlan_databases/"}

    input:
    tuple val(accession), path(reads)
    
    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate humann3

    if [[ "${reads[0]}" == *.gz ]]; then
        cat ${reads[0]} ${reads[1]} > ${accession}_merged.fastq.gz
        humann -i ${accession}_merged.fastq.gz -o ${accession}_humann3 --threads ${task.cpus} \
        --input-format fastq.gz --search-mode uniref90 --output-basename ${accession} --output-format tsv \
        --pathways metacyc --nucleotide-database ${params.db_baseDir}/humann3/chocophlan/ --protein-database ${params.db_baseDir}/humann3/uniref/
        rm ${accession}_merged.fastq.gz
    else
        cat ${reads[0]} ${reads[1]} > ${accession}_merged.fastq
        humann -i ${accession}_merged.fastq -o ${accession}_humann3 --threads ${task.cpus} \
        --input-format fastq --search-mode uniref90 --output-basename ${accession} --output-format tsv \
        --pathways metacyc --nucleotide-database ${params.db_baseDir}/humann3/chocophlan/ --protein-database ${params.db_baseDir}/humann3/uniref/
        rm ${accession}_merged.fastq
    fi


    """
}

/*

    cat ${reads[0]} ${reads[1]}  > ${accession}_merged.fastq.gz
    humann -i ${accession}_merged.fastq.gz  -o ${accession}_humann3 --threads ${task.cpus} \
    --input-format fastq.gz --search-mode uniref90 --output-basename ${accession} --output-format tsv \
    --pathways metacyc --nucleotide-database  ${params.db_baseDir}/humann3/chocophlan/ --protein-database ${params.db_baseDir}/humann3/uniref/

    rm ${accession}_merged.fastq.gz

*/

//https://github.com/DerrickWood/kraken2/issues/100

process humann3_parsing{
    publishDir "${params.outdir}/humann3_combined", mode: 'copy'
    tag "humann3 result parssing"

    input:
    path(humann_output)
    
    output:
    //path "*"
    path("humann_split_stratified_table/*")
    //path("${accession}_bracken.out")
    //awk -v sample=$sample '{if ($1 == sample "fastp_hg38_1") print $(NF - 2)}' ${params.multiqcFile}
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate humann3
    humann_join_tables -i ./ -o pathabund_join.tsv --file_name  _pathabundance.tsv -s 
    humann_renorm_table -i pathabund_join.tsv -o pathabund_join_renorm_cpm.tsv -u cpm
    humann_split_stratified_table -i pathabund_join_renorm_cpm.tsv  -o humann_split_stratified_table
    sed -i 's/_Abundance//g' humann_split_stratified_table/pathabund_join_renorm_cpm_unstratified.tsv

    """
}

process function_analysis {
    publishDir "${params.outdir}", mode: 'copy'

    input: 
    path(humann3_parsing_result)
    path(metadata)

    output:
    path("*")
    when: params.analysiscolumn >= 0

    //3 custome script : 2 for bracken relab data manage, 1 for phyloseq object creation 
    //1 database information file : taxID to species of kraken database 
    //   BrackentaxID_GTDBtaxa_completed.csv

    script: 
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    
    Rscript ${params.scripts_baseDir}/humann3_visualization.R \
    -i pathabund_join_renorm_cpm_unstratified.tsv \
    -out WMS_Function_result -m ${metadata} -mc ${params.analysiscolumn} -sc  ${params.sampleIDcolumn} > humann_analysis_inspect.log
    echo "loo at humann_analysis_inspect.log file  " 

    """
}
//-out WMS_Function_result -m ${launchDir}/${params.metadata} -mc ${params.analysiscolumn} -sc  ${params.sampleIDcolumn}


workflow{
    reads_ch = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fq}{.gz,}", checkIfExists: true)
                .map { accession, reads -> 
                    // Split the accession name to remove the '_fastp' part
                    //def sam_accession = accession.split('_')[0]
                    def sam_accession = accession.replaceAll(/_fastp.*$/, '')
                    tuple(sam_accession, reads) 
                }                
    reads_ch.view()
    metadata_ch=Channel.fromPath("${params.metadata}")
    humann3_result = humann3_run(reads_ch)
    humann3_combined = humann3_result.collect()
    humann3_parsing_result = humann3_parsing(humann3_combined)
    //function_analysis(humann3_parsing_result, "${params.metadata}")
    function_analysis(humann3_parsing_result, metadata_ch)
    // without below line, accession information lost 
    //phyloseq_input = bracken_result.collect()
    //phyloseq_input.view()

    //metadata_ch = Channel.fromPath("${params.metadata}")
    //phyloseq_creation(phyloseq_input,metadata_ch)

    //reads_ch = Channel.fromFilePairs("", checkIfExists : true )
}


