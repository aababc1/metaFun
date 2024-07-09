#!/usr/bin/env nextflow

nextflow.enable.dsl=2 

params.db_baseDir = "/opt/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.inputDir = "${launchDir}/results/metagenome/RAWREAD_QC/read_filtered"

params.outdir =  "${launchDir}/results/metagenome/WMS_TAXONOMY"
params.metadata = ""
params.sampleIDcolumn = 0
params.analysiscolumn = 0
params.relab_filter = 0.000001 
params.cpus = 15
params.kraken_method = 'default'

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
log.info """\
    De novo assembly, binning and binning refinement 
    =================================================
    MEGAHIT option              : {params.megahit_presets} # available values : default, meta-large, meta-sensitive

    SemiBin2 model              : {params.semibin2_mode} # defulat is human_gut
                                Available models in SemiBin2 : human_gut dog_gut ocean soil cat_gut human_oral mouse_gut pig_gut
                                built_environment wastewater chicken_caecum global
                                If you choose other values  : SemiBin2 will be executed with self-train mode.

    selectd SemiBin2 mode       : {chosenFilter}

    inputDir                    : {params.inputDir}
    outdir                      : {params.outdir}
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
    if(params.kraken_method == 'default')
    """
    kraken2 --db ${params.db_baseDir}/kraken2_GTDBr220/ --output ${accession}_kraken2.out --report ${accession}_kraken2.report \
    --confidence 0.25 --gzip-compressed --paired  \
    --threads ${task.cpus}  ${reads[0]} ${reads[1]} 

    """
    else(params.kraken_method == 'memory-mapping')
    """
    kraken2 --db /dev/shm/kraken2_GTDB/ --output ${accession}_kraken2.out --report ${accession}_kraken2.report \
    --confidence 0.25 --gzip-compressed --paired  --memory-mapping \
    --threads ${task.cpus}  ${reads[0]} ${reads[1]} 
    """

}


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

def adjusted_analysiscolumn = params.analysiscolumn -1 
process statistical_analysis {
    publishDir "${params.outdir}/stats_analysis", mode: 'copy'
    tag "stats_analysis"
    when: params.analysiscolumn >= 0

    input:
    path(phyloseq_object)

    output:
    path "*" // Adjust output as necessary

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment 
    Rscript ${params.scripts_baseDir}/phyloseq_WMS_taxa_analysis.R ${phyloseq_object}\
    WMS_taxonomy_analysis ${adjusted_analysiscolumn}
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
    //    reads_ch = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.fastq.gz", checkIfExists: true)
    //            .map { accession, reads -> 
    //                 // Split the accession name to remove the '_fastp' part
    //                def sam_accession = accession.split('_')[0]
    //                tuple(sam_accession, reads) 
    //            }
    reads_ch.view()
    kraken2_result = kraken2_run(reads_ch)
    bracken_result = bracken_run(kraken2_result)
    // without below line, accession information lost 
    phyloseq_input = bracken_result.collect()
    phyloseq_input.view()

    metadata_ch = Channel.fromPath("${params.metadata}")
    phyloseq_ch = phyloseq_creation(phyloseq_input,metadata_ch)
    statistical_analysis(phyloseq_ch)

    //reads_ch = Channel.fromFilePairs("", checkIfExists : true )
}