#!/usr/bin/env nextflow
// this script was optimized for K-BDS hosted by KISTI in South Korea. 
nextflow.enable.dsl=2
params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

params.outdir_Base = "$HOME" // modify this if you want output files stored any other directory. 
params.inputDir = "${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT/bins_quality_passedFinal"
params.outdir = "${params.outdir_Base}/results/metagenome/GENOMAD"

params.cpus = 8

if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir."
} else {
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with .fa files."
    }
}


process run_genomad {
    publishDir "${params.outdir}" , mode : 'copy'
    conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    maxForks 4

    input:
    tuple val(id), path(genome)
    output:
    path "*"

    script:
    """
    source activate genomad
    genomad end-to-end  ${genome} ${id}_genomad ${params.db_baseDir}/genomad/genomad_db/ --threads ${params.cpus} --cleanup --splits 4
    """
}


workflow{
    bins_ch = Channel.fromFilePairs("${params.inputDir}/*.fa",size:1)
    run_genomad(bins_ch)
}
