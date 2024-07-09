#!/usr/bin/env nextflow
// this script was optimized apptainer in any compurational environment. 
nextflow.enable.dsl=2
//params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"
params.db_baseDir = "/opt/database"
//params.scripts_baseDir ="/data2/leehg/OMD3/Generated_pipeline/scratch/tools/microbiome_analysis/scripts/"
////////////////////////
params.run_id = UUID.randomUUID().toString() // if you don't wnat random id,  set --run_id
////////////////////////
params.cpus = 76
params.outputDirCheckM2 = 'bin_assess_checkm2'
params.outputDirGUNC = 'bin_assess_gunc'

params.outdir_Base = "${launchDir}" // modify this if you want output files stored any other directory. 
params.inputDir = "${params.outdir_Base}/results/metagenome/ASSEMBLY_BINNING/final_bins"
def getAbsolutePath(String path) {
    def file = new File(path)
    return file.getAbsolutePath()
}
//params.inputDir = getAbsolutePath(params.inputDir)



params.outdir = "${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT"

params.GUNCdb = "${params.db_baseDir}/gunc/gunc_db_progenomes2.1.dmnd"
params.checkm2db="${params.db_baseDir}/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
//params.cpus = 64

if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir."
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : assembled genomic fna files."
    }
}
log.info """\
    BIN_ASSESSMENT : assess genome quality and assign taxonomy
    =================================================
    CheckM2                        : No option , database 1.0.2
    GUNC                           : No option , database progenomes2.1 diamond 
    GTDB r220 version              : 
    Final metadata file is created : quality_taxonomy_combined_${params.run_id}.csv
                                   : contains GUNC CheckM2 and GTDB information 
    cpus used                      : ${params.cpus}
    
    
    """
    .stripIndent(true)


process runCheckM2 {
    publishDir "${params.outdir}/checkm2_${params.run_id}", mode: 'copy'   
    //conda "$HOME/miniforge3/envs/checkm2"
    cpus params.cpus
    errorStrategy 'retry'
    maxRetries 3    

    input:
    path bin_dir

    output:
    path "${params.outputDirCheckM2}"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate checkm2 
    time checkm2 predict --threads ${params.cpus} --input ${params.inputDir} \
    --output-directory ${params.outputDirCheckM2} -x fa --database_path ${params.checkm2db}
    """
}

process runGUNC {
    publishDir "${params.outdir}/gunc_${params.run_id}", mode: 'copy'
    //conda "$HOME/miniforge3/envs/checkm2"
    cpus params.cpus

    input:
    path checkm2faa_dir 

    output:
    path "${params.outputDirGUNC}"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate checkm2     
    mkdir ${params.outputDirGUNC}
    time gunc run -d ${checkm2faa_dir}/protein_files -g -t ${task.cpus} -o ${params.outputDirGUNC} \
    --db_file ${params.GUNCdb} -e .faa
    """
}

process combineFiles {
    publishDir "${params.outdir}", mode: 'copy'
    //containerOptions = "--bind ${params.inputDir}"
    //runOptions = "--bind ${launchDir}/${params.inputDir}"
    cpus params.cpus

    //conda '$HOME/miniforge3/envs/checkm2'
    //apptainer.runOptions = "--bind ${params.inputDir}:${params.inputDir}"
    input:
    path checkm2Files
    path guncFiles
    path inputdir
    output:
    path "combined_report_${params.run_id}.tsv"// , emit : combined_report
    path "bins_quality_passed_${params.run_id}"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate checkm2     
    cp ${checkm2Files}/quality_report.tsv ./ 
    cp ${guncFiles}/GUNC.progenomes_2.1.maxCSS_level.tsv ./ 
    python ${params.scripts_baseDir}/Checkm2_GUNC_combine_quality_pass.py quality_report.tsv GUNC.progenomes_2.1.maxCSS_level.tsv combined_report_${params.run_id}.tsv
    mkdir bins_quality_passed_${params.run_id}
    # get only at least medium quality genomes based on MIMAG and CheckM2
    col_num=\$(head -1 combined_report_${params.run_id}.tsv | awk -F'\\t' '{for(i=1;i<=NF;i++) if (\$i=="medium_quality.pass") print i}')
    awk -v col="\$col_num" -F'\\t' '\$col=="True" {print \$1}' combined_report_${params.run_id}.tsv | sed -e 's/\$/.fa/g' > passlist
    #cat passlist
    echo \$PWD
    cp \$(grep -Ff passlist <(ls ${inputdir}/* )) bins_quality_passed_${params.run_id}
    
    """
}
    //awk -F'\t' '\$NF=="True" {print \$1}' combined_report.tsv  | sed -e 's/\$/.fa/g' > passlist

process gtdbtk {
    publishDir "${params.outdir}", mode: 'copy'
    //conda '$HOME/miniforge3/envs/gtdbtk-2.3.2'
    cpus params.cpus

    input:
    path combined_report
    path bins_quality_passed
    
    output:
    
    path "gtdb_outdir_${params.run_id}"
    path "quality_taxonomy_combined_${params.run_id}.csv"
    path "gtdbtk_r220_${params.run_id}.tsv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate gtdbtk-2.4.0
    export GTDBTK_DATA_PATH=${params.db_baseDir}/gtdb/release220/release220
    time gtdbtk classify_wf \
    --mash_db ${params.db_baseDir}/gtdb/release220/mash_220.db.msh --genome_dir ${bins_quality_passed}  -x fa  \
    --cpus ${task.cpus} --pplacer_cpus ${task.cpus} --out_dir gtdb_outdir_${params.run_id}
    awk -F'\\t' 'NR==1 || \$4=="True" {print \$0}'  ${combined_report} > quality_tmp.tsv
    if [[ -f gtdb_outdir_${params.run_id}/gtdbtk.bac120.summary.tsv ]]; then
        cat gtdb_outdir_${params.run_id}/gtdbtk.bac120.summary.tsv | grep classification | head -n 1 > gtdb_combined_tmp.tsv
    elif [[ -f gtdb_outdir_${params.run_id}/gtdbtk.ar53.summary.tsv ]]; then
        cat gtdb_outdir_${params.run_id}/gtdbtk.ar53.summary.tsv | grep classification | head -n 1 > gtdb_combined_tmp.tsv
    fi
    for file in gtdb_outdir_${params.run_id}/gtdbtk.*.summary.tsv; do
        if [[ -f "\$file" ]]; then
            cat "\$file" | tail -n+2  >> gtdb_combined_tmp.tsv
        fi
    done
    cut -f1,2  gtdb_combined_tmp.tsv > gtdb_classification_only.tsv
    mv gtdb_combined_tmp.tsv gtdbtk_r220_${params.run_id}.tsv

    micromamba activate checkm2 
    python ${params.scripts_baseDir}/GTDB_add2_check2gunc.py quality_tmp.tsv gtdb_classification_only.tsv quality_taxonomy_combined_${params.run_id}.csv
    """
}
//TBU :
//TBU : merge GTDB file together to output of GUNC +  Checkm2 
// cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep classification  | head -n 1  |cut -f1,2 > gtdb_combined_tmp.tsv
//     cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep -v classification |sort -u |cut -f1,2 >> gtdb_combined_tmp.tsv
//     cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep classification | head -n1  > gtdbtk_r214_${params.run_id}.tsv
//     cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep -v classification | sort -u >> gtdbtk_r214_${params.run_id}.tsv
    process createFinalFile {
    publishDir "${launchDir}", mode: 'copy'

    input:
    path A
    path B
    path C

    output:
    path "quality_taxonomy_combined.csv"

    script:
    """
    if [ ! -f "${launchDir}/quality_taxonomy_combined.csv" ]; then
        head -n 1 quality_taxonomy_combined_${params.run_id}.csv > quality_taxonomy_combined.csv
    fi
    tail -n +2 quality_taxonomy_combined_${params.run_id}.csv >> quality_taxonomy_combined.csv
    """
}

workflow {
    input_ch = Channel.fromPath(params.inputDir, type: 'dir', checkIfExists: true)
    checkm2_result = runCheckM2(input_ch)
    gunc2_result = runGUNC(checkm2_result)
    qc_filter_result = combineFiles(checkm2_result, gunc2_result, input_ch)
    result_ch=gtdbtk(qc_filter_result)
    createFinalFile(result_ch)
}
