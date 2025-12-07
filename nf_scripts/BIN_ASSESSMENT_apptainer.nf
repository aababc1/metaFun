#!/usr/bin/env nextflow
// this script was optimized apptainer in any compurational environment. 
nextflow.enable.dsl=2

/*
params.workflow_name="BIN_ASSESSMENT"
manifest {
    name = 'BIN_ASSESSMENT'
}
*/
//params.db_baseDir = "/scratch/tools/microbiome_analysis/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"
params.db_baseDir = "/opt/database"
//params.scripts_baseDir ="/data2/leehg/OMD3/Generated_pipeline/scratch/tools/microbiome_analysis/scripts/"
////////////////////////
def runId = "${new Date().format('yyyyMMddHHmmss')}_${workflow.runName}"
params.run_id = runId


//params.run_id = UUID.randomUUID().toString() // if you don't wnat random id,  set --run_id
////////////////////////
params.cpus = 20
params.outputDirCheckM2 = 'bin_assess_checkm2'
params.outputDirGUNC = 'bin_assess_gunc'

params.outdir_Base = "${launchDir}" // modify this if you want output files stored any other directory. 
params.inputDir = "${params.outdir_Base}/results/metagenome/ASSEMBLY_BINNING/final_bins"

params.metadata = null
params.accession_column = 1
params.data_type='MAG'

params.pass_quality = "QS50.pass"
def valid_quality_filters = ['medium_quality.pass', 
                           'high_quality.pass', 
                           'medium_quality_gunc.pass', 
                           'high_quality_gunc.pass', 
                           'QS50.pass', 
                           'QS50_gunc.pass',
                           'all']  // 

// filter check . 
if (!valid_quality_filters.contains(params.pass_quality)) {
    error """Invalid quality filter: ${params.pass_quality}
            Valid options are: ${valid_quality_filters.join(', ')}"""
}

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
println("https://metafun-doc.readthedocs.io/en/latest/workflows/BIN_ASSESSMENT.html\n")


def getAbsolutePath(String path) {
    def file = new File(path)
    return file.getAbsolutePath()
}


//params.inputDir = getAbsolutePath(params.inputDir)


params.outdir = "${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT"

params.GUNCdb = "${params.db_baseDir}/gunc/gunc_db_progenomes2.1.dmnd"
params.checkm2db="${params.db_baseDir}/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
//params.cpus = 64

def ANSI_RESET = "\u001B[0m"
def ANSI_YELLOW = "\u001B[33m"
def ANSI_GREEN = "\u001B[32m"
def ANSI_BLUE = "\u001B[34m"
def ANSI_RED= "\u001B[31m"

params.help=null 
//${custom_logo}
//Visit documentation to find out description of this workflow: 
//https://metafun-doc.readthedocs.io/en/latest/workflows/BIN_ASSESSMENT.html

if (params.help) {
log.info """
Usage : 
    nextflow run nf_scripts/BIN_ASSESSMENT_apptainer.nf --metadata "\${your metadata}" --accession_column "\${your accessioncolumn}"
${ANSI_YELLOW}BIN_ASSESSMENT : assess genome quality and assign taxonomy${ANSI_RESET}
=========================================================
 ${ANSI_GREEN}CheckM2${ANSI_RESET}            : No option, database 1.0.2
 ${ANSI_GREEN}GUNC${ANSI_RESET}               : No option, database progenomes2.1 diamond
 ${ANSI_GREEN}GTDB${ANSI_RESET}               : r220 version

This analysis contains ${ANSI_GREEN}GUNC${ANSI_RESET}, ${ANSI_GREEN}CheckM2${ANSI_RESET}, and ${ANSI_GREEN}GTDB${ANSI_RESET} information.

Required Parameters:
------------------
${ANSI_YELLOW}--metadata${ANSI_RESET}         : Path to metadata file
                     Required for ${ANSI_RED}COMPARATIVE_ANNOTATION${ANSI_RESET}

${ANSI_YELLOW}--accession_column${ANSI_RESET} : Column index in metadata matching genome filenames
                     (e.g., columns containing accession number of metadata)

Optional Parameters:
------------------
${ANSI_YELLOW}--pass_quality${ANSI_RESET}     : Quality filter for selecting genomes for GTDB-Tk analysis
                     Available options:
                     ${ANSI_GREEN}medium_quality.pass${ANSI_RESET}      : Completeness ≥ 50%, Contamination < 10%
                     ${ANSI_GREEN}high_quality.pass${ANSI_RESET}       : Completeness > 90%, Contamination < 5%
                     ${ANSI_GREEN}medium_quality_gunc.pass${ANSI_RESET} : medium_quality + GUNC pass
                     ${ANSI_GREEN}high_quality_gunc.pass${ANSI_RESET}  : high_quality + GUNC pass
                     ${ANSI_GREEN}QS50.pass${ANSI_RESET}               : QS50 score ≥ 50 (default)
                     ${ANSI_GREEN}QS50_gunc.pass${ANSI_RESET}          : QS50 pass + GUNC pass
                     ${ANSI_GREEN}all${ANSI_RESET}                     : No quality filtering

${ANSI_YELLOW}--inputDir${ANSI_RESET}        : Input directory containing genome bins (default: ${params.outdir_Base}/results/metagenome/ASSEMBLY_BINNING/final_bins)
                     Files must have .fa, .fna, or .fasta extension

${ANSI_YELLOW}--outdir_Base${ANSI_RESET}     : Base output directory (default: ${launchDir})
                     Final output will be in: ${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT

${ANSI_YELLOW}--cpus${ANSI_RESET}             : Number of CPUs to use (default: ${params.cpus})

Output Files:
-----------
- quality_taxonomy_combined_[TIMESTAMP].csv : Combined quality and taxonomy information
- bins_quality_passed/                     : Directory containing genomes passing quality filter
                                           These bins will be used in ${ANSI_RED}COMPARATIVE_ANNOTATION${ANSI_RESET}

Note: Selected quality filter (--pass_quality) determines which genomes will be 
      available for downstream analysis in ${ANSI_RED}COMPARATIVE_ANNOTATION${ANSI_RESET}
"""
exit 0
}


if (!new File(params.inputDir).exists()) {
    error """
${ANSI_RED}Input directory does not exist: ${params.inputDir}
Please specify a valid directory with --inputDir that contains genome bins.
Expected location: ${params.outdir_Base}/results/metagenome/ASSEMBLY_BINNING/final_bins
Or provide a custom path with --inputDir option.${ANSI_RESET}
"""
} else {
    if (new File(params.inputDir).list().length == 0) {
        error """
${ANSI_RED}Input directory is empty: ${params.inputDir}
Directory should contain genome bins with .fa, .fna, or .fasta extensions.
Please check the path or specify a different directory with --inputDir.${ANSI_RESET}
"""
    }
}

"""
if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir."
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir : assembled genomic fna files."
    }
}
"""

//${custom_logo}
//Visit documentation to find out description of this workflow: 
//https://metafun-doc.readthedocs.io/en/latest/workflows/BIN_ASSESSMENT.html
log.info """
=========================================================
${ANSI_GREEN}BIN_ASSESSMENT${ANSI_RESET} - Assess Genome Quality and Assign Taxonomy
=========================================================
(Apptainer Container with Nextflow)

Key Features:
  ✓ CheckM2-based quality assessment
  ✓ GUNC-based contamination check
  ✓ GTDB-Tk taxonomy classification
  ✓ Quality filtering for downstream analysis


This analysis contains ${ANSI_GREEN}GUNC${ANSI_RESET}, ${ANSI_GREEN}CheckM2${ANSI_RESET}, and ${ANSI_GREEN}GTDB${ANSI_RESET}.

\033[31mRequired Parameters:\033[0m 
------------------
${ANSI_RED}--metadata${ANSI_RESET}         : Path to metadata file
                     Required for genome metadata file generation. Metagenome metadata file is needed for automatic genome metadata file generation in this module.

${ANSI_RED}--accession_column${ANSI_RESET} : Column index in metadata matching prefix of genome filenames.
                     (e.g., columns containing accession number of metadata)


Optional Parameters:
------------------
${ANSI_YELLOW}--pass_quality${ANSI_RESET}     : Quality filter for selecting genomes for GTDB-Tk analysis
                     Available options:
                     ${ANSI_GREEN}medium_quality.pass${ANSI_RESET}      : Completeness ≥ 50%, Contamination < 10%
                     ${ANSI_GREEN}high_quality.pass${ANSI_RESET}       : Completeness > 90%, Contamination < 5%
                     ${ANSI_GREEN}medium_quality_gunc.pass${ANSI_RESET} : medium_quality + GUNC pass
                     ${ANSI_GREEN}high_quality_gunc.pass${ANSI_RESET}  : high_quality + GUNC pass
                     ${ANSI_GREEN}QS50.pass${ANSI_RESET}               : QS50 score ≥ 50 (default)
                     ${ANSI_GREEN}QS50_gunc.pass${ANSI_RESET}          : QS50 pass + GUNC pass
                     ${ANSI_GREEN}all${ANSI_RESET}                     : No quality filtering

${ANSI_YELLOW}--inputDir${ANSI_RESET}        : Input directory containing genome bins 
                    (default: ${params.outdir_Base}/results/metagenome/ASSEMBLY_BINNING/final_bins)
                     Files must have .fa, .fna, or .fasta extension

${ANSI_YELLOW}--outdir_Base${ANSI_RESET}     : Base output directory (default: ${launchDir})
                     Final output will be in: ${params.outdir_Base}/results/metagenome/BIN_ASSESSMENT

${ANSI_YELLOW}--cpus${ANSI_RESET}             : Number of CPUs to use (default: ${params.cpus})

${ANSI_YELLOW}--data_type${ANSI_RESET}       : Type of input data (default: 'MAG')

${ANSI_YELLOW}--run_id${ANSI_RESET}          : Custom run identifier (default: timestamp_workflowName)


\033[38;2;0;255;255mOutput Files:\033[0m
-----------
\033[38;2;0;255;255m- quality_taxonomy_combined_[TIMESTAMP].csv\033[0m : Combined quality and taxonomy information
\033[38;2;0;255;255m- bins_quality_passed/\033[0m                     : Directory containing genomes passing quality filter
                                           These bins will be used in \033[38;2;78;149;217mCOMPARATIVE_ANNOTATION${ANSI_RESET}
\033[38;2;0;255;255m- checkm2_[RUN_ID]/\033[0m : CheckM2 output directory
\033[38;2;0;255;255m- gunc_[RUN_ID]/\033[0m : GUNC output directory
\033[38;2;0;255;255m- gtdb_outdir_[RUN_ID]/\033[0m : GTDB-Tk output directory


Note: Selected quality filter (--pass_quality) determines which genomes will be 
      available for downstream analysis in \033[38;2;78;149;217mCOMPARATIVE_ANNOTATION${ANSI_RESET}.
      Using GENOME SELECTOR,  intersted species are selected and used for \033[38;2;78;149;217mCOMPARATIVE_ANNOTATION${ANSI_RESET}..

${params.metadata ? "" : ANSI_RED + "You did not provide any metadata. Please provide genome metadata using --metadata" + ANSI_RESET}
${params.metadata ? "" : ANSI_RED + "Any kind of genome metadata would be needed in COMPARATIVE_ANNOTATION analysis" + ANSI_RESET}
${params.metadata ? "" : ANSI_RED + "Created metadata file will contain only genome quality and taxonomy information." + ANSI_RESET}



Final metadata file is created and saved into your launch directory : 
${ANSI_BLUE}combined_metadata_quality_taxonomy_${params.run_id}.csv${ANSI_RESET}
"""
//:when you use --metadata option

process prepareInputFiles {
    publishDir "${params.outdir}/prepared_bins_${params.run_id}", mode: 'copy'
    
    input:
    path inputDir

    output:
    path "renamed_bins"

    script:
    """
    mkdir -p renamed_bins
    
    for file in ${inputDir}/*.{fa,fna,fasta,fa.gz,fna.gz,fasta.gz}; do
        if [[ -f "\$file" ]]; then
            basename=\$(basename "\$file")
            if [[ "\$file" == *.gz ]]; then
                gunzip -c "\$file" > "renamed_bins/\${basename%.gz}"
                basename=\${basename%.gz}
            else
                cp "\$file" "renamed_bins/\$basename"
            fi
            if [[ "\$basename" == *.fna || "\$basename" == *.fasta ]]; then
                mv "renamed_bins/\$basename" "renamed_bins/\${basename%.*}.fa"
            fi
        fi
    done
    """
}






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
    time checkm2 predict --threads ${task.cpus} --input ${bin_dir} \
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
    publishDir "${params.outdir}/checkm_gunc_combined_${params.run_id}", mode: 'copy'   , pattern: 'combined_report.tsv'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'bins_quality_passed'

    //publishDir "${params.outdir}", mode: 'copy'
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
    path "combined_report.tsv"// , emit : combined_report
    path "bins_quality_passed"

    script:
    def filter_command = params.pass_quality == 'all' ? 
        "cp \$(ls ${inputdir}/*.fa) bins_quality_passed/" :
        """
        col_num=\$(head -1 combined_report.tsv | tr '\\t' '\\n' | grep -n "^${params.pass_quality}\$" | cut -d: -f1)
        awk -v col="\$col_num" -F'\\t' '\$col=="True" {print \$1}' combined_report.tsv | sed -e 's/\$/.fa/g' > passlist
        cp \$(grep -Ff passlist <(ls ${inputdir}/*.fa)) bins_quality_passed/
        """
    
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate checkm2     
    cp ${checkm2Files}/quality_report.tsv ./ 
    cp ${guncFiles}/GUNC.progenomes_2.1.maxCSS_level.tsv ./ 
    python ${params.scripts_baseDir}/Checkm2_GUNC_combine_quality_pass.py quality_report.tsv GUNC.progenomes_2.1.maxCSS_level.tsv combined_report.tsv
    mkdir bins_quality_passed
    ${filter_command}
    """    
    /* """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate checkm2     
    cp ${checkm2Files}/quality_report.tsv ./ 
    cp ${guncFiles}/GUNC.progenomes_2.1.maxCSS_level.tsv ./ 
    python ${params.scripts_baseDir}/Checkm2_GUNC_combine_quality_pass.py quality_report.tsv GUNC.progenomes_2.1.maxCSS_level.tsv combined_report.tsv
    mkdir bins_quality_passed
    # get only at least medium quality genomes based on MIMAG and CheckM2
    col_num=\$(head -1 combined_report.tsv | awk -F'\\t' '{for(i=1;i<=NF;i++) if (\$i=="medium_quality.pass") print i}')
    awk -v col="\$col_num" -F'\\t' '\$col=="True" {print \$1}' combined_report.tsv | sed -e 's/\$/.fa/g' > passlist
    #cat passlist
    echo \$PWD
    cp \$(grep -Ff passlist <(ls ${inputdir}/* )) bins_quality_passed
    
    """*/
}
    //awk -F'\t' '\$NF=="True" {print \$1}' combined_report.tsv  | sed -e 's/\$/.fa/g' > passlist

process gtdbtk {
    publishDir "${params.outdir}/gtdb_outdir_${params.run_id}", mode: 'copy'
    //conda '$HOME/miniforge3/envs/gtdbtk-2.3.2'
    cpus params.cpus

    input:
    path combined_report
    path bins_quality_passed
    
    output:
    
    path "gtdb_outdir"
    path "quality_taxonomy_combined.csv"
    path "gtdbtk_r220.tsv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate gtdbtk-2.4.0
    export GTDBTK_DATA_PATH=${params.db_baseDir}/gtdb/release220/release220
    time gtdbtk classify_wf \
    --mash_db ${params.db_baseDir}/gtdb/release220/mash_220.db.msh --genome_dir ${bins_quality_passed}  -x fa  \
    --cpus ${task.cpus} --pplacer_cpus ${task.cpus} --out_dir gtdb_outdir
    awk -F'\\t' 'NR==1 || \$4=="True" {print \$0}'  ${combined_report} > quality_tmp.tsv
    if [[ -f gtdb_outdir/gtdbtk.bac120.summary.tsv ]]; then
        cat gtdb_outdir/gtdbtk.bac120.summary.tsv | grep classification | head -n 1 > gtdb_combined_tmp.tsv
    elif [[ -f gtdb_outdir/gtdbtk.ar53.summary.tsv ]]; then
        cat gtdb_outdir/gtdbtk.ar53.summary.tsv | grep classification | head -n 1 > gtdb_combined_tmp.tsv
    fi
    for file in gtdb_outdir/gtdbtk.*.summary.tsv; do
        if [[ -f "\$file" ]]; then
            cat "\$file" | tail -n+2  >> gtdb_combined_tmp.tsv
        fi
    done
    cut -f1,2  gtdb_combined_tmp.tsv > gtdb_classification_only.tsv
    mv gtdb_combined_tmp.tsv gtdbtk_r220.tsv

    micromamba activate checkm2 
    python ${params.scripts_baseDir}/GTDB_add2_check2gunc.py quality_tmp.tsv gtdb_classification_only.tsv quality_taxonomy_combined.csv
    """
}
//TBU :
//TBU : merge GTDB file together to output of GUNC +  Checkm2 
// cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep classification  | head -n 1  |cut -f1,2 > gtdb_combined_tmp.tsv
//     cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep -v classification |sort -u |cut -f1,2 >> gtdb_combined_tmp.tsv
//     cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep classification | head -n1  > gtdbtk_r214_${params.run_id}.tsv
//     cat gtdb_outdir_${params.run_id}/gtdbtk.{bac120,ar53}.summary.tsv | grep -v classification | sort -u >> gtdbtk_r214_${params.run_id}.tsv
    
    process createFinalFile {
    publishDir "${params.outdir}", mode: 'copy',  overwrite: true

    input:
    path A
    path quality_taxonomy_combined
    path C
    val existing_file

    output:
    path "quality_taxonomy_combined_final.csv"

    script:
    """
    if [ -f "${existing_file}" ]; then
        head -n 1 ${quality_taxonomy_combined} > quality_taxonomy_combined_final.csv
        (tail -n +2 "${existing_file}"; tail -n +2 ${quality_taxonomy_combined}) | sort | uniq >> quality_taxonomy_combined_final.csv
    else
        cp ${quality_taxonomy_combined} quality_taxonomy_combined_final.csv
    fi
    """
}

//   if [ ! -f "${launchDir}/quality_taxonomy_combined.csv" ]; then
//        head -n 1 quality_taxonomy_combined.csv > ${launchDir}/quality_taxonomy_combined.csv
//    fi
//    tail -n +2 quality_taxonomy_combined.csv >> ${launchDir}/quality_taxonomy_combined.csv

// automatic detecion of metadata file format
def detectFileFormat(file_path) {
    if (file_path.endsWith('.tsv')) {
        return 'tsv'
    } else if (file_path.endsWith('.csv')) {
        return 'csv'
    } else {
        // if no extension, check the first line
        def firstLine = new File(file_path).withReader { reader -> reader.readLine() }
        return firstLine.contains('\t') ? 'tsv' : 'csv'
    }
}


process combineMetadata {
    publishDir "${launchDir}", mode: 'copy'

    input:
    path quality_taxonomy_file
    path metadata_file

    output:
    path "combined_metadata_quality_taxonomy_${params.run_id}.csv"
    path "BIN_ASSESSMENT_metadata_merge_summary_${params.run_id}.log"

    when:
    params.metadata != null

    script:
    def format = detectFileFormat(metadata_file.toString())
    def separator = format == 'tsv' ? 'tab' : 'comma'
    """
   
    python ${params.scripts_baseDir}/combine_metadata_WMS_genome.py \
        -g ${quality_taxonomy_file} \
        -m ${metadata_file} \
        -o combined_metadata_quality_taxonomy_${params.run_id}.csv \
        -a ${params.accession_column}\
        -d ${params.data_type} > BIN_ASSESSMENT_metadata_merge_summary_${params.run_id}.log 2>& 1 

    lines=\$(wc -l < combined_metadata_quality_taxonomy_${params.run_id}.csv)
    if [ \$lines -le 1 ]; then
        echo "
    ╔═══════════════════════════════════ ERROR ═══════════════════════════════════╗
    ║                                                                             ║
    ║  No data in your specified metadata matched with genome names.              ║
    ║  Please check:                                                              ║
    ║  1. The file format (detected as ${format})                                 ║
    ║  2. The column index (--accession_column ${params.accession_column})        ║
    ║                                                                             ║
    ╚═════════════════════════════════════════════════════════════════════════════╝
    " >> BIN_ASSESSMENT_metadata_merge_summary_${params.run_id}.log
        cat BIN_ASSESSMENT_metadata_merge_summary_${params.run_id}.log
        exit 1
    fi

    """
}

process create_metadata_summary {
    publishDir "${launchDir}", mode: 'copy'

    input:
    path combined_metadata_file
    path log

    output:
    path "metadata_column_BIN_ASSESSMENT_summary.tsv"

    when:
    params.metadata != null

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv("${combined_metadata_file}", sep=',', header=0)

    summary = pd.DataFrame({
        'Column_Index': range(1,len(df.columns)+1),
        'Column_Name': df.columns
    })

    summary.to_csv('metadata_column_BIN_ASSESSMENT_summary.tsv', sep='\t', index=False)
    """
}




workflow {
    input_ch = Channel.fromPath(params.inputDir, type: 'dir', checkIfExists: true)
    if (params.metadata) {
        
        def metadata_file = file(params.metadata)
        if (!metadata_file.exists()) {
            error """
    ${ANSI_RED}ERROR: Metadata file does not exist: ${params.metadata}
    Please verify the file path is correct.${ANSI_RESET}
    """
        } else {
            log.info "Using metadata file: ${metadata_file.toAbsolutePath()}"
        }
    }


    prepared_bins = prepareInputFiles(input_ch)

    checkm2_result = runCheckM2(prepared_bins)
    gunc2_result = runGUNC(checkm2_result)
    qc_filter_result = combineFiles(checkm2_result, gunc2_result, prepared_bins)
    result_ch=gtdbtk(qc_filter_result)
    genome_metadata = createFinalFile(result_ch,"${params.outdir}/quality_taxonomy_combined_final.csv")

    if (params.metadata) {
        def metadata_file = file(params.metadata)
        if (!metadata_file.exists()) {
            error "ERROR: Metadata file does not exist: ${params.metadata}"
        } else {
            // Use the basename and stageInMode 'copy' to properly stage the file
            metadata= file(params.metadata)
            //metadata_ch = Channel.fromPath(metadata, checkIfExists: true)
      
            combined_metadata = combineMetadata(genome_metadata, metadata)
            create_metadata_summary(combined_metadata)
        }
    }

    // if (params.metadata) {
    //     //def metadata_file = file(params.metadata)
    //     metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
    //     //metadata_ch = Channel.fromPath(metadata_file, checkIfExists: true)
    //     combined_metadata = combineMetadata(genome_metadata, metadata_ch)
    //     create_metadata_summary(combined_metadata)
    // }    
}
