#!/usr/bin/env nextflow
// this script was optimized apptainer in any compurational environment. 
nextflow.enable.dsl=2 

params.db_baseDir = "/opt/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"


params.accessions = [] // List of accessions
params.cpus = 8 // Default CPU count

//params.outdir_Base = "/data2/leehg/OMD3/01pipeline/IBD_test/" // modify this if you want output files stored any other directory. 
params.inputDir = "${launchDir}/results/metagenome/RAWREAD_QC/read_filtered" // Input directory
params.outdir =  "${launchDir}/results/metagenome/ASSEMBLY_BINNING"
params.semibin2_mode = "self" 
params.megahit_presets = "default"
if (!new File(params.inputDir).exists()) {
    error "Input directory does not exist: ${params.inputDir}. Please specify a valid directory with --inputDir.\n Quality controlled paired-end fastq files are needed. "
} else {
    // Check if the directory is empty
    if (new File(params.inputDir).list().length == 0) {
        error "Input directory is empty: ${params.inputDir}. Please specify a directory with --inputDir . \n \n Quality controlled paired-end fastq files are needed"
    }
}

def validenvs = ['human_gut','dog_gut','ocean','soil','cat_gut','human_oral','mouse_gut','pig_gut','built_environment','wastewater','chicken_caecum','global']
def chosenFilter = validenvs.contains(params.semibin2_mode) ? params.semibin2_mode : 'self'
def ANSI_RESET = "\u001B[0m"
def ANSI_ORANGE = "\033[38;2;240;152;55m"
def ANSI_YELLOW = "\u001B[33m"
def ANSI_GREEN = "\u001B[32m"
def ANSI_CYAN = "\033[38;2;0;255;255m"
def ANSI_RED= "\u001B[31m"

log.info """\

===========================================
${ANSI_ORANGE}ASSEMBLY_BINNING${ANSI_RESET} - Metagenome Assembly & Binning
===========================================
(Apptainer Container with Nextflow)

Key Features:
  ✓ MEGAHIT assembly with multiple presets
  ✓ SemiBin2 and MetaBAT2 for single-sample binning
  ✓ DAS_Tool for bin refinement


This analysis contains ${ANSI_ORANGE}MEGAHIT${ANSI_RESET}, ${ANSI_ORANGE}MetaBAT2${ANSI_RESET },  ${ANSI_ORANGE}SemiBin2${ANSI_RESET} and ${ANSI_ORANGE}DAS_Tool${ANSI_RESET}.

${ANSI_RED}Required Parameters:${ANSI_RESET}
------------------
${ANSI_RED}--inputDir${ANSI_RESET}         : Directory containing quality-controlled reads
                     Should contain paired-end fastq files from RAWREAD_QC


Optional Parameters:
------------------
${ANSI_YELLOW}--megahit_presets${ANSI_RESET}  : Assembly presets for MEGAHIT
                     Available options:
                     ${ANSI_GREEN}default${ANSI_RESET}          : Default parameters (recommended for most data)
                     ${ANSI_GREEN}meta-large${ANSI_RESET}       : For large & complex metagenomes
                     ${ANSI_GREEN}meta-sensitive${ANSI_RESET}  : More sensitive but slower

${ANSI_YELLOW}--semibin2_mode${ANSI_RESET}    : self-supervised method or prebuilt environment SemiBin2 model 
                     Available models:
                     ${ANSI_GREEN}human_gut${ANSI_RESET}, ${ANSI_GREEN}dog_gut${ANSI_RESET}, ${ANSI_GREEN}ocean${ANSI_RESET}, 
                     ${ANSI_GREEN}soil${ANSI_RESET}, ${ANSI_GREEN}cat_gut${ANSI_RESET}, ${ANSI_GREEN}human_oral${ANSI_RESET},
                     ${ANSI_GREEN}mouse_gut${ANSI_RESET}, ${ANSI_GREEN}pig_gut${ANSI_RESET}, ${ANSI_GREEN}built_environment${ANSI_RESET},
                     ${ANSI_GREEN}wastewater${ANSI_RESET}, ${ANSI_GREEN}chicken_caecum${ANSI_RESET}, ${ANSI_GREEN}global${ANSI_RESET},
                     ${ANSI_GREEN}self${ANSI_RESET}
                     (default: self)

${ANSI_YELLOW}--cpus${ANSI_RESET}             : Number of CPUs to use (default: ${params.cpus})


${ANSI_CYAN}Output Directories:${ANSI_RESET}
------------------
${ANSI_CYAN}- assembled_contigs/${ANSI_RESET}      : Final assembled contigs
${ANSI_CYAN}- metabat2_bins/${ANSI_RESET}          : MetaBat2 binning results
${ANSI_CYAN}- semibin2_bins/${ANSI_RESET}          : SemiBin2 binning results
${ANSI_CYAN}- final_bins/${ANSI_RESET}             : Refined bins from DAS_Tool

${params.semibin2_mode == 'self' ? ANSI_YELLOW + 'NOTE: Running SemiBin2 in self-training mode' + ANSI_RESET : ''}

\033[38;2;180;180;180mNote:\033[0m Output bins will be used in: 
${ANSI_GREEN}BIN_ASSESSMENT${ANSI_RESET}
\033[38;2;78;149;217mCOMPARATIVE_ANNOTATION${ANSI_RESET}
    """
    .stripIndent(true)

// Process for assembly and renaming contigs
process AssemblyAndRename {
    tag "assembly megahit of ${accession}"
    publishDir "${params.outdir}/assembled_contigs", mode : 'copy'
    //maxForks params.accessions.size()
    cpus params.cpus
    errorStrategy 'Retry'
    //conda "$HOME/miniforge3/envs/ASSEMBLY_BINNING"
    input:
    //val accession 
    tuple val(accession), path(reads)

    output:
    tuple val(accession) , path("results/assembly/${accession}/${accession}_renamed_MH.contigs.fa")
    script:    
    if(params.megahit_presets == "default")
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    mkdir -p work/assembly
    megahit -1 ${reads[0]} -2 ${reads[1]} -o work/assembly/${accession}  -t ${task.cpus}    
    mkdir -p results/assembly/${accession}
    cat work/assembly/${accession}/final.contigs.fa | sed 's/^>/>${accession}_/g' | awk '{print \$1}' > results/assembly/${accession}/${accession}_renamed_MH.contigs.fa
    """
    else
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    mkdir -p work/assembly
    megahit -1 ${reads[0]} -2 ${reads[1]} -o work/assembly/${accession} --presets ${params.megahit_presets} -t ${task.cpus}    
    mkdir -p results/assembly/${accession} 
    cat work/assembly/${accession}/final.contigs.fa | sed 's/^>/>${accession}_/g' | awk '{print \$1}' > results/assembly/${accession}/${accession}_renamed_MH.contigs.fa
    """
    //replace to this after completion 
}
// Additional processes for Bowtie2 mapping, sorting, indexing, Metabat2, SemiBin2, and DAS Tool
// ...

process Bowtie2IndexBuild {
    tag "${accession}"
    //maxForks params.accessions.size()
    cpus 4
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 

    input : 
    tuple val(accession) , path(contigs)

    output : 
    tuple val(accession), path("${accession}_MH_bt2_index.*.bt2")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    bowtie2-build --threads ${task.cpus} ${contigs} ${accession}_MH_bt2_index
    """
    //publishDir "${params.outdir}/quant", mode:'copy'
}

process MHcontig2sortedbam {
    tag "${accession}"
    maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 
    cpus 4

    input: 
    tuple val(accession), path(reads) ,path(index)  // accession and read pair
      // bowtie2 accession and bowtie2 index files 
    // first is reac_ch second is from bowtie2index, acceiion2 is not used here
    output:
    tuple val(accession)  ,path("${accession}_sorted.bam") ,path("${accession}_sorted.bam.bai")
     //\$index_base

    script:
    def idxPrefix = index[0].toString().replaceAll(/\.(\d+)\.bt2$/, '')
    """
    echo ${idxPrefix}
    eval "\$(micromamba shell hook --shell bash)"
    #micromamba activate ASSEMBLY_BINNING
    micromamba activate RAWREAD_QC

    index_base=${accession}_MH_bt2_index
    bowtie2 --sensitive -p ${task.cpus} -x ${idxPrefix} -1 ${reads[0]} -2${reads[1]}  -S ${accession}.sam  
    samtools sort -@ ${task.cpus} -o ${accession}_sorted.bam ${accession}.sam 
    samtools index -@ ${task.cpus}  ${accession}_sorted.bam -o ${accession}_sorted.bam.bai
    """
// bowtie2 --sensitive -p ${params.cpus} -x ${accession}_MH_bt2_index -1 ${reads[0]} -2${reads[1]}  -S ${accession}.sam

}

process MB2_binning {
    tag "${accession}"
    maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 
    cpus 4
    publishDir "${params.outdir}/metabat2_bins", mode : 'copy'

    input: 
    //tuple val(accession),path(sorted_bam),path(bai)
    //tuple val(accession1), path(contig)
    tuple val(accession),path(sorted_bam),path(bai), path(contig)


    output:
    tuple val(accession), path ("${accession}_MB2*.fa"), optional: true, emit: mb2_output

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    jgi_summarize_bam_contig_depths --outputDepth ${accession}_jgi_summerize_depth.txt ${sorted_bam}
    metabat2 -i ${contig} -m 1500 -t ${task.cpus} -o ${accession}_MB2 -a ${accession}_jgi_summerize_depth.txt 
    """
}
//def validenvs = ['human_gut','dog_gut','ocean','soil','cat_gut','human_oral','mouse_gut','pig_gut','built_environment','wastewater','chicken_caecum','global']

//def chosenFilter = validenvs.contains(params.semibin2_mode) ? params.semibin2_mode : 'self'



process SB2_binning {
    tag "${accession}"
    maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 
    publishDir "${params.outdir}/semibin2_bins", mode : 'copy'
    cpus 4

    input: 
    //tuple val(accession),path(sorted_bam),path(bai)
    //tuple val(accession1), path(contig)
    tuple val(accession),path(sorted_bam),path(bai), path(contig)

    output:
    tuple val(accession), path ("*.fa"), optional: true, emit : sb2_output
    //path ("${accession}_SB2_[0-9]*.fa"), optional: true, emit : sb2_output
// semibin options will be added further . 
    script:
    if(chosenFilter == 'self')
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    SemiBin2  single_easy_bin  --self-supervised   \
    -i ${contig}  -b ${sorted_bam} -o ./ --compression none --tag-output ${accession}_${chosenFilter}_SB2 -p ${task.cpus}
    mv output_bins/*.fa ./ 

    """
    else
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    SemiBin2  single_easy_bin   --environment ${params.semibin2_mode}  \
    -i ${contig}  -b ${sorted_bam} -o ./ --compression none --tag-output ${accession}_${chosenFilter}_SB2 -p ${task.cpus}
    mv output_bins/*.fa ./ 
    """

}

process Contigs2bin_prep_mb2 {
    tag "${accession}"
    maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 
    cpus 1
    input:
    
    tuple val(accession), path(bins)

    output:
    tuple val(accession),path ("${accession}_mb2_aggregated_contigs2bin.tsv") ,optional: true


    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    Fasta_to_Contig2Bin.sh -e fa -i . > ${accession}_mb2_aggregated_contigs2bin.tsv
    """
}
process Contigs2bin_prep_sb2 {
    tag "${accession}"
    maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING'
    cpus 1 
 
    input:
    tuple val(accession), path(bins)


    output:
    tuple val(accession), path ("${accession}_sb2_aggregated_contigs2bin.tsv") ,  optional: true 
    
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    Fasta_to_Contig2Bin.sh -e fa -i . > ${accession}_sb2_aggregated_contigs2bin.tsv
    """
}

process Dastool { 
    publishDir "${params.outdir}/dastool_bins", mode : 'copy'
    tag "${accession}"
    //maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 
    cpus 4
        
        input:
    tuple val(accession),path(mb2) ,path(sb2),path(contig)
   
    //tuple val(accession), path(contig)

    output:
    path("*.fa"), optional: true
    
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate ASSEMBLY_BINNING
    DAS_Tool -i ${accession}_mb2_aggregated_contigs2bin.tsv,${accession}_sb2_aggregated_contigs2bin.tsv \
    -l metabat2,semibin2 -c ${contig} -o ${accession}_dastool -t ${task.cpus} --write_bins --score_threshold=0
    for i in ${accession}_dastool_DASTool_bins/*.fa; do 
        if [ -e \$i ]; then
            mv \$i ./
        fi
    done

    mv ${accession}_renamed_MH.contigs.fa ${accession}_dastool_DASTool_bins/
    """
}
// #mv ${accession}_dastool_DASTool_bins/*.fa ./ 
process get_bins {
    publishDir "${params.outdir}", mode : 'copy'
    tag "${accession}"
    maxForks params.accessions.size()
    //conda '$HOME/miniforge3/envs/ASSEMBLY_BINNING' 
            input:
    path(bins)
        output:
    path("final_bins/*.fa"), optional: true
       
        script:
    """
    
    mkdir final_bins
    mv *.fa final_bins
    """

}

workflow {
    // Define channels
    //read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )

    reads_ch = Channel.fromFilePairs("${params.inputDir}/*_{1,2}.{fastq,fq}{.gz,}", checkIfExists: true)
                .map { accession, reads -> 
                    // Split the accession name to remove the '_fastp' part
                    //def sam_accession = accession.split('_')[0]
                    def sam_accession = accession.replaceAll(/_fastp.*$/, '')
                    tuple(sam_accession, reads) 
                }
    reads_ch.view()

    // Process for assembly and renaming contigs
    assembly_and_rename = AssemblyAndRename(reads_ch)
    assembly_and_rename.view()
    bowtie2_index_build = Bowtie2IndexBuild(assembly_and_rename)
    
    acc_bam_bai =reads_ch\
        |join(bowtie2_index_build)\
        |MHcontig2sortedbam
    acc_bam_bai.view()

    acc_bam_bai.map { tuple(accession, sorted_bam, bai) -> tuple(accession, sorted_bam) }.view()
    // MB2 input accession, sorted_bam, MH_contig
            //MB2_input.view()
    //mb2_bins =  MB2_binning(accession_ch,sorted_bam_ch,contig_ch )
    //mb2_bins =  MB2_binning(acc_bam_bai,assembly_and_rename)
    mb2_bins = acc_bam_bai\
        |join(assembly_and_rename)\
        | MB2_binning

    mb2_bins.view()

    //sb2_bins =  SB2_binning(acc_bam_bai,assembly_and_rename)
    sb2_bins = acc_bam_bai\
        |join(assembly_and_rename)\
        | SB2_binning

    //sb2_bins.view()

        // Define a prefix for each use case
    def mb2_prefix = "mb2"
    def sb2_prefix = "sb2"
    mb2_bins.view()
    sb2_bins.view()
    mb2_contig2bin = Contigs2bin_prep_mb2(mb2_bins)
    sb2_contig2bin = Contigs2bin_prep_sb2(sb2_bins)
    dastool_bins = mb2_contig2bin\
        |join(sb2_contig2bin)\
        |join(assembly_and_rename)\
        |Dastool
    get_bins(dastool_bins)

    // if (mb2_count > 0 && sb2_count > 0) {
    //     mb2_contig2bin = Contigs2bin_prep_mb2(mb2_prefix, mb2_bins)
    //     sb2_contig2bin = Contigs2bin_prep_sb2(sb2_prefix, sb2_bins)
    //     dastool_bins = Dastool(mb2_contig2bin, sb2_contig2bin, assembly_and_rename)
    //     get_bins(dastool_bins)
    // } else {
    //     // If only one of MB2 or SB2 results exists, move it to final_bins directory
    //     if (mb2_count > 0) {
    //         get_bins(mb2_bins)
    //     }
    //     else (sb2_count> 0) {
    //         get_bins(sb2_bins)
    //     }
    //     }
    }    


    //nextflow run script7.nf -resume -with-report -with-trace -with-timeline -with-dag dag.png

    // Additional workflow logic for Bowtie2 mapping, sorting, indexing, Metabat2, SemiBin2, and DAS Tool


