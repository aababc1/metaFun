#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// =============================================================================
// WMS_STRAIN_apptainer.nf
// =============================================================================
// Strain-level Microbial Analysis Pipeline using InStrain
//
// Pipeline Overview:
//   1. Prevalence filter  → Select prevalent taxa from phyloseq object (sylph/kraken2)
//   2. Genome preparation → Fetch from GTDB, concatenate, build indices
//   3. Gene annotation    → Prodigal + eggNOG-mapper
//   4. Read mapping       → Bowtie2 alignment
//   5. InStrain analysis  → Profile + Compare
//   6. Data aggregation   → Shiny-ready RDS files
//
// Input: phyloseq object from WMS_TAXONOMY (sylph or kraken2 based)
//
// Author: MGSSB, Yonsei University
// =============================================================================

// ========================= Base Directories ==================================
params.db_baseDir = "/opt/database"
params.scripts_baseDir = "/scratch/tools/microbiome_analysis/scripts"

// ========================= Parameters ========================================

// Input/Output
params.input_dir = ""                    // Paired-end reads (*_{1,2}.fastq.gz)
params.outdir = "${launchDir}/results/metagenome/WMS_STRAIN"
params.metadata = ""                      // Sample metadata CSV
params.sampleIDcolumn = 1                 // Column number for sample IDs

// Taxonomy input (from WMS_TAXONOMY)
params.phyloseq_object = ""               // Phyloseq object RDS file (from WMS_TAXONOMY)
params.abundance_table = ""               // Alternative: merged abundance table (TSV)

// Prevalence filtering
params.prevalence_threshold = 5           // Minimum % of samples
params.min_abundance = 0.001              // Minimum relative abundance (0.1%)
params.gtdb_metadata = "${projectDir}/../scripts/gtdbr220_taxonomy_list_TaxID.tsv"
// genome_paths.tsv contains relative paths (filename + relative path)
params.gtdb_genome_paths_tsv = "${params.db_baseDir}/gtdb/release220/release220/skani/genome_paths.tsv"
// Base directory for genome files inside container
params.gtdb_genome_base = "${params.db_baseDir}/gtdb/release220/release220/skani"

// eggNOG-mapper
params.eggnog_db = "${params.db_baseDir}/eggNOG5"
// eggnog_cpus set below with auto-detection

// InStrain parameters
params.min_coverage = 5
params.min_freq = 0.05
params.min_read_ani = 0.92                // 0.92=strain, 0.95=species, 0.99=clonal
params.min_mapq = 2
params.fdr = 0.05                         // FDR for SNV calling
params.min_snp = 20                       // Minimum SNPs for genome-wide calculations
params.skip_mm = false                    // Skip mismatch profiling (faster)
params.skip_plot_generation = true        // Skip plot generation (faster)
params.database_mode = true               // Database mode for large references

// Resource - Auto-detect available CPUs
def availableCpus = Runtime.getRuntime().availableProcessors()
params.cpus = availableCpus                           // Max CPUs for single heavy tasks (bowtie2_build)
params.eggnog_cpus = availableCpus                    // Max CPUs for eggNOG-mapper (single heavy task)
params.mapping_cpus = Math.max(4, (availableCpus / 10) as int)   // CPUs per sample for bowtie2_mapping (default ~10 parallel jobs)
params.instrain_cpus = Math.max(4, (availableCpus / 10) as int)  // CPUs per sample for inStrain profile
params.instrain_compare_cpus = Math.max(8, (availableCpus / 2) as int)  // CPUs for inStrain compare (single job)

// Skip flags
params.skip_prevalence = false
params.skip_genome_prep = false
params.skip_annotation = false
params.skip_instrain = false

// Pre-existing files (skip stages)
params.reference_fasta = ""
params.stb_file = ""
params.gene_file = ""
params.protein_file = ""                  // Pre-existing protein file for eggNOG (genes.faa)
params.bowtie2_index = ""
params.bam_dir = ""
params.eggnog_annotations = ""
params.profile_dir = ""                   // Pre-existing InStrain profile directories (skip instrain_profile)

params.help = false

// ========================= Help ==============================================

def helpMessage() {
    def ANSI_RESET = "\u001B[0m"
    def ANSI_VIOLET = "\033[38;2;138;43;226m"
    def ANSI_YELLOW = "\u001B[33m"
    def ANSI_RED = "\u001B[31m"
    def ANSI_CYAN = "\033[38;2;0;255;255m"
    def ANSI_GREEN = "\u001B[32m"

    log.info """
${ANSI_VIOLET}╔═══════════════════════════════════════════════════════════════════════════╗
║     WMS_STRAIN - Strain-level Analysis with InStrain                        ║
╚═══════════════════════════════════════════════════════════════════════════╝${ANSI_RESET}

${ANSI_RED}USAGE:${ANSI_RESET}
    nextflow run WMS_STRAIN_apptainer.nf \\
        --input_dir /path/to/fastq \\
        --phyloseq_object /path/to/phyloseq_object.RDS \\
        --metadata /path/to/metadata.csv

${ANSI_RED}Required Parameters:${ANSI_RESET}
------------------
${ANSI_RED}--input_dir${ANSI_RESET}           Directory containing paired-end reads (*_{1,2}.fastq.gz)
                      OR use --bam_dir for pre-aligned BAM files
${ANSI_RED}--phyloseq_object${ANSI_RESET}     Phyloseq object RDS from WMS_TAXONOMY
                      (e.g., results/metagenome/WMS_TAXONOMY/phyloseq/phyloseq_object_sylph.RDS)

${ANSI_YELLOW}--metadata${ANSI_RESET}            (Optional) Sample metadata CSV file
                      If not provided, extracts from phyloseq sample_data()
                      If provided, replaces phyloseq metadata (must have same sample IDs)

${ANSI_YELLOW}Optional Parameters:${ANSI_RESET}
------------------
${ANSI_YELLOW}--abundance_table${ANSI_RESET}     Alternative: merged abundance TSV (instead of phyloseq)
${ANSI_YELLOW}--prevalence_threshold${ANSI_RESET} Minimum % of samples for prevalence (default: ${params.prevalence_threshold})
${ANSI_YELLOW}--min_abundance${ANSI_RESET}       Minimum relative abundance threshold (default: ${params.min_abundance})
${ANSI_YELLOW}--sampleIDcolumn${ANSI_RESET}      Column number for sample IDs in metadata (default: ${params.sampleIDcolumn})

${ANSI_YELLOW}Skip Options:${ANSI_RESET}
------------------
${ANSI_YELLOW}--skip_prevalence${ANSI_RESET}     Use existing prevalent taxa results
${ANSI_YELLOW}--skip_genome_prep${ANSI_RESET}    Use existing reference genome and index
${ANSI_YELLOW}--skip_annotation${ANSI_RESET}     Use existing gene annotations
${ANSI_YELLOW}--skip_instrain${ANSI_RESET}       Skip InStrain analysis

${ANSI_YELLOW}Resume from Intermediate Files:${ANSI_RESET}
------------------
${ANSI_YELLOW}--profile_dir${ANSI_RESET}         Use existing InStrain profile directories
                      (Skip instrain_profile, run compare and aggregate only)
                      Example: --profile_dir results/metagenome/WMS_STRAIN/05_instrain_profiles

${ANSI_CYAN}InStrain Parameters:${ANSI_RESET}
------------------
${ANSI_CYAN}--min_read_ani${ANSI_RESET}        ANI threshold (default: 0.92)
                      0.92 = strain-level, 0.95 = species-level, 0.99 = clonal
${ANSI_CYAN}--min_coverage${ANSI_RESET}        Minimum coverage (default: ${params.min_coverage})
${ANSI_CYAN}--min_freq${ANSI_RESET}            Minimum SNP frequency (default: ${params.min_freq})
${ANSI_CYAN}--database_mode${ANSI_RESET}       Database mode for large refs (default: ${params.database_mode})

${ANSI_CYAN}Resources:${ANSI_RESET}
------------------
${ANSI_CYAN}--cpus${ANSI_RESET}                Max CPUs for heavy single tasks: bowtie2_build, eggnog (default: ${params.cpus})
${ANSI_CYAN}--mapping_cpus${ANSI_RESET}        CPUs per sample for bowtie2 mapping (default: ${params.mapping_cpus})
${ANSI_CYAN}--instrain_cpus${ANSI_RESET}       CPUs per sample for inStrain profile (default: ${params.instrain_cpus})
${ANSI_CYAN}--instrain_compare_cpus${ANSI_RESET}  CPUs for inStrain compare (default: ${params.instrain_compare_cpus})

${ANSI_GREEN}Typical Workflow:${ANSI_RESET}
    RAWREAD_QC → WMS_TAXONOMY → ${ANSI_VIOLET}WMS_STRAIN${ANSI_RESET} → INTERACTIVE_STRAIN

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ========================= Logging ===========================================
def ANSI_RESET = "\u001B[0m"
def ANSI_VIOLET = "\033[38;2;138;43;226m"

log.info """
${ANSI_VIOLET}╔═══════════════════════════════════════════════════════════════════════════╗
║     WMS_STRAIN - Strain-level Microbial Analysis                            ║
╚═══════════════════════════════════════════════════════════════════════════╝${ANSI_RESET}

Configuration:
  Input reads    : ${params.input_dir ?: params.bam_dir ?: 'not specified'}
  Output         : ${params.outdir}
  Metadata       : ${params.metadata}
  Phyloseq       : ${params.phyloseq_object ?: 'not specified'}
  Abundance TSV  : ${params.abundance_table ?: 'not specified'}

Prevalence Filtering:
  Threshold      : ${params.prevalence_threshold}%
  Min abundance  : ${params.min_abundance}

InStrain Parameters:
  Min read ANI   : ${params.min_read_ani}
  Min coverage   : ${params.min_coverage}
  Min frequency  : ${params.min_freq}
  Database mode  : ${params.database_mode}

Resources:
  Max CPUs (build/eggnog) : ${params.cpus}
  Mapping CPUs per sample : ${params.mapping_cpus}
  InStrain CPUs per sample: ${params.instrain_cpus}
  InStrain Compare CPUs   : ${params.instrain_compare_cpus}
  eggNOG CPUs             : ${params.eggnog_cpus}

Skip Flags:
  Prevalence     : ${params.skip_prevalence}
  Genome prep    : ${params.skip_genome_prep}
  Annotation     : ${params.skip_annotation}
  InStrain       : ${params.skip_instrain}
"""

// =============================================================================
// STAGE 1: PREVALENCE FILTERING FROM PHYLOSEQ
// =============================================================================

process prevalence_filter_phyloseq {
    publishDir "${params.outdir}/01_prevalent_taxa", mode: 'copy'
    tag "prevalence_filter"

    input:
    path(phyloseq_rds)
    path(external_metadata)  // Optional: external metadata file or 'NO_METADATA'

    output:
    path("prevalent_taxa_taxonomy_ids.txt"), emit: taxonomy_ids
    path("prevalent_taxa_metadata.tsv"), emit: taxa_metadata
    path("prevalent_taxa_genome_paths.txt"), emit: genome_paths
    path("prevalence_summary.tsv"), emit: summary
    path("sample_metadata.csv"), emit: sample_metadata  // Extracted or provided metadata
    path("prevalence_filter.log"), emit: log

    when:
    !params.skip_prevalence

    script:
    def use_external_meta = external_metadata.name != 'NO_METADATA'
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(phyloseq)
        library(dplyr)
        library(data.table)
    })

    log_file <- file("prevalence_filter.log", open = "wt")
    sink(log_file, type = "message")
    sink(log_file, type = "output", append = TRUE)

    cat("=== Prevalence Filtering from Phyloseq ===\\n")
    cat("Threshold: ${params.prevalence_threshold}%\\n")
    cat("Min abundance: ${params.min_abundance}\\n\\n")

    # Load phyloseq object
    ps <- readRDS("${phyloseq_rds}")
    cat("Loaded phyloseq object\\n")
    cat("  Samples:", nsamples(ps), "\\n")
    cat("  Taxa:", ntaxa(ps), "\\n\\n")

    # Extract or use external metadata
    use_external <- ${use_external_meta ? 'TRUE' : 'FALSE'}

    if (use_external) {
        cat("Using external metadata file: ${external_metadata}\\n")
        sample_meta <- read.csv("${external_metadata}", stringsAsFactors = FALSE)
    } else {
        cat("Extracting metadata from phyloseq sample_data()\\n")
        if (!is.null(sample_data(ps))) {
            sample_meta <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
            sample_meta <- data.frame(sample_meta, stringsAsFactors = FALSE)  # Ensure plain data.frame
            sample_meta\$sample_id <- rownames(as.data.frame(sample_data(ps)))
            # Move sample_id to first column
            sample_meta <- sample_meta[, c("sample_id", setdiff(names(sample_meta), "sample_id")), drop = FALSE]
        } else {
            # Create minimal metadata from sample names
            sample_meta <- data.frame(sample_id = sample_names(ps), stringsAsFactors = FALSE)
        }
    }

    cat("Metadata samples:", nrow(sample_meta), "\\n")
    cat("Metadata columns:", paste(colnames(sample_meta), collapse = ", "), "\\n\\n")

    # Save sample metadata
    write.csv(sample_meta, "sample_metadata.csv", row.names = FALSE)

    # Get OTU table and calculate prevalence
    otu <- as.data.frame(otu_table(ps))
    if (!taxa_are_rows(ps)) {
        otu <- t(otu)
    }

    # Calculate prevalence (% of samples where taxon is present)
    prevalence_threshold <- ${params.prevalence_threshold} / 100
    min_abundance <- ${params.min_abundance}

    # For each taxon, count samples with abundance > min_abundance
    prevalence_df <- data.frame(
        taxon = rownames(otu),
        prevalence = apply(otu, 1, function(x) sum(x > min_abundance) / length(x)),
        mean_abundance = apply(otu, 1, mean),
        max_abundance = apply(otu, 1, max)
    )

    # Filter by prevalence threshold
    prevalent_taxa <- prevalence_df %>%
        filter(prevalence >= prevalence_threshold) %>%
        arrange(desc(prevalence))

    cat("Prevalent taxa (>=", ${params.prevalence_threshold}, "% samples):", nrow(prevalent_taxa), "\\n\\n")

    # Get taxonomy information
    if (!is.null(tax_table(ps))) {
        tax_df <- as.data.frame(tax_table(ps))
        tax_df\$taxon <- rownames(tax_df)
        prevalent_taxa <- left_join(prevalent_taxa, tax_df, by = "taxon")
    }

    # Save summary
    write.table(prevalent_taxa, "prevalence_summary.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    # Extract taxonomy IDs (for GTDB matching)
    write.table(prevalent_taxa\$taxon, "prevalent_taxa_taxonomy_ids.txt",
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Load GTDB metadata and match
    gtdb_meta <- fread("/scratch/tools/microbiome_analysis/scripts/gtdbr220_taxonomy_list_TaxID.tsv", sep = "\\t")
    cat("Loaded GTDB metadata:", nrow(gtdb_meta), "entries\\n")

    # Match prevalent taxa to GTDB accessions
    matched_taxa <- data.frame()

    if ("taxonomy_id" %in% colnames(gtdb_meta)) {
        matched_taxa <- gtdb_meta %>%
            filter(taxonomy_id %in% prevalent_taxa\$taxon)
    }

    if (nrow(matched_taxa) == 0 && "Species" %in% colnames(prevalent_taxa)) {
        matched_taxa <- gtdb_meta %>%
            filter(grepl(paste(prevalent_taxa\$Species, collapse = "|"), gtdb_taxonomy))
    }

    cat("Matched GTDB entries:", nrow(matched_taxa), "\\n")

    # Save matched metadata
    if (nrow(matched_taxa) > 0) {
        write.table(matched_taxa, "prevalent_taxa_metadata.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
    } else {
        write.table(data.frame(accession = character(), gtdb_taxonomy = character(), taxonomy_id = character()),
                    "prevalent_taxa_metadata.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
    }

    # Get genome paths from GTDB genome_paths.tsv (contains relative paths)
    genome_paths_tsv <- "${params.gtdb_genome_paths_tsv}"
    genome_base <- "${params.gtdb_genome_base}"

    if (file.exists(genome_paths_tsv)) {
        # genome_paths.tsv format: "filename relative_path/" (space separated, may have trailing space)
        # e.g., "GCA_000008085.1_genomic.fna.gz database/GCA/000/008/085/ "
        genome_df <- fread(genome_paths_tsv, header = FALSE, fill = TRUE)
        # Only keep first two columns (filename, relative_path)
        genome_df <- genome_df[, 1:2]
        colnames(genome_df) <- c("filename", "relative_path")

        # Extract accession from filename (e.g., GCA_000008085.1 from GCA_000008085.1_genomic.fna.gz)
        genome_df\$accession <- gsub("_genomic\\\\.fna(\\\\.gz)?\$", "", genome_df\$filename)

        cat("Loaded genome_paths.tsv:", nrow(genome_df), "entries\\n")

        if (nrow(matched_taxa) > 0) {
            # Match accessions
            matched_genomes <- genome_df %>%
                filter(accession %in% matched_taxa\$accession)

            # Build full container paths
            full_paths <- paste0(genome_base, "/", matched_genomes\$relative_path, matched_genomes\$filename)

            writeLines(full_paths, "prevalent_taxa_genome_paths.txt")
            cat("Genome paths found:", length(full_paths), "\\n")
        } else {
            file.create("prevalent_taxa_genome_paths.txt")
            cat("No genome paths found (no matched taxa)\\n")
        }
    } else {
        file.create("prevalent_taxa_genome_paths.txt")
        cat("Warning: Genome paths TSV not found:", genome_paths_tsv, "\\n")
    }

    sink()
    cat("=== Prevalence filtering complete ===\\n")
    """
}

// Alternative: Filter from abundance table (TSV) - for backward compatibility
process prevalence_filter_table {
    publishDir "${params.outdir}/01_prevalent_taxa", mode: 'copy'
    tag "prevalence_filter"

    input:
    path(abundance_table)

    output:
    path("prevalent_taxa_taxonomy_ids.txt"), emit: taxonomy_ids
    path("prevalent_taxa_metadata.tsv"), emit: taxa_metadata
    path("prevalent_taxa_genome_paths.txt"), emit: genome_paths
    path("prevalence_filter.log"), emit: log

    when:
    !params.skip_prevalence

    script:
    """
    #!/bin/bash
    set -e

    PREVALENCE_DECIMAL=\$(echo "scale=4; ${params.prevalence_threshold} / 100" | bc)

    echo "=== Prevalence Filtering ===" | tee prevalence_filter.log
    echo "Threshold: ${params.prevalence_threshold}%" | tee -a prevalence_filter.log

    # Extract prevalent taxa from abundance table
    awk -v threshold="\$PREVALENCE_DECIMAL" '
    NR==1 {
        total_samples = NF - 1
        min_samples = total_samples * threshold
        next
    }
    {
        non_zero_count = 0
        for(i=2; i<=NF; i++) {
            if(\$i > 0) non_zero_count++
        }
        if(non_zero_count >= min_samples) print \$1
    }' ${abundance_table} > prevalent_taxa_taxonomy_ids.txt

    echo "Prevalent taxa: \$(wc -l < prevalent_taxa_taxonomy_ids.txt)" | tee -a prevalence_filter.log

    # Extract metadata from GTDB
    awk '
    BEGIN { FS="\\t" }
    NR==FNR { prevalent_ids[\$1] = 1; next }
    FNR==1 { print "accession\\tgtdb_taxonomy\\ttaxonomy_id"; next }
    { if (\$3 in prevalent_ids) print \$1 "\\t" \$2 "\\t" \$3 }
    ' prevalent_taxa_taxonomy_ids.txt ${params.gtdb_metadata} > prevalent_taxa_metadata.tsv

    # Extract genome paths using genome_paths.tsv (relative paths format)
    # genome_paths.tsv format: "filename relative_path/"
    GENOME_PATHS_TSV="${params.gtdb_genome_paths_tsv}"
    GENOME_BASE="${params.gtdb_genome_base}"

    # Simple approach: extract accession from metadata, match in genome_paths.tsv
    > prevalent_taxa_genome_paths.txt
    while IFS=\$'\\t' read -r acc rest; do
        [[ "\$acc" == "accession" ]] && continue
        # Find matching line in genome_paths.tsv
        grep "^\${acc}_genomic" "\$GENOME_PATHS_TSV" | while read filename relpath extra; do
            echo "\${GENOME_BASE}/\${relpath}\${filename}"
        done
    done < prevalent_taxa_metadata.tsv >> prevalent_taxa_genome_paths.txt

    echo "Genome paths: \$(wc -l < prevalent_taxa_genome_paths.txt)" | tee -a prevalence_filter.log
    """
}

// =============================================================================
// STAGE 2: GENOME PREPARATION
// =============================================================================

process concat_genomes {
    publishDir "${params.outdir}/02_genome_prep", mode: 'copy'
    tag "concat_genomes"

    input:
    path(genome_paths)

    output:
    path("all_genomes_combined.fa"), emit: combined_fasta
    path("genome_list.txt"), emit: genome_list

    when:
    !params.skip_genome_prep && params.reference_fasta.isEmpty()

    script:
    """
    #!/bin/bash
    set -e
    > genome_list.txt

    while read genome_path; do
        if [[ -f "\$genome_path" ]]; then
            if [[ "\$genome_path" == *.gz ]]; then
                zcat "\$genome_path" >> all_genomes_combined.fa
            else
                cat "\$genome_path" >> all_genomes_combined.fa
            fi
            basename "\$genome_path" >> genome_list.txt
        fi
    done < ${genome_paths}

    echo "Combined \$(wc -l < genome_list.txt) genomes"
    """
}

process generate_stb {
    publishDir "${params.outdir}/02_genome_prep", mode: 'copy'
    tag "generate_stb"

    input:
    path(genome_paths)

    output:
    path("prevalent_taxa.stb"), emit: stb_file

    when:
    !params.skip_genome_prep && params.stb_file.isEmpty()

    script:
    """
    #!/usr/bin/env python3
    import os, gzip

    stb = {}
    with open("${genome_paths}", 'r') as f:
        for line in f:
            genome_path = line.strip()
            if not genome_path or not os.path.exists(genome_path):
                continue
            genome_name = os.path.basename(genome_path)
            opener = gzip.open if genome_path.endswith('.gz') else open
            mode = 'rt' if genome_path.endswith('.gz') else 'r'
            with opener(genome_path, mode) as handle:
                for fasta_line in handle:
                    if fasta_line.startswith('>'):
                        scaffold_id = fasta_line[1:].strip().split()[0]
                        stb[scaffold_id] = genome_name

    with open('prevalent_taxa.stb', 'w') as out:
        for scaffold, genome in sorted(stb.items(), key=lambda x: x[1]):
            out.write(f"{scaffold}\\t{genome}\\n")
    print(f"Generated STB with {len(stb)} scaffolds")
    """
}

process bowtie2_build {
    publishDir "${params.outdir}/02_genome_prep/bowtie2_index", mode: 'copy'
    tag "bowtie2_build"
    cpus params.cpus

    input:
    path(reference_fasta)

    output:
    path("reference_index*"), emit: index

    when:
    !params.skip_genome_prep && params.bowtie2_index.isEmpty()

    script:
    """
    bowtie2-build --threads ${task.cpus} ${reference_fasta} reference_index
    """
}

// =============================================================================
// STAGE 3: GENE ANNOTATION (Parallelized per-genome prodigal)
// =============================================================================

// Run prodigal on individual genomes in parallel - much faster than on concatenated file
process prodigal_per_genome {
    tag "${genome_id}"
    cpus 1

    input:
    tuple val(genome_id), val(genome_path)

    output:
    tuple val(genome_id), path("${genome_id}.faa"), path("${genome_id}.fna"), path("${genome_id}.gff")

    when:
    !params.skip_annotation && params.gene_file.isEmpty()

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate base

    # genome_path is container path (e.g., /opt/database/gtdb/...)
    # Decompress if needed
    if [[ "${genome_path}" == *.gz ]]; then
        zcat "${genome_path}" > genome_temp.fa
        INPUT_FILE=genome_temp.fa
    else
        INPUT_FILE="${genome_path}"
    fi

    # Run prodigal on single genome
    prodigal -i \$INPUT_FILE -a "${genome_id}.faa" -d "${genome_id}.fna" -f gff -o "${genome_id}.gff" -p single
    """
}

// Concatenate all per-genome prodigal results
process concat_gene_predictions {
    publishDir "${params.outdir}/03_gene_annotation", mode: 'copy'
    tag "concat_genes"

    input:
    path(faa_files)
    path(fna_files)
    path(gff_files)

    output:
    path("genes.faa"), emit: proteins
    path("genes.fna"), emit: genes
    path("genes.gff"), emit: gff

    script:
    """
    # Concatenate all protein sequences
    cat *.faa > genes.faa

    # Concatenate all gene sequences
    cat *.fna > genes.fna

    # Concatenate all GFF files (handle headers)
    # First file with header, rest without header lines
    first=true
    for gff in *.gff; do
        if [ "\$first" = true ]; then
            cat "\$gff" >> genes.gff
            first=false
        else
            grep -v "^#" "\$gff" >> genes.gff || true
        fi
    done

    echo "Combined \$(ls *.faa | wc -l) genome annotations"
    echo "Total proteins: \$(grep -c '^>' genes.faa)"
    echo "Total genes: \$(grep -c '^>' genes.fna)"
    """
}

// Legacy single-file prodigal (kept for backward compatibility if user provides pre-combined fasta)
process prodigal {
    publishDir "${params.outdir}/03_gene_annotation", mode: 'copy'
    tag "prodigal"

    input:
    path(reference_fasta)

    output:
    path("genes.faa"), emit: proteins
    path("genes.fna"), emit: genes
    path("genes.gff"), emit: gff

    when:
    !params.skip_annotation && params.gene_file.isEmpty()

    script:
    """
    prodigal -i ${reference_fasta} -a genes.faa -d genes.fna -f gff -o genes.gff -p meta
    """
}

process eggnog_mapper {
    publishDir "${params.outdir}/03_gene_annotation", mode: 'copy'
    tag "eggnog_mapper"
    cpus params.eggnog_cpus

    input:
    path(proteins)

    output:
    path("eggnog_results.emapper.annotations"), emit: annotations

    when:
    !params.skip_annotation && params.eggnog_annotations.isEmpty()

    script:
    """
    emapper.py -i ${proteins} --output eggnog_results \\
        --data_dir ${params.eggnog_db} --cpu ${task.cpus} -m mmseqs --override
    """
}

// =============================================================================
// STAGE 4: READ MAPPING
// =============================================================================

process bowtie2_mapping {
    publishDir "${params.outdir}/04_bam_files", mode: 'copy'
    tag "mapping ${accession}"
    cpus params.mapping_cpus

    input:
    tuple val(accession), path(reads)
    path(bowtie2_index)
    path(reference_fasta)

    output:
    tuple val(accession), path("${accession}.sorted.bam"), path("${accession}.sorted.bam.bai")

    script:
    // Use mapping_cpus (default 8) for per-sample parallel mapping
    // samtools uses half the cpus for sorting (4 for sorting is optimal)
    def samtools_threads = Math.max(1, (task.cpus / 2) as int)
    """
    bowtie2 -p ${task.cpus} -x reference_index \\
        -1 ${reads[0]} -2 ${reads[1]} --sensitive-local --no-unal \\
        2> ${accession}_bowtie2.log \\
        | samtools sort -@ ${samtools_threads} -o ${accession}.sorted.bam -
    samtools index -@ ${samtools_threads} ${accession}.sorted.bam
    """
}

// =============================================================================
// STAGE 5: INSTRAIN PROFILE
// =============================================================================

process instrain_profile {
    publishDir "${params.outdir}/05_instrain_profiles", mode: 'copy'
    tag "instrain ${accession}"
    cpus params.instrain_cpus

    input:
    tuple val(accession), path(bam_file), path(bam_index)
    path(reference_fasta)
    path(stb_file)
    path(gene_file)

    output:
    tuple val(accession), path("${accession}_instrain_profile")

    when:
    !params.skip_instrain

    script:
    def stb_arg = stb_file.name != 'NO_STB' ? "-s ${stb_file}" : ""
    def gene_arg = gene_file.name != 'NO_GENES' ? "-g ${gene_file}" : ""
    def skip_mm_arg = params.skip_mm ? "--skip_mm" : ""
    def skip_plot_arg = params.skip_plot_generation ? "--skip_plot_generation" : ""
    def db_mode_arg = params.database_mode ? "--database_mode" : ""
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate instrain_env

    inStrain profile ${bam_file} ${reference_fasta} \\
        -o ${accession}_instrain_profile \\
        -p ${task.cpus} \\
        -c ${params.min_coverage} \\
        -f ${params.min_freq} \\
        --min_read_ani ${params.min_read_ani} \\
        --min_mapq ${params.min_mapq} \\
        --fdr ${params.fdr} \\
        --min_snp ${params.min_snp} \\
        ${stb_arg} ${gene_arg} ${skip_mm_arg} ${skip_plot_arg} ${db_mode_arg}
    """
}

// =============================================================================
// STAGE 5.5: VALIDATE INSTRAIN PROFILES
// =============================================================================

process validate_instrain_profiles {
    publishDir "${params.outdir}/05_instrain_profiles", mode: 'copy', pattern: "*.txt"
    tag "validate_profiles"

    input:
    path(profile_dirs)

    output:
    path("valid_profiles/*"), emit: valid_dirs, optional: true
    path("valid_profiles.txt"), emit: valid_list
    path("invalid_profiles.txt"), emit: invalid_list
    path("validation_summary.txt"), emit: summary

    script:
    """
    #!/bin/bash
    mkdir -p valid_profiles
    > valid_profiles.txt
    > invalid_profiles.txt

    echo "=== InStrain Profile Validation ===" > validation_summary.txt
    echo "Date: \$(date)" >> validation_summary.txt
    echo "" >> validation_summary.txt

    for profile_dir in ${profile_dirs.join(' ')}; do
        if [ -d "\$profile_dir" ]; then
            sample_name=\$(basename "\$profile_dir")
            genome_info=\$(find -L "\$profile_dir" -name "*genome_info.tsv" 2>/dev/null | head -1)

            if [ -f "\$genome_info" ]; then
                genome_count=\$(((\$(wc -l < "\$genome_info") - 1)))
                if [ \$genome_count -gt 0 ]; then
                    # Copy valid profile to output
                    cp -rL "\$profile_dir" valid_profiles/
                    echo "\$sample_name" >> valid_profiles.txt
                    echo "[VALID] \$sample_name: \$genome_count genomes detected" >> validation_summary.txt
                else
                    echo "\$sample_name" >> invalid_profiles.txt
                    echo "[INVALID] \$sample_name: 0 genomes detected (low coverage)" >> validation_summary.txt
                fi
            else
                echo "\$sample_name" >> invalid_profiles.txt
                echo "[INVALID] \$sample_name: No genome_info.tsv found" >> validation_summary.txt
            fi
        fi
    done

    valid_count=\$(wc -l < valid_profiles.txt)
    invalid_count=\$(wc -l < invalid_profiles.txt)

    echo "" >> validation_summary.txt
    echo "=== Summary ===" >> validation_summary.txt
    echo "Valid profiles: \$valid_count" >> validation_summary.txt
    echo "Invalid profiles: \$invalid_count" >> validation_summary.txt

    if [ \$valid_count -lt 2 ]; then
        echo "" >> validation_summary.txt
        echo "WARNING: Less than 2 valid profiles. inStrain compare requires at least 2 profiles." >> validation_summary.txt
    fi

    cat validation_summary.txt
    """
}

// =============================================================================
// STAGE 6: INSTRAIN COMPARE
// =============================================================================

process instrain_compare {
    publishDir "${params.outdir}/06_instrain_compare", mode: 'copy'
    tag "instrain_compare"
    cpus params.instrain_compare_cpus  // Single job - use more CPUs

    input:
    path(profile_dirs)
    path(stb_file)

    output:
    path("instrainComparer_output")

    when:
    !params.skip_instrain

    script:
    def stb_arg = stb_file.name != 'NO_STB' ? "-s ${stb_file}" : ""
    def db_mode_arg = params.database_mode ? "--database_mode" : ""
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate instrain_env

    inStrain compare -i ${profile_dirs.join(' ')} \\
        -o instrainComparer_output -p ${task.cpus} ${stb_arg} ${db_mode_arg}
    """
}

// =============================================================================
// STAGE 7: DATA AGGREGATION FOR SHINY
// =============================================================================

process aggregate_for_shiny {
    publishDir "${params.outdir}/07_shiny_data", mode: 'copy'
    tag "aggregate_shiny"

    input:
    path(profile_dirs)
    path(compare_dir)
    path(metadata)
    path(stb_file)
    path(eggnog_annotations)
    path(gtdb_metadata)

    output:
    path("integrated_microbiome_data.rds")
    path("pN_pS_gene_level.rds")
    path("pN_pS_genome_wide.rds")
    path("eggnog_annotations_subset.rds")

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr); library(tidyr); library(data.table)

    cat("=== Aggregating for Shiny ===\\n")

    metadata <- read.csv("${metadata}", stringsAsFactors = FALSE)
    sample_id_col <- colnames(metadata)[${params.sampleIDcolumn}]
    cat("Sample ID column:", sample_id_col, "\\n")

    gtdb_data <- fread("${gtdb_metadata}", sep = "\\t") %>% as.data.frame()
    stb_data <- fread("${stb_file}", header = FALSE, col.names = c("scaffold", "genome")) %>% as.data.frame()

    # Combine genome_info
    genome_files <- list.files(".", "genome_info\\\\.tsv\$", recursive = TRUE, full.names = TRUE)
    cat("Found", length(genome_files), "genome_info files\\n")
    genome_list <- lapply(genome_files, function(f) {
        df <- fread(f, showProgress = FALSE)
        df\$sample_id <- gsub("_instrain_profile", "", basename(dirname(dirname(f))))
        df
    })
    genome_combined <- rbindlist(genome_list, fill = TRUE) %>% as.data.frame()
    cat("Combined genome_info rows:", nrow(genome_combined), "\\n")

    # Combine gene_info.tsv files (filename includes sample prefix: {sample}_instrain_profile_gene_info.tsv)
    gene_files <- list.files(".", "_gene_info\\\\.tsv\$", recursive = TRUE, full.names = TRUE)
    cat("Found", length(gene_files), "gene_info files\\n")

    # Also find genes_SNP_count.csv files (for N_sites, S_sites - required for pN/pS)
    snp_files <- c(
        list.files(".", "genes_SNP_count\\\\.csv\$", recursive = TRUE, full.names = TRUE),
        list.files(".", "genes_SNP_count\\\\.csv\\\\.gz\$", recursive = TRUE, full.names = TRUE)
    )
    cat("Found", length(snp_files), "genes_SNP_count files\\n")

    # File path: {sample}_instrain_profile/output/{sample}_instrain_profile_gene_info.tsv
    # dirname(f) = {sample}_instrain_profile/output
    # dirname(dirname(f)) = {sample}_instrain_profile
    if (length(gene_files) > 0) {
        cat("First gene_info file path:", gene_files[1], "\\n")
        gene_list <- lapply(gene_files, function(f) {
            df <- fread(f, showProgress = FALSE)
            df\$sample_id <- gsub("_instrain_profile", "", basename(dirname(dirname(f))))
            df
        })
        gene_combined <- rbindlist(gene_list, fill = TRUE) %>% as.data.frame()
        cat("Combined gene_info rows:", nrow(gene_combined), "\\n")
        cat("Gene_info sample_ids:", paste(unique(gene_combined\$sample_id), collapse=", "), "\\n")
    } else {
        gene_combined <- data.frame()
    }

    # Read genes_SNP_count files for N_sites and S_sites
    # File path: {sample}_instrain_profile/raw_data/genes_SNP_count.csv.gz
    # dirname(f) = {sample}_instrain_profile/raw_data
    # dirname(dirname(f)) = {sample}_instrain_profile
    if (length(snp_files) > 0) {
        cat("First genes_SNP_count file path:", snp_files[1], "\\n")
        snp_list <- lapply(snp_files, function(f) {
            df <- fread(f, showProgress = FALSE)
            df\$sample_id <- gsub("_instrain_profile", "", basename(dirname(dirname(f))))
            # Select only columns needed for pN/pS
            cols_needed <- c("gene", "SNV_N_count", "SNV_S_count", "N_sites", "S_sites", "sample_id")
            cols_available <- intersect(cols_needed, colnames(df))
            df[, ..cols_available]
        })
        snp_combined <- rbindlist(snp_list, fill = TRUE) %>% as.data.frame()
        cat("Combined genes_SNP_count rows:", nrow(snp_combined), "\\n")
        cat("SNP_count sample_ids:", paste(unique(snp_combined\$sample_id), collapse=", "), "\\n")
    } else {
        snp_combined <- data.frame()
    }

    # Gene-level pN/pS - join gene_info with genes_SNP_count
    if (nrow(gene_combined) > 0 && nrow(snp_combined) > 0) {
        cat("Joining gene_info with genes_SNP_count for pN/pS calculation...\\n")

        # Join by gene and sample_id
        gene_pnps <- gene_combined %>%
            select(gene, scaffold, sample_id, coverage, breadth, breadth_minCov) %>%
            inner_join(snp_combined, by = c("gene", "sample_id")) %>%
            inner_join(stb_data, by = "scaffold") %>%
            filter(!is.na(N_sites) & !is.na(S_sites)) %>%
            mutate(
                clean_genome = gsub("_genomic\\\\.fna(\\\\.gz)?\$", "", genome),
                pN = ifelse(N_sites > 0, SNV_N_count / N_sites, NA),
                pS = ifelse(S_sites > 0, SNV_S_count / S_sites, NA),
                pNpS = ifelse(!is.na(pN) & !is.na(pS) & pS > 0, pN / pS, NA)
            )
        cat("Gene pN/pS rows after join:", nrow(gene_pnps), "\\n")
        cat("pN/pS range:", round(min(gene_pnps\$pNpS, na.rm = TRUE), 4), "to",
            round(max(gene_pnps\$pNpS, na.rm = TRUE), 4), "\\n")

        # Join sample metadata to gene_pnps for Shiny visualization
        cat("Joining sample metadata to gene pN/pS data...\\n")
        cat("Sample ID column:", sample_id_col, "\\n")
        cat("Metadata sample IDs:", paste(head(metadata[[sample_id_col]], 3), collapse=", "), "\\n")
        cat("Gene data sample IDs:", paste(head(unique(gene_pnps\$sample_id), 3), collapse=", "), "\\n")

        # Rename sample_id_col to match gene_pnps sample_id
        metadata_for_join <- metadata
        colnames(metadata_for_join)[which(colnames(metadata_for_join) == sample_id_col)] <- "sample_id"
        gene_pnps <- gene_pnps %>%
            left_join(metadata_for_join, by = "sample_id")

        n_with_meta <- sum(!is.na(gene_pnps\$disease_group))
        if ("disease_group" %in% names(gene_pnps)) {
            cat("Genes with metadata (disease_group):", n_with_meta, "of", nrow(gene_pnps), "\\n")
        }
        cat("Final gene_pnps columns:", paste(names(gene_pnps), collapse=", "), "\\n")

        saveRDS(gene_pnps, "pN_pS_gene_level.rds", compress = "xz")
    } else if (nrow(gene_combined) > 0) {
        # Fallback: check if gene_info already has pN/pS columns (some versions)
        pnps_cols <- c("SNV_N_count", "SNV_S_count", "N_sites", "S_sites")
        if (all(pnps_cols %in% colnames(gene_combined))) {
            gene_pnps <- gene_combined %>%
                select(gene, scaffold, sample_id, coverage, breadth, breadth_minCov,
                       SNV_N_count, SNV_S_count, N_sites, S_sites) %>%
                inner_join(stb_data, by = "scaffold") %>%
                mutate(
                    clean_genome = gsub("_genomic\\\\.fna(\\\\.gz)?\$", "", genome),
                    pN = ifelse(N_sites > 0, SNV_N_count / N_sites, NA),
                    pS = ifelse(S_sites > 0, SNV_S_count / S_sites, NA),
                    pNpS = ifelse(!is.na(pN) & !is.na(pS) & pS > 0, pN / pS, NA)
                )

            # Join sample metadata (fallback case)
            metadata_for_join <- metadata
            colnames(metadata_for_join)[which(colnames(metadata_for_join) == sample_id_col)] <- "sample_id"
            gene_pnps <- gene_pnps %>%
                left_join(metadata_for_join, by = "sample_id")

            saveRDS(gene_pnps, "pN_pS_gene_level.rds", compress = "xz")
        } else {
            cat("WARNING: genes_SNP_count.csv not found and gene_info lacks N_sites/S_sites\\n")
            cat("pN/pS calculation not possible - saving empty data\\n")
            saveRDS(data.frame(), "pN_pS_gene_level.rds")
            gene_pnps <- data.frame()
        }
    } else {
        saveRDS(data.frame(), "pN_pS_gene_level.rds")
        gene_pnps <- data.frame()
    }

    # Calculate genome-wide pN/pS (aggregate from gene-level)
    genome_wide_pnps <- data.frame()
    if (nrow(gene_pnps) > 0) {
        cat("Calculating genome-wide pN/pS...\\n")
        genome_wide_pnps <- gene_pnps %>%
            group_by(clean_genome, genome, sample_id) %>%
            summarise(
                total_SNV_N = sum(SNV_N_count, na.rm = TRUE),
                total_SNV_S = sum(SNV_S_count, na.rm = TRUE),
                total_N_sites = sum(N_sites, na.rm = TRUE),
                total_S_sites = sum(S_sites, na.rm = TRUE),
                n_genes = n(),
                mean_coverage = mean(coverage, na.rm = TRUE),
                mean_breadth = mean(breadth, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            mutate(
                pN = ifelse(total_N_sites > 0, total_SNV_N / total_N_sites, NA),
                pS = ifelse(total_S_sites > 0, total_SNV_S / total_S_sites, NA),
                pNpS = ifelse(!is.na(pN) & !is.na(pS) & pS > 0, pN / pS, NA)
            )
        cat("Genome-wide pN/pS rows:", nrow(genome_wide_pnps), "\\n")
        cat("Genome-wide pN/pS range:", round(min(genome_wide_pnps\$pNpS, na.rm = TRUE), 4), "to",
            round(max(genome_wide_pnps\$pNpS, na.rm = TRUE), 4), "\\n")
    } else {
        cat("WARNING: No gene-level pN/pS data - genome-wide pN/pS will be empty\\n")
    }
    saveRDS(genome_wide_pnps, "pN_pS_genome_wide.rds", compress = "xz")

    # eggNOG subset
    eggnog_file <- "${eggnog_annotations}"
    if (file.exists(eggnog_file) && file.size(eggnog_file) > 0) {
        con <- file(eggnog_file, "r")
        skip <- 0
        while(grepl("^##", readLines(con, 1))) skip <- skip + 1
        close(con)

        eggnog <- fread(eggnog_file, sep = "\\t", header = TRUE, skip = skip,
                        select = c("#query", "COG_category", "Description", "KEGG_ko", "KEGG_Pathway")) %>%
            as.data.frame() %>% rename(gene = `#query`)

        if (nrow(gene_pnps) > 0) {
            eggnog <- eggnog %>% filter(gene %in% unique(gene_pnps\$gene))
        }

        eggnog <- eggnog %>%
            mutate(COG_category_name = case_when(
                COG_category == "J" ~ "Translation",
                COG_category == "K" ~ "Transcription",
                COG_category == "L" ~ "Replication and repair",
                COG_category == "C" ~ "Energy production",
                COG_category == "E" ~ "Amino acid metabolism",
                COG_category == "G" ~ "Carbohydrate metabolism",
                COG_category == "M" ~ "Cell wall biogenesis",
                COG_category == "S" ~ "Function unknown",
                TRUE ~ "Other"
            ))
        saveRDS(eggnog, "eggnog_annotations_subset.rds", compress = "xz")
    } else {
        saveRDS(data.frame(), "eggnog_annotations_subset.rds")
    }

    # Integrated data
    genome_combined <- genome_combined %>%
        mutate(clean_genome = gsub("_genomic\\\\.fna(\\\\.gz)?\$", "", genome))

    join_by_cols <- setNames(sample_id_col, "sample_id")
    diversity_data <- genome_combined %>%
        left_join(metadata, by = join_by_cols)
    cat("Diversity data rows after metadata join:", nrow(diversity_data), "\\n")

    # Parse GTDB taxonomy - like preprocessing_scripts.R
    if ("gtdb_taxonomy" %in% names(gtdb_data)) {
        cat("Parsing GTDB taxonomy...\\n")

        # Helper function for root accession extraction
        extract_root <- function(x) {
            x <- as.character(x)
            m <- regexpr("GC[FA]_[0-9]+\\\\.[0-9]+", x, perl = TRUE)
            v <- regmatches(x, m)
            ifelse(length(v) == 0 | is.na(v) | v == "", NA_character_, v)
        }
        extract_root_vec <- Vectorize(extract_root)

        gtdb_parsed <- gtdb_data %>%
            mutate(
                clean_accession = gsub("^(RS_|GB_)", "", accession),
                domain = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=1) gsub("^d__", "", x[1]) else NA),
                phylum = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=2) gsub("^p__", "", x[2]) else NA),
                class = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=3) gsub("^c__", "", x[3]) else NA),
                order = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=4) gsub("^o__", "", x[4]) else NA),
                family = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=5) gsub("^f__", "", x[5]) else NA),
                genus = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=6) gsub("^g__", "", x[6]) else NA),
                species = sapply(strsplit(gtdb_taxonomy, ";"), function(x) if(length(x)>=7) gsub("^s__", "", x[7]) else NA),
                root_accession = extract_root_vec(clean_accession),
                accession_core = sub("^GC[FA]_", "", root_accession)
            )

        # Join GTDB data to diversity_data (include genome_size, gtdb_taxonomy, accession, etc.)
        gtdb_cols_for_diversity <- c("clean_accession", "accession", "genome_size", "gtdb_taxonomy",
                                     "domain", "phylum", "class", "order", "family", "genus", "species",
                                     "root_accession", "accession_core")
        gtdb_cols_available <- intersect(gtdb_cols_for_diversity, names(gtdb_parsed))

        diversity_data <- diversity_data %>%
            left_join(gtdb_parsed %>% select(all_of(gtdb_cols_available)),
                      by = c("clean_genome" = "clean_accession"))

        # Add calculated columns like preprocessing_scripts.R
        if ("length" %in% names(diversity_data) && "breadth" %in% names(diversity_data)) {
            diversity_data <- diversity_data %>%
                mutate(
                    analyzed_genome_length = length * breadth,
                    snvs_per_kbp = ifelse(
                        !is.na(SNV_count) & !is.na(analyzed_genome_length) & analyzed_genome_length > 0,
                        (SNV_count / analyzed_genome_length) * 1000, NA
                    )
                )
        }

        cat("GTDB taxonomy joined, diversity_data cols:", ncol(diversity_data), "\\n")
    }

    # ANI comparison data - Full processing like preprocessing_scripts.R
    compare_file <- list.files("instrainComparer_output/output", "genomeWide_compare.tsv", full.names = TRUE)
    ani_data <- data.frame()
    if (length(compare_file) > 0) {
        instrain_data <- fread(compare_file[1]) %>% as.data.frame()
        cat("Raw ANI comparison rows:", nrow(instrain_data), "\\n")

        if (nrow(instrain_data) > 0) {
            # Helper function to extract root accession
            extract_root <- function(x) {
                x <- as.character(x)
                m <- regexpr("GC[FA]_[0-9]+\\\\.[0-9]+", x, perl = TRUE)
                v <- regmatches(x, m)
                ifelse(length(v) == 0 | is.na(v) | v == "", NA_character_, v)
            }
            extract_root_vec <- Vectorize(extract_root)

            # Add clean columns like preprocessing_scripts.R
            instrain_data <- instrain_data %>%
                mutate(
                    clean_genome = gsub("_genomic\\\\.fna(\\\\.gz)?\$", "", genome),
                    clean_name1 = gsub("\\\\.sorted\\\\.bam\$", "", name1),
                    clean_name2 = gsub("\\\\.sorted\\\\.bam\$", "", name2),
                    root_acc = extract_root_vec(clean_genome),
                    core_acc = sub("^GC[FA]_", "", root_acc)
                )

            # Add root_accession and accession_core to gtdb_parsed for joining
            gtdb_parsed <- gtdb_parsed %>%
                mutate(
                    root_accession = extract_root_vec(clean_accession),
                    accession_core = sub("^GC[FA]_", "", root_accession)
                )

            # First pass: accession_core based matching
            ani_data <- instrain_data %>%
                left_join(
                    gtdb_parsed %>% select(accession_core, domain, phylum, class, order, family, genus, species),
                    by = c("core_acc" = "accession_core")
                )

            # Second pass: retry with root_accession if still NA
            na_idx <- is.na(ani_data\$species)
            if (any(na_idx)) {
                ani_fallback <- ani_data[na_idx, ] %>%
                    select(-domain, -phylum, -class, -order, -family, -genus, -species) %>%
                    left_join(
                        gtdb_parsed %>% select(root_accession, domain, phylum, class, order, family, genus, species),
                        by = c("root_acc" = "root_accession")
                    )
                ani_data[na_idx, c("domain", "phylum", "class", "order", "family", "genus", "species")] <-
                    ani_fallback[, c("domain", "phylum", "class", "order", "family", "genus", "species")]
            }

            # Join metadata for sample1 and sample2 (critical for Shiny app!)
            join_col <- sample_id_col
            ani_data <- ani_data %>%
                left_join(metadata, by = setNames(join_col, "clean_name1"), suffix = c("", "_sample1")) %>%
                left_join(metadata, by = setNames(join_col, "clean_name2"), suffix = c("_sample1", "_sample2"))

            cat("ANI data after full processing:", nrow(ani_data), "rows,", ncol(ani_data), "cols\\n")
        }
    } else {
        cat("No ANI comparison data found\\n")
    }

    # Add taxonomy and metadata to genome_wide_pnps (like preprocessing_scripts.R)
    if (nrow(genome_wide_pnps) > 0 && exists("gtdb_parsed")) {
        cat("Adding taxonomy and metadata to genome_wide_pnps...\\n")
        genome_wide_pnps <- genome_wide_pnps %>%
            left_join(
                gtdb_parsed %>% select(clean_accession, domain, phylum, class, order, family, genus, species),
                by = c("clean_genome" = "clean_accession")
            ) %>%
            left_join(metadata, by = setNames(sample_id_col, "sample_id"))
        cat("genome_wide_pnps after join:", nrow(genome_wide_pnps), "rows,", ncol(genome_wide_pnps), "cols\\n")

        # Re-save with taxonomy/metadata
        saveRDS(genome_wide_pnps, "pN_pS_genome_wide.rds", compress = "xz")
    }

    # Build pnps_data structure (compatible with shiny app)
    pnps_data <- list(
        genome_wide_pnps = genome_wide_pnps,
        gene_level_pnps = NULL,  # Saved separately to pN_pS_gene_level.rds
        gene_level_file = "pN_pS_gene_level.rds"
    )

    # Use gtdb_parsed (with taxonomy columns) instead of raw gtdb_data
    # Filter to only genomes present in diversity_data and add genome_size from diversity_data
    gtdb_final <- if (exists("gtdb_parsed")) gtdb_parsed else gtdb_data

    # Get unique genomes from diversity_data with their genome_size (length)
    detected_genomes <- diversity_data %>%
        select(clean_genome, length) %>%
        distinct(clean_genome, .keep_all = TRUE) %>%
        rename(genome_size = length)

    cat("Detected unique genomes:", nrow(detected_genomes), "\\n")

    # Filter gtdb_final to only detected genomes and add genome_size
    gtdb_final <- gtdb_final %>%
        inner_join(detected_genomes, by = c("clean_accession" = "clean_genome"))

    cat("Filtered gtdb_data rows:", nrow(gtdb_final), "\\n")

    # Rename sample ID column to 'accession_used_in_analysis' for Shiny app compatibility
    if (sample_id_col %in% names(metadata) && sample_id_col != "accession_used_in_analysis") {
        names(metadata)[names(metadata) == sample_id_col] <- "accession_used_in_analysis"
        cat("Renamed metadata column", sample_id_col, "to accession_used_in_analysis\\n")
    }

    integrated <- list(
        diversity_data = diversity_data,
        ani_data = ani_data,
        gtdb_data = gtdb_final,
        metadata = metadata,
        pnps_data = pnps_data,
        sample_id_column = "accession_used_in_analysis",  # Always use standard name
        processing_time = Sys.time()
    )
    saveRDS(integrated, "integrated_microbiome_data.rds", compress = "xz")

    cat("\\n=== Aggregation Complete! ===\\n")
    cat("Output files:\\n")
    cat("  - integrated_microbiome_data.rds\\n")
    cat("  - pN_pS_gene_level.rds\\n")
    cat("  - pN_pS_genome_wide.rds\\n")
    cat("  - eggnog_annotations_subset.rds\\n")
    """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {
    // GTDB metadata (always required for genome matching)
    gtdb_metadata_ch = Channel.fromPath(params.gtdb_metadata, checkIfExists: true)

    // External metadata is optional - will extract from phyloseq if not provided
    external_metadata_ch = params.metadata.isEmpty()
        ? Channel.value(file('NO_METADATA'))
        : Channel.fromPath(params.metadata, checkIfExists: true)

    // STAGE 1: Prevalence filtering
    if (!params.skip_prevalence) {
        // Check input type: phyloseq object or abundance table
        if (!params.phyloseq_object.isEmpty()) {
            phyloseq_ch = Channel.fromPath(params.phyloseq_object, checkIfExists: true)
            prevalence_filter_phyloseq(phyloseq_ch, external_metadata_ch)
            genome_paths_ch = prevalence_filter_phyloseq.out.genome_paths
            metadata_ch = prevalence_filter_phyloseq.out.sample_metadata  // Use extracted/provided metadata
        } else if (!params.abundance_table.isEmpty()) {
            abundance_ch = Channel.fromPath(params.abundance_table, checkIfExists: true)
            prevalence_filter_table(abundance_ch)
            genome_paths_ch = prevalence_filter_table.out.genome_paths
            // For table-based filtering, metadata is required
            if (params.metadata.isEmpty()) {
                error "ERROR: --metadata is required when using --abundance_table"
            }
            metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
        } else {
            error """ERROR: Either --phyloseq_object or --abundance_table is required.

Typical usage with WMS_TAXONOMY output:
  --phyloseq_object results/metagenome/WMS_TAXONOMY/phyloseq/phyloseq_object_sylph.RDS

Or with abundance table:
  --abundance_table results/metagenome/WMS_TAXONOMY/sylph/merged_sylph_species.tsv
"""
        }
    } else {
        // When skipping, check if files exist
        def prevalence_file = file("${params.outdir}/01_prevalent_taxa/prevalent_taxa_genome_paths.txt")
        def metadata_file = file("${params.outdir}/01_prevalent_taxa/sample_metadata.csv")
        if (!prevalence_file.exists()) {
            error "ERROR: --skip_prevalence requires existing file: ${prevalence_file}"
        }
        genome_paths_ch = Channel.fromPath(prevalence_file)
        metadata_ch = metadata_file.exists()
            ? Channel.fromPath(metadata_file)
            : (params.metadata.isEmpty() ? error("ERROR: --metadata required when skipping prevalence without existing sample_metadata.csv") : Channel.fromPath(params.metadata))
    }

    // STAGE 2: Genome preparation
    if (!params.skip_genome_prep && params.reference_fasta.isEmpty()) {
        concat_genomes(genome_paths_ch)
        generate_stb(genome_paths_ch)
        bowtie2_build(concat_genomes.out.combined_fasta)
        reference_ch = concat_genomes.out.combined_fasta
        stb_ch = generate_stb.out.stb_file
        bowtie2_index_ch = bowtie2_build.out.index.collect()

        // Create channel of individual genomes for parallel prodigal
        // Read genome paths file and create tuple(genome_id, genome_path) for each
        // Note: paths are container paths (/opt/database/...), existence check done inside container
        individual_genomes_ch = genome_paths_ch
            .splitText()
            .map { line ->
                def path = line.trim()
                if (path) {
                    def genome_name = path.tokenize('/').last().replaceAll(/(_genomic)?\.fn?a?(\.gz)?$/, '')
                    return tuple(genome_name, path)
                }
                return null
            }
            .filter { it != null }
    } else {
        // Using pre-existing files
        if (!params.reference_fasta.isEmpty()) {
            reference_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true)
        } else {
            def ref_file = file("${params.outdir}/02_genome_prep/all_genomes_combined.fa")
            if (!ref_file.exists()) {
                error "ERROR: --skip_genome_prep requires existing reference FASTA: ${ref_file}"
            }
            reference_ch = Channel.fromPath(ref_file)
        }

        stb_ch = params.stb_file.isEmpty()
            ? Channel.value(file('NO_STB'))
            : Channel.fromPath(params.stb_file, checkIfExists: true)

        if (!params.bowtie2_index.isEmpty()) {
            bowtie2_index_ch = Channel.fromPath("${params.bowtie2_index}*").collect()
        } else {
            def idx_pattern = "${params.outdir}/02_genome_prep/bowtie2_index/reference_index*"
            bowtie2_index_ch = Channel.fromPath(idx_pattern).collect()
        }

        // Create individual_genomes_ch from genome_paths_ch for parallel prodigal (if annotation needed)
        // Note: paths are container paths (/opt/database/...), existence check done inside container
        individual_genomes_ch = genome_paths_ch
            .splitText()
            .map { line ->
                def path = line.trim()
                if (path) {
                    def genome_name = path.tokenize('/').last().replaceAll(/(_genomic)?\.fn?a?(\.gz)?$/, '')
                    return tuple(genome_name, path)
                }
                return null
            }
            .filter { it != null }
    }

    // STAGE 3: Gene annotation (parallelized per-genome prodigal)
    if (!params.skip_annotation && params.gene_file.isEmpty()) {
        // Run prodigal on individual genomes in parallel (much faster than on combined fasta)
        prodigal_per_genome(individual_genomes_ch)

        // Collect all prodigal outputs and concatenate
        all_faa = prodigal_per_genome.out.map { genome_id, faa, fna, gff -> faa }.collect()
        all_fna = prodigal_per_genome.out.map { genome_id, faa, fna, gff -> fna }.collect()
        all_gff = prodigal_per_genome.out.map { genome_id, faa, fna, gff -> gff }.collect()

        concat_gene_predictions(all_faa, all_fna, all_gff)

        // Run eggNOG on combined proteins
        eggnog_mapper(concat_gene_predictions.out.proteins)
        gene_ch = concat_gene_predictions.out.genes
        eggnog_ch = eggnog_mapper.out.annotations
    } else if (!params.skip_annotation && !params.gene_file.isEmpty() && !params.protein_file.isEmpty() && params.eggnog_annotations.isEmpty()) {
        // Pre-existing gene file + protein file provided -> run eggNOG with protein file
        gene_ch = Channel.fromPath(params.gene_file, checkIfExists: true)
        protein_ch = Channel.fromPath(params.protein_file, checkIfExists: true)
        eggnog_mapper(protein_ch)
        eggnog_ch = eggnog_mapper.out.annotations
    } else {
        // Skip annotation or use pre-existing files
        gene_ch = params.gene_file.isEmpty()
            ? Channel.value(file('NO_GENES'))
            : Channel.fromPath(params.gene_file, checkIfExists: true)

        if (!params.eggnog_annotations.isEmpty()) {
            eggnog_ch = Channel.fromPath(params.eggnog_annotations, checkIfExists: true)
        } else {
            def eggnog_file = file("${params.outdir}/03_gene_annotation/eggnog_results.emapper.annotations")
            eggnog_ch = eggnog_file.exists()
                ? Channel.fromPath(eggnog_file)
                : Channel.value(file('NO_EGGNOG'))
        }
    }

    // STAGE 4: Read mapping
    if (params.bam_dir.isEmpty()) {
        if (params.input_dir.isEmpty()) {
            error "ERROR: Either --input_dir or --bam_dir is required"
        }
        reads_ch = Channel.fromFilePairs("${params.input_dir}/*_{1,2}.{fastq,fq}{.gz,}", checkIfExists: true)
            .map { acc, reads -> tuple(acc.replaceAll(/\_fastp.*$/, ''), reads) }
        bowtie2_mapping(reads_ch, bowtie2_index_ch, reference_ch.first())
        bam_ch = bowtie2_mapping.out
    } else {
        bam_ch = Channel.fromFilePairs("${params.bam_dir}/*.sorted.{bam,bam.bai}", checkIfExists: true)
            .map { basename, files ->
                def acc = basename.replaceAll(/\.sorted$/, '')
                def bam = files.find { it.name.endsWith('.bam') && !it.name.endsWith('.bam.bai') }
                def bai = files.find { it.name.endsWith('.bam.bai') }
                tuple(acc, bam, bai)
            }
    }

    // STAGE 5-6: InStrain
    if (!params.skip_instrain) {
        // Check if using pre-existing profile directories
        if (!params.profile_dir.isEmpty()) {
            // Use existing profile directories (skip instrain_profile step)
            log.info "Using existing InStrain profiles from: ${params.profile_dir}"
            all_profiles = Channel.fromPath("${params.profile_dir}/*_instrain_profile", type: 'dir')
                .collect()

            // Validate existing profiles
            validate_instrain_profiles(all_profiles)

            // Get valid profiles directly from validation (already copied)
            valid_profiles_ch = validate_instrain_profiles.out.valid_dirs
                .collect()
                .filter { it.size() >= 2 }
        } else {
            // Run instrain_profile from scratch
            instrain_profile(bam_ch, reference_ch.first(), stb_ch.first(), gene_ch.first())
            all_profiles = instrain_profile.out.map { acc, dir -> dir }.collect()

            // STAGE 5.5: Validate profiles before compare
            validate_instrain_profiles(all_profiles)

            // Get valid profiles directly from validation (already copied)
            valid_profiles_ch = validate_instrain_profiles.out.valid_dirs
                .collect()
                .filter { it.size() >= 2 }
        }

        // Run compare only with valid profiles (skip if < 2)
        valid_profiles_ch.ifEmpty {
            log.warn "Less than 2 valid profiles detected. Skipping inStrain compare."
        }

        instrain_compare(valid_profiles_ch, stb_ch.first())

        // STAGE 7: Aggregate for Shiny
        // Handle case when compare didn't run (use empty channel)
        compare_out_ch = instrain_compare.out.ifEmpty { file('NO_COMPARE') }
        aggregate_for_shiny(all_profiles, compare_out_ch, metadata_ch, stb_ch.first(), eggnog_ch.first(), gtdb_metadata_ch.first())
    }
}
