#!/usr/bin/env Rscript

# =============================================================================
# Preprocessing Script - Exact Matching Version
# =============================================================================

run_full_preprocessing <- function(
  gtdb_file = Sys.getenv("GTDB_METADATA_PATH", "/data/gtdb/bac_ar_combined_metadata.tsv"),
  instrain_file = NULL,  # Auto-detect from instrain_base_path
  metadata_file = Sys.getenv("SAMPLE_METADATA_PATH", "/data/metadata/sample_metadata.txt"),
  instrain_base_path = Sys.getenv("INSTRAIN_DATA_PATH", "/data/instrain/")
) {
  cat("========================================\n")
  cat("Starting Preprocessing\n")
  cat("========================================\n")

  start_time <- Sys.time()

  dir.create("data", showWarnings = FALSE)
  dir.create("temp", showWarnings = FALSE)

  base_path <- instrain_base_path
  
  # ==========================================================================
  # 1. Consolidate Genome Info
  # ==========================================================================
  cat("\n1. Consolidating genome info...\n")
  
  genome_dirs <- list.dirs(base_path, recursive = FALSE)
  genome_dirs <- genome_dirs[grepl("_0\\.92$", genome_dirs)]
  
  genome_files <- file.path(genome_dirs, "output", 
                            paste0(basename(genome_dirs), "_genome_info.tsv"))
  genome_files <- genome_files[file.exists(genome_files)]
  
  cat("   Number of files:", length(genome_files), "\n")
  
  combined_genome_info <- do.call(rbind, lapply(genome_files, function(f) {
    d <- read.delim(f, stringsAsFactors = FALSE)
    d$sample_id <- gsub("_fastp_hg38\\.sorted_instrain_profile.*", "", basename(dirname(dirname(f))))
    d$clean_genome <- gsub("_genomic\\.fna(\\.gz)?$", "", d$genome)
    
    return(d[!is.na(d$coverage) & d$coverage >= 1 & 
             !is.na(d$breadth) & d$breadth >= 0.1 & 
             !is.na(d$nucl_diversity) & d$nucl_diversity > 0, ])
  }))
  
  required_genomes <- unique(combined_genome_info$clean_genome)
  
  saveRDS(combined_genome_info, "data/combined_genome_info.rds")
  cat("   ✓ Genome info:", nrow(combined_genome_info), "records\n")
  cat("   ✓ Unique genomes:", length(required_genomes), "\n")

  # ==========================================================================
  # 2. GTDB Filtering (Exact matching in R)
  # ==========================================================================
  cat("\n2. Filtering GTDB data...\n")

  # Auto-detect GTDB file format (metafun 3-column vs full 113-column)
  # First read header to detect format
  header_line <- readLines(gtdb_file, n = 1)
  n_cols <- length(strsplit(header_line, "\t")[[1]])
  cat("   Detected GTDB format:", n_cols, "columns\n")

  if(n_cols <= 5) {
    # metafun simplified format: accession, gtdb_taxonomy, taxonomy_id
    gtdb_data <- read.delim(gtdb_file, stringsAsFactors = FALSE)
    cat("   Using metafun simplified GTDB format\n")
  } else {
    # Full GTDB format (113 columns) - read only necessary columns
    col_classes <- c("character", rep("NULL", 15), "numeric",
                     rep("NULL", 2), "character", rep("NULL", 93))
    gtdb_data <- read.delim(gtdb_file, stringsAsFactors = FALSE,
                           colClasses = col_classes)
    cat("   Using full GTDB format\n")
  }

  cat("   Total GTDB records:", nrow(gtdb_data), "\n")

  # Exact matching: remove RS_/GB_ prefix then exact match
  gtdb_data$clean_accession <- gsub("^(RS_|GB_)", "", gtdb_data$accession)
  gtdb_data <- gtdb_data[gtdb_data$clean_accession %in% required_genomes, ]

  cat("   ✓ After filtering:", nrow(gtdb_data), "records\n")
  cat("   ✓ Match rate:", round(nrow(gtdb_data)/length(required_genomes)*100, 1), "%\n")

  # ==========================================================================
  # 3. Parse GTDB Taxonomy
  # ==========================================================================
  cat("\n3. Parsing GTDB taxonomy...\n")
  
  taxonomy_split <- strsplit(gtdb_data$gtdb_taxonomy, ";", fixed = TRUE)
  tax_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  for(i in 1:7) {
    gtdb_data[, tax_levels[i]] <- sapply(taxonomy_split, function(x) {
      if(length(x) >= i && !is.na(x[i]) && x[i] != "" && 
         !grepl("unclassified", x[i], ignore.case = TRUE)) {
        return(x[i])
      }
      return(NA)
    })
  }
  # Add accession transformation for joining
  extract_root <- function(x) {
    x <- as.character(x)
    m <- regexpr("GC[FA]_\\d+\\.\\d+", x, perl = TRUE)
    v <- regmatches(x, m)
    v[is.na(v) | v == ""] <- NA_character_
    return(v)
  }
  gtdb_data$root_accession <- extract_root(gtdb_data$clean_accession)
  gtdb_data$accession_core <- sub("^GC[FA]_", "", gtdb_data$root_accession)

  saveRDS(gtdb_data, "data/filtered_gtdb_metadata.rds")
  cat("   ✓ Final GTDB:", nrow(gtdb_data), "records (",
      length(unique(gtdb_data$clean_accession)), "unique genomes)\n")

  # ==========================================================================
  # 4. InStrain Comparison & Metadata
  # ==========================================================================
  cat("\n4. Processing InStrain comparison & metadata...\n")

  # Auto-detect InStrain comparison file if not specified
  if(is.null(instrain_file) || !file.exists(instrain_file)) {
    cat("   Searching for InStrain comparison files...\n")
    compare_files <- list.files(instrain_base_path,
                                pattern = "genomeWide_compare\\.tsv$",
                                recursive = TRUE, full.names = TRUE)
    if(length(compare_files) > 0) {
      # Use the largest file (most comprehensive comparison)
      file_sizes <- file.info(compare_files)$size
      instrain_file <- compare_files[which.max(file_sizes)]
      cat("   Found:", basename(instrain_file), "\n")
    }
  }

  instrain_data <- NULL
  if(!is.null(instrain_file) && file.exists(instrain_file)) {
    instrain_data <- read.delim(instrain_file, stringsAsFactors = FALSE)
    instrain_data$clean_genome <- gsub("_genomic\\.fna(\\.gz)?$", "", instrain_data$genome)
    instrain_data$clean_name1 <- gsub("_fastp_hg38\\.sorted\\.bam$", "", instrain_data$name1)
    instrain_data$clean_name2 <- gsub("_fastp_hg38\\.sorted\\.bam$", "", instrain_data$name2)
    saveRDS(instrain_data, "data/processed_instrain_comparison.rds")
    cat("   ✓ InStrain comparison:", nrow(instrain_data), "records\n")
  } else {
    cat("   ⚠ No InStrain comparison file found\n")
  }

  # Read metadata - support both CSV and TSV
  if(grepl("\\.csv$", metadata_file, ignore.case = TRUE)) {
    metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  } else {
    metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
  }
  saveRDS(metadata, "data/processed_metadata.rds")
  cat("   ✓ Metadata:", nrow(metadata), "samples\n")

  # Detect the sample ID column in metadata
  sample_id_col <- NULL
  possible_cols <- c("sample_id", "accession_used_in_analysis", "Sample_ID", "SampleID", "accession")
  for(col in possible_cols) {
    if(col %in% names(metadata)) {
      sample_id_col <- col
      break
    }
  }
  if(is.null(sample_id_col)) {
    sample_id_col <- names(metadata)[1]
    cat("   ⚠ Using first column as sample ID:", sample_id_col, "\n")
  } else {
    cat("   ✓ Detected sample ID column:", sample_id_col, "\n")
  }

  # ==========================================================================
  # 5. Final Integration
  # ==========================================================================
  cat("\n5. Final integration...\n")

  diversity_data <- merge(combined_genome_info, gtdb_data,
                         by.x = "clean_genome", by.y = "clean_accession", all.x = TRUE)
  diversity_data <- merge(diversity_data, metadata,
                         by.x = "sample_id", by.y = sample_id_col, all.x = TRUE)
  
  if("length" %in% colnames(diversity_data) && "breadth" %in% colnames(diversity_data)) {
    diversity_data$analyzed_genome_length <- diversity_data$length * diversity_data$breadth
    diversity_data$snvs_per_kbp <- ifelse(
      !is.na(diversity_data$SNV_count) & !is.na(diversity_data$analyzed_genome_length) & 
      diversity_data$analyzed_genome_length > 0,
      (diversity_data$SNV_count / diversity_data$analyzed_genome_length) * 1000, NA
    )
  }
  
  # Data integration for ANI analysis (enhanced version)
  ani_data <- NULL
  if(!is.null(instrain_data)) {
    # Add accession_core to instrain_data
    instrain_data$root_acc <- extract_root(instrain_data$clean_genome)
    instrain_data$core_acc <- sub("^GC[FA]_", "", instrain_data$root_acc)

    # First pass: accession_core based matching
    ani_data <- merge(
      instrain_data,
      gtdb_data[, c("accession_core","domain","phylum","class","order","family","genus","species")],
      by.x = "core_acc", by.y = "accession_core", all.x = TRUE
    )

    # Second pass: retry with root_accession if still NA
    na_idx <- is.na(ani_data$species)
    if(any(na_idx)) {
      ani_fallback_in <- ani_data[na_idx, !(names(ani_data) %in% c("domain","phylum","class","order","family","genus","species"))]
      ani_fallback_out <- merge(
        ani_fallback_in,
        gtdb_data[, c("root_accession","domain","phylum","class","order","family","genus","species")],
        by.x = "root_acc", by.y = "root_accession", all.x = TRUE
      )
      ani_data[na_idx, c("domain","phylum","class","order","family","genus","species")] <-
        ani_fallback_out[, c("domain","phylum","class","order","family","genus","species")]
    }

    # Join metadata (per sample) using detected sample_id column
    ani_data <- merge(ani_data, metadata,
                      by.x = "clean_name1", by.y = sample_id_col,
                      all.x = TRUE, suffixes = c("", "_sample1"))
    ani_data <- merge(ani_data, metadata,
                      by.x = "clean_name2", by.y = sample_id_col,
                      all.x = TRUE, suffixes = c("_sample1", "_sample2"))
  }
  
  # ==========================================================================
  # 6. pN/pS Data Loading (OPTIMIZED)
  # ==========================================================================
  cat("\n6. Loading pN/pS data (optimized)...\n")

  # Check if separate pN/pS file already exists
  pnps_rds_file <- "data/pN_pS_integrated_data.rds"

  if(file.exists(pnps_rds_file)) {
    cat("  Loading existing pN/pS data from:", pnps_rds_file, "\n")
    pnps_data <- readRDS(pnps_rds_file)
    cat("  ✓ Loaded pN/pS data from cache\n")
  } else {
    cat("  No cached pN/pS data found, running preprocessing...\n")

    # Source optimized pN/pS preprocessing function
    source("modules/pnps_preprocessing_optimized.R")

    pnps_data <- tryCatch({
      load_pnps_data(
        base_path = base_path,
        gtdb_data = gtdb_data,
        sample_metadata = metadata,
        min_breadth_minCov = 0.5,
        n_cores = 1  # Sequential mode (parallel doesn't work in Shiny)
      )
    }, error = function(e) {
      cat("  WARNING: pN/pS data loading failed:", e$message, "\n")
      list(genome_wide_pnps = NULL, gene_level_pnps = NULL)
    })

    # Save pN/pS data separately
    saveRDS(pnps_data, pnps_rds_file)
    cat("  ✓ Saved pN/pS data to:", pnps_rds_file, "\n")
  }

  final_data <- list(
    diversity_data = diversity_data,
    ani_data = ani_data,
    gtdb_data = gtdb_data,
    metadata = metadata,
    pnps_data = pnps_data,
    sample_id_column = sample_id_col,
    processing_time = Sys.time()
  )

  saveRDS(final_data, "data/integrated_microbiome_data.rds")
  
  end_time <- Sys.time()
  total_time <- round(as.numeric(end_time - start_time, units = "mins"), 1)

  cat("\n========================================\n")
  cat("Preprocessing completed -", total_time, "minutes\n")
  cat("   ✓ Diversity:", nrow(diversity_data), "records\n")
  if(!is.null(ani_data)) cat("   ✓ ANI:", nrow(ani_data), "records\n")
  cat("   ✓ Saved: data/integrated_microbiome_data.rds\n")
  cat("========================================\n")
  
  return(final_data)
}

# Only run preprocessing in non-interactive mode if cached data doesn't exist
# This prevents automatic preprocessing when Shiny app is launched via Rscript
if(!interactive() && !file.exists("data/integrated_microbiome_data.rds")) {
  message("No cached data found. Running preprocessing...")
  run_full_preprocessing()
} else if(!interactive()) {
  message("Cached data found at data/integrated_microbiome_data.rds. Skipping preprocessing.")
}