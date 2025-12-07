# pnps_preprocessing_optimized.R
# Optimized version using parallel processing and fast file reading

library(dplyr)
library(tidyr)
library(data.table)  # For fast file reading
library(parallel)    # For parallel processing

#' Load and process pN/pS data from InStrain gene files (OPTIMIZED)
#'
#' @param base_path Base directory containing InStrain output folders
#' @param gtdb_data GTDB taxonomy data for genome annotation
#' @param sample_metadata Sample metadata for integration
#' @param min_breadth_minCov Minimum breadth coverage threshold (default 0.5)
#' @param n_cores Number of CPU cores to use (default: detectCores() - 1)
#' @param stb_file Path to scaffold-to-bin (.stb) file for scaffold-to-genome mapping
#' @return List containing genome_wide_pnps and gene_level_pnps data frames
load_pnps_data <- function(base_path, gtdb_data, sample_metadata,
                           min_breadth_minCov = 0.5,
                           n_cores = NULL,
                           stb_file = Sys.getenv("STB_FILE_PATH", "/data/instrain/prevalent_taxa.stb")) {

  cat("=== Loading pN/pS Data (OPTIMIZED) ===\n")
  cat("Base path:", base_path, "\n")
  cat("Minimum breadth_minCov:", min_breadth_minCov, "\n")

  # Load scaffold-to-genome mapping
  cat("Loading scaffold-to-genome mapping from .stb file...\n")
  if(!file.exists(stb_file)) {
    stop("Scaffold-to-genome mapping file not found: ", stb_file)
  }

  scaffold_to_genome <- tryCatch({
    fread(stb_file, header = FALSE, showProgress = FALSE) %>%
      as.data.frame() %>%
      setNames(c("scaffold", "genome"))
  }, error = function(e) {
    stop("Error reading .stb file: ", e$message)
  })

  cat("  Loaded", nrow(scaffold_to_genome), "scaffold-to-genome mappings\n")

  # Determine number of cores
  if(is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)
  }
  cat("Using", n_cores, "CPU cores for parallel processing\n")

  # Find all InStrain directories
  instrain_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
  instrain_dirs <- instrain_dirs[grepl("_0.92$", instrain_dirs)]

  cat("Found", length(instrain_dirs), "InStrain directories\n")

  if(length(instrain_dirs) == 0) {
    stop("No InStrain directories found at: ", base_path)
  }

  # Define function to process a single sample
  process_sample <- function(instrain_dir, scaffold_to_genome, show_progress = FALSE) {

    sample_name <- basename(instrain_dir)
    sample_name <- gsub("_0.92$", "", sample_name)
    sample_name <- gsub("_fastp_hg38.sorted_instrain_profile", "", sample_name)

    if(show_progress) {
      cat("  Processing:", sample_name, "...")
    }

    # Paths to required files
    gene_info_file <- file.path(instrain_dir, "output",
                                paste0(basename(instrain_dir), "_gene_info.tsv"))

    # Check for genes_SNP_count.csv (uncompressed or compressed)
    gene_snp_file <- file.path(instrain_dir, "raw_data", "genes_SNP_count.csv")
    gene_snp_file_gz <- file.path(instrain_dir, "raw_data", "genes_SNP_count.csv.gz")

    if(file.exists(gene_snp_file_gz)) {
      gene_snp_file <- gene_snp_file_gz
    }

    # Check if files exist
    if(!file.exists(gene_info_file) || !file.exists(gene_snp_file)) {
      if(show_progress) cat(" [SKIP: files not found]\n")
      return(list(genome_wide = NULL, gene_level = NULL, sample_name = sample_name))
    }

    # Read gene_info.tsv using fread (MUCH FASTER)
    gene_info <- tryCatch({
      fread(gene_info_file, showProgress = FALSE) %>%
        as.data.frame() %>%
        select(gene, scaffold, coverage, breadth, breadth_minCov) %>%
        mutate(
          breadth_minCov = as.numeric(breadth_minCov),
          coverage = as.numeric(coverage),
          breadth = as.numeric(breadth)
        )
    }, error = function(e) {
      if(show_progress) cat(" [ERROR reading gene_info:", e$message, "]\n")
      return(NULL)
    })

    if(is.null(gene_info) || nrow(gene_info) == 0) {
      if(show_progress) cat(" [SKIP: no gene_info data]\n")
      return(list(genome_wide = NULL, gene_level = NULL, sample_name = sample_name))
    }

    # Read genes_SNP_count.csv using fread (MUCH FASTER)
    gene_snp <- tryCatch({
      fread(gene_snp_file, showProgress = FALSE) %>%
        as.data.frame() %>%
        select(gene, SNV_N_count, SNV_S_count, N_sites, S_sites) %>%
        mutate(
          SNV_N_count = as.numeric(SNV_N_count),
          SNV_S_count = as.numeric(SNV_S_count),
          N_sites = as.numeric(N_sites),
          S_sites = as.numeric(S_sites)
        )
    }, error = function(e) {
      if(show_progress) cat(" [ERROR reading SNP data:", e$message, "]\n")
      return(NULL)
    })

    if(is.null(gene_snp) || nrow(gene_snp) == 0) {
      if(show_progress) cat(" [SKIP: no SNP count data]\n")
      return(list(genome_wide = NULL, gene_level = NULL, sample_name = sample_name))
    }

    # Join gene_info with genes_SNP_count and scaffold-to-genome mapping
    gene_data <- gene_info %>%
      inner_join(gene_snp, by = "gene") %>%
      inner_join(scaffold_to_genome, by = "scaffold") %>%
      filter(breadth_minCov >= min_breadth_minCov)

    if(nrow(gene_data) == 0) {
      if(show_progress) cat(" [SKIP: no genes pass breadth filter]\n")
      return(list(genome_wide = NULL, gene_level = NULL, sample_name = sample_name))
    }

    # Create clean_genome for GTDB taxonomy matching
    gene_data <- gene_data %>%
      mutate(clean_genome = gsub("_genomic\\.fna(\\.gz)?$", "", genome))

    # Calculate genome-wide pN/pS (grouped by genome)
    genome_wide <- gene_data %>%
      group_by(clean_genome, genome) %>%
      summarise(
        total_SNV_N = sum(SNV_N_count, na.rm = TRUE),
        total_SNV_S = sum(SNV_S_count, na.rm = TRUE),
        total_N_sites = sum(N_sites, na.rm = TRUE),
        total_S_sites = sum(S_sites, na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
      ) %>%
      mutate(
        pN = ifelse(total_N_sites > 0, total_SNV_N / total_N_sites, NA),
        pS = ifelse(total_S_sites > 0, total_SNV_S / total_S_sites, NA),
        pNpS = ifelse(!is.na(pN) & !is.na(pS) & pS > 0, pN / pS, NA),
        sample_id = sample_name
      )

    # Calculate gene-level pN/pS (kept for downstream analysis)
    gene_level <- gene_data %>%
      mutate(
        pN = ifelse(N_sites > 0, SNV_N_count / N_sites, NA),
        pS = ifelse(S_sites > 0, SNV_S_count / S_sites, NA),
        pNpS = ifelse(!is.na(pN) & !is.na(pS) & pS > 0, pN / pS, NA),
        sample_id = sample_name
      ) %>%
      select(sample_id, clean_genome, genome, scaffold, gene, coverage, breadth, breadth_minCov,
             SNV_N_count, SNV_S_count, N_sites, S_sites, pN, pS, pNpS)

    if(show_progress) {
      cat(" [OK:", nrow(genome_wide), "genomes,",
          nrow(gene_level), "genes, avg pN/pS=",
          round(mean(genome_wide$pNpS, na.rm = TRUE), 4), "]\n", sep=" ")
    }

    return(list(
      genome_wide = genome_wide,
      gene_level = gene_level,
      sample_name = sample_name
    ))
  }

  # Process samples (parallel or sequential)
  if(n_cores == 1) {
    cat("\nProcessing samples sequentially (safer for Shiny)...\n")
    results <- lapply(instrain_dirs, function(dir) {
      process_sample(dir, scaffold_to_genome, show_progress = TRUE)
    })
  } else {
    cat("\nProcessing samples in parallel using", n_cores, "cores...\n")
    results <- mclapply(instrain_dirs, function(dir) {
      process_sample(dir, scaffold_to_genome, show_progress = FALSE)
    }, mc.cores = n_cores)
  }

  # Extract results
  genome_wide_list <- list()
  gene_level_list <- list()

  successful <- 0
  for(result in results) {
    if(!is.null(result$genome_wide)) {
      genome_wide_list[[result$sample_name]] <- result$genome_wide
      gene_level_list[[result$sample_name]] <- result$gene_level
      successful <- successful + 1
    }
  }

  cat("\nSuccessfully processed:", successful, "samples\n")

  if(successful == 0) {
    stop("No genome-wide pN/pS data loaded!")
  }

  # Combine all samples
  genome_wide_pnps <- bind_rows(genome_wide_list)
  gene_level_pnps <- bind_rows(gene_level_list)

  cat("Total samples with genome-wide pN/pS:", nrow(genome_wide_pnps), "\n")
  cat("Total gene-level records:", nrow(gene_level_pnps), "\n")

  # ===== Save gene-level data BEFORE metadata/taxonomy integration =====
  # (For on-demand processing - keeps file small and loading fast)
  cat("\n=== Saving Gene-Level Data (Raw) ===\n")
  gene_level_file <- "data/pN_pS_gene_level.rds"
  saveRDS(gene_level_pnps, gene_level_file, compress = "xz")  # High compression
  cat("  ✓ Saved", nrow(gene_level_pnps), "gene-level records to:", gene_level_file, "\n")
  cat("  ✓ File size:", format(file.size(gene_level_file) / 1024^2, digits = 1), "MB\n")

  # Clear gene-level from memory (not needed for genome-wide analysis)
  rm(gene_level_pnps)
  gc()

  # ===== Integrate Genome-Wide with Metadata =====
  cat("\n=== Integrating Genome-Wide Data ===\n")

  if(!is.null(sample_metadata)) {
    # Join genome-wide data with metadata
    genome_wide_pnps <- genome_wide_pnps %>%
      left_join(sample_metadata, by = c("sample_id" = "accession_used_in_analysis"))

    cat("Genome-wide data columns after metadata join:", ncol(genome_wide_pnps), "\n")
  }

  # ===== Join with GTDB taxonomy =====
  if(!is.null(gtdb_data)) {
    cat("\n=== Integrating GTDB Taxonomy ===\n")

    # Join genome-wide pN/pS with GTDB taxonomy
    genome_wide_pnps <- genome_wide_pnps %>%
      left_join(
        gtdb_data %>% select(clean_accession, domain, phylum, class, order, family, genus, species),
        by = c("clean_genome" = "clean_accession")
      )

    cat("Genome-wide data with taxonomy:", sum(!is.na(genome_wide_pnps$genus)), "records\n")
  }

  # Summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat("Total genome-sample combinations:", nrow(genome_wide_pnps), "\n")
  cat("Unique genomes:", length(unique(genome_wide_pnps$clean_genome)), "\n")
  cat("Unique samples:", length(unique(genome_wide_pnps$sample_id)), "\n")
  cat("Genome-wide pN/pS range:",
      round(min(genome_wide_pnps$pNpS, na.rm = TRUE), 4), "to",
      round(max(genome_wide_pnps$pNpS, na.rm = TRUE), 4), "\n")
  cat("Genome-wide pN/pS median:",
      round(median(genome_wide_pnps$pNpS, na.rm = TRUE), 4), "\n")

  return(list(
    genome_wide_pnps = genome_wide_pnps,
    gene_level_pnps = NULL,  # Saved separately to data/pN_pS_gene_level.rds
    gene_level_file = gene_level_file  # Path to gene-level data for on-demand loading
  ))
}
