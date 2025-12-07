# gene_level_pnps_preprocessing_optimized.R
# OPTIMIZED: Faster loading using pre-filtered data

library(dplyr)

#' Load gene-level pN/pS data with COG annotations (OPTIMIZED)
#'
#' @param use_filtered Use pre-filtered data (recommended, much faster)
#' @param min_coverage Minimum coverage (only used if use_filtered = FALSE)
#' @param min_breadth Minimum breadth (only used if use_filtered = FALSE)
#' @return Data frame with gene-level pN/pS + COG annotations
load_gene_pnps_optimized <- function(
  use_filtered = TRUE,
  min_coverage = 5,
  min_breadth = 0.5
) {

  cat("=== Loading Gene-Level pN/pS Data (OPTIMIZED) ===\n")

  # === FASTEST PATH: Use expanded COG file (multi-letter categories split) ===
  expanded_file <- "data/pN_pS_gene_expanded_cog.rds"

  if (use_filtered && file.exists(expanded_file)) {
    cat("\n✓ Using expanded COG file with individual categories (FASTEST!)\n")
    cat("  Loading from:", expanded_file, "\n")

    data <- readRDS(expanded_file)

    cat("  ✓ Loaded", nrow(data), "rows (genes with multi-COG are expanded)\n")
    cat("  ✓ Individual COG categories:", length(unique(data$COG_category[!is.na(data$COG_category) & data$COG_category != "-"])), "\n")

    return(data)
  }

  # === FAST PATH: Use pre-integrated file (original with multi-letter COG) ===
  integrated_filtered_file <- "data/pN_pS_gene_filtered_with_cog.rds"

  if (use_filtered && file.exists(integrated_filtered_file)) {
    cat("\n✓ Using pre-filtered + integrated file (FAST)\n")
    cat("  Loading from:", integrated_filtered_file, "\n")
    cat("  NOTE: This file has multi-letter COG categories (e.g., 'ADL')\n")
    cat("  TIP: Run scripts/expand_cog_categories.R for individual category analysis\n")

    data <- readRDS(integrated_filtered_file)

    cat("  ✓ Loaded", nrow(data), "genes\n")
    cat("  ✓ With COG annotations:", sum(!is.na(data$COG_category) & data$COG_category != "-"), "\n")

    return(data)
  }

  # === MEDIUM PATH: Use pre-filtered file without COG ===
  filtered_file <- "data/pN_pS_gene_filtered.rds"

  if (use_filtered && file.exists(filtered_file)) {
    cat("\n✓ Using pre-filtered file (coverage >=5, breadth >=0.5)\n")
    cat("  Loading from:", filtered_file, "\n")

    gene_data <- readRDS(filtered_file)
    cat("  ✓ Loaded", nrow(gene_data), "genes\n")

  } else {
    # === SLOW PATH: Load full file and filter ===
    cat("\nℹ Pre-filtered file not found, using full gene file\n")
    cat("  TIP: Run scripts/create_filtered_gene_pnps.R for faster loading\n\n")

    full_file <- "data/pN_pS_gene_level.rds"
    if (!file.exists(full_file)) {
      stop("Gene-level data not found: ", full_file)
    }

    cat("  Loading full gene data (this is slow)...\n")
    gene_data_full <- readRDS(full_file)
    cat("  Loaded", nrow(gene_data_full), "genes\n")

    cat("  Filtering by quality (coverage >=", min_coverage,
        ", breadth >=", min_breadth, ")...\n")

    gene_data <- gene_data_full %>%
      filter(
        coverage >= min_coverage,
        breadth_minCov >= min_breadth,
        !is.na(pNpS),
        is.finite(pNpS)
      )

    cat("  Kept", nrow(gene_data), "genes after filtering\n")
    rm(gene_data_full)
    gc()
  }

  # === Load and join eggNOG annotations ===
  cat("\n  Loading eggNOG annotations...\n")

  eggnog_subset_file <- "data/eggnog_annotations_subset.rds"

  if (file.exists(eggnog_subset_file)) {
    cat("  Using eggNOG subset file...\n")
    eggnog_data <- readRDS(eggnog_subset_file)

  } else {
    cat("  WARNING: eggNOG subset not found!\n")
    cat("  TIP: Run scripts/preprocess_eggnog_once.R first\n")
    cat("  Returning data without COG annotations\n")
    return(gene_data)
  }

  cat("  Joining gene pN/pS with eggNOG...\n")
  integrated_data <- gene_data %>%
    left_join(eggnog_data, by = "gene")

  n_annotated <- sum(!is.na(integrated_data$COG_category) &
                     integrated_data$COG_category != "-")
  pct_annotated <- 100 * n_annotated / nrow(integrated_data)

  cat("\n=== Summary ===\n")
  cat("Total genes:", nrow(integrated_data), "\n")
  cat("With COG annotation:", n_annotated, paste0("(", round(pct_annotated, 1), "%)"), "\n")

  # Optionally save for future use
  if (!file.exists(integrated_filtered_file) && use_filtered) {
    cat("\nSaving integrated data for future use...\n")
    saveRDS(integrated_data, integrated_filtered_file, compress = "xz")
    cat("✓ Saved to:", integrated_filtered_file, "\n")
  }

  return(integrated_data)
}

# Convenience wrapper for backward compatibility
integrate_eggnog_with_gene_pnps <- function(...) {
  load_gene_pnps_optimized(use_filtered = TRUE)
}
