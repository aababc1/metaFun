# ==============================================================================
# FastSpar Wrapper Functions - SIF Container Version
# FastSpar is installed in /opt/conda/bin (in PATH)
# ==============================================================================

#' Run FastSpar correlation analysis
#'
#' @param ps Phyloseq object
#' @param iterations Number of iterations for FastSpar (default: 20)
#' @param bootstraps Number of bootstrap replicates for p-values (default: 100)
#' @param threads Number of threads to use (default: 4)
#' @param fastspar_path Path to fastspar binary (default: auto-detect)
#' @return List containing correlation matrix and p-value matrix
run_fastspar <- function(ps, iterations = 20, bootstraps = 100, threads = 4,
                         fastspar_path = NULL) {

  # Auto-detect fastspar if not provided
  if (is.null(fastspar_path)) {
    # In container, fastspar is in /opt/conda/bin which is in PATH
    fastspar_path <- Sys.which("fastspar")
    if (fastspar_path == "") {
      # Fallback to explicit container path
      fastspar_path <- "/opt/conda/bin/fastspar"
      if (!file.exists(fastspar_path)) {
        stop("FastSpar not found. Please check container installation.")
      }
    }
  }

  # Create temporary directory
  tmp_dir <- tempdir()
  otu_file <- file.path(tmp_dir, "fastspar_input.tsv")
  cor_file <- file.path(tmp_dir, "fastspar_cor.tsv")
  cov_file <- file.path(tmp_dir, "fastspar_cov.tsv")
  boot_dir <- file.path(tmp_dir, "fastspar_bootstrap")
  pval_file <- file.path(tmp_dir, "fastspar_pvals.tsv")

  # Prepare OTU table (taxa as rows, samples as columns)
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # =========================================================================
  # CRITICAL: FastSpar requires INTEGER count data
  # If data is relative abundance (decimals), convert to pseudo-counts
  # =========================================================================

  # Check if data contains non-integer values
  is_integer_data <- all(otu_mat == floor(otu_mat), na.rm = TRUE)

  if (!is_integer_data) {
    cat("  [FastSpar] Detected non-integer data (relative abundance). Converting to pseudo-counts...\n")

    # Find the smallest non-zero value to determine scaling factor
    non_zero_values <- otu_mat[otu_mat > 0]
    if (length(non_zero_values) == 0) {
      stop("OTU table contains only zeros!")
    }

    min_nonzero <- min(non_zero_values)
    max_value <- max(otu_mat)

    cat(sprintf("    Min non-zero value: %.10f\n", min_nonzero))
    cat(sprintf("    Max value: %.4f\n", max_value))

    # Calculate scaling factor: 10^n where n is chosen so that
    # min_nonzero * 10^n >= 1 (to preserve smallest values as at least 1)
    # But also ensure max_value * 10^n doesn't overflow
    n_decimals <- ceiling(-log10(min_nonzero))
    scaling_factor <- 10^n_decimals

    # Cap scaling factor to avoid overflow (max reasonable integer ~2 billion)
    max_safe_scaling <- floor(2e9 / max_value)
    if (scaling_factor > max_safe_scaling) {
      scaling_factor <- max_safe_scaling
      cat(sprintf("    Capped scaling factor to prevent overflow: %.0f\n", scaling_factor))
    }

    cat(sprintf("    Scaling factor: 10^%d = %.0f\n", n_decimals, scaling_factor))

    # Scale and round to integers
    otu_mat <- round(otu_mat * scaling_factor)

    # Ensure all values are non-negative integers
    otu_mat[otu_mat < 0] <- 0
    otu_mat <- as.matrix(apply(otu_mat, 2, as.integer))
    rownames(otu_mat) <- rownames(as(otu_table(ps), "matrix"))
    if (!taxa_are_rows(ps)) {
      rownames(otu_mat) <- colnames(as(otu_table(ps), "matrix"))
    }

    cat(sprintf("    Converted range: %d - %d\n", min(otu_mat), max(otu_mat)))
    cat(sprintf("    Non-zero entries: %d / %d\n", sum(otu_mat > 0), length(otu_mat)))
  } else {
    cat("  [FastSpar] Data is already integer counts. Using directly.\n")
  }

  # Write OTU table in FastSpar format (tab-delimited, taxa names in first column)
  otu_df <- data.frame(
    OTU_ID = rownames(otu_mat),
    otu_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(otu_df, otu_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)

  cat("Running FastSpar correlation estimation...\n")

  # Run FastSpar correlation (direct call - no conda run in container)
  cmd_cor <- sprintf(
    "%s --otu_table %s --correlation %s --covariance %s --iterations %d --threads %d",
    fastspar_path, otu_file, cor_file, cov_file, iterations, threads
  )

  result_code <- system(cmd_cor)

  if (result_code != 0 || !file.exists(cor_file)) {
    stop("FastSpar correlation failed. Check input data.")
  }

  cat("Running FastSpar bootstrap for p-values...\n")

  # Create bootstrap directories
  boot_counts_dir <- file.path(boot_dir, "counts")
  boot_cor_dir <- file.path(boot_dir, "cors")
  dir.create(boot_counts_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(boot_cor_dir, showWarnings = FALSE, recursive = TRUE)

  # Get fastspar_bootstrap path (same directory as fastspar)
  fastspar_bootstrap <- paste0(fastspar_path, "_bootstrap")
  fastspar_pvalues <- paste0(fastspar_path, "_pvalues")

  # Step 1: Generate bootstrap count tables (with threading)
  cat("  Generating", bootstraps, "bootstrap samples...\n")
  cmd_boot <- sprintf(
    "%s --otu_table %s --number %d --prefix %s/boot --threads %d",
    fastspar_bootstrap, otu_file, bootstraps, boot_counts_dir, threads
  )
  system(cmd_boot)

  # Step 2: Run fastspar on each bootstrap sample (PARALLELIZED via shell)
  cat("  Running correlations on bootstrap samples (parallel with", threads, "cores)...\n")
  boot_files <- list.files(boot_counts_dir, pattern = "^boot", full.names = TRUE)

  if (length(boot_files) == 0) {
    stop("No bootstrap files generated")
  }

  # Use shell-level parallelization with xargs (more efficient for system calls)
  # Create a script file with all the commands to run
  parallel_script <- file.path(tmp_dir, "run_boot_correlations.sh")

  boot_commands <- sapply(boot_files, function(boot_file) {
    boot_id <- basename(boot_file)
    boot_cor_file <- file.path(boot_cor_dir, paste0("cor_", boot_id))
    boot_cov_file <- file.path(boot_cor_dir, paste0("cov_", boot_id))
    sprintf("%s --otu_table %s --correlation %s --covariance %s --iterations %d --threads 1",
            fastspar_path, boot_file, boot_cor_file, boot_cov_file, iterations)
  })

  # Write commands to file and run with xargs for parallel execution
  writeLines(boot_commands, parallel_script)

  cmd_parallel <- sprintf("cat %s | xargs -P %d -I {} sh -c '{} 2>/dev/null'",
                          parallel_script, threads)
  system(cmd_parallel, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Clean up script
  unlink(parallel_script)

  # Step 3: Calculate p-values from bootstrap correlations
  # NOTE: --prefix must match the actual filenames created in Step 2
  # Files are named: cor_boot_0.tsv, cor_boot_1.tsv, etc. (with underscore)
  # So prefix should be: boot_cor_dir/cor_boot_ (with trailing underscore)
  cat("  Calculating p-values from", length(boot_files), "bootstrap correlations...\n")
  cmd_pval <- sprintf(
    "%s --otu_table %s --correlation %s --prefix %s/cor_boot_ --permutations %d --outfile %s --threads %d",
    fastspar_pvalues, otu_file, cor_file, boot_cor_dir, bootstraps, pval_file, threads
  )

  system(cmd_pval)

  if (!file.exists(pval_file)) {
    stop("FastSpar p-value calculation failed.")
  }

  # Read results
  cor_matrix <- as.matrix(read.table(cor_file, sep = "\t", header = TRUE,
                                     row.names = 1, check.names = FALSE))
  p_matrix <- as.matrix(read.table(pval_file, sep = "\t", header = TRUE,
                                   row.names = 1, check.names = FALSE))

  # Handle taxa that FastSpar may have filtered out (BEFORE validation)
  n_taxa_expected <- nrow(otu_mat)
  n_taxa_result <- nrow(cor_matrix)

  if (n_taxa_result != n_taxa_expected) {
    warning(sprintf("FastSpar filtered some taxa: input=%d, output=%d. Padding matrices.",
                   n_taxa_expected, n_taxa_result))

    # Get taxa names
    fastspar_taxa <- rownames(cor_matrix)
    original_taxa <- rownames(otu_mat)

    # Create full-size matrices with zeros for filtered taxa
    full_cor <- matrix(0, nrow = n_taxa_expected, ncol = n_taxa_expected)
    full_p <- matrix(1, nrow = n_taxa_expected, ncol = n_taxa_expected)
    rownames(full_cor) <- colnames(full_cor) <- original_taxa
    rownames(full_p) <- colnames(full_p) <- original_taxa

    # Fill in the values for taxa that FastSpar kept
    for (i in seq_along(fastspar_taxa)) {
      for (j in seq_along(fastspar_taxa)) {
        taxa_i <- fastspar_taxa[i]
        taxa_j <- fastspar_taxa[j]
        full_cor[taxa_i, taxa_j] <- cor_matrix[i, j]
        full_p[taxa_i, taxa_j] <- p_matrix[i, j]
      }
    }

    cor_matrix <- full_cor
    p_matrix <- full_p
  }

  # Final validation - matrices should now be square
  if (nrow(cor_matrix) != ncol(cor_matrix)) {
    stop(sprintf("FastSpar correlation matrix is STILL not square after padding: %d x %d",
                nrow(cor_matrix), ncol(cor_matrix)))
  }

  if (nrow(p_matrix) != ncol(p_matrix)) {
    stop(sprintf("FastSpar p-value matrix is STILL not square after padding: %d x %d",
                nrow(p_matrix), ncol(p_matrix)))
  }

  # Ensure matrices are symmetric
  if (!isSymmetric(cor_matrix)) {
    warning("FastSpar correlation matrix is not symmetric - forcing symmetry")
    cor_matrix <- (cor_matrix + t(cor_matrix)) / 2
  }

  # Clean up temp files
  unlink(otu_file)
  unlink(cor_file)
  unlink(cov_file)
  unlink(pval_file)
  unlink(boot_dir, recursive = TRUE)

  cat("FastSpar completed!\n")

  return(list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    method = "fastspar"
  ))
}


#' Check if FastSpar is available
#' @return Logical indicating if FastSpar is installed
check_fastspar <- function() {
  # Check if fastspar is in PATH (container)
  fastspar_path <- Sys.which("fastspar")

  if (fastspar_path != "") {
    # Test if it runs
    result <- tryCatch({
      output <- system(sprintf("%s --version", fastspar_path),
             intern = TRUE, ignore.stderr = TRUE)
      TRUE
    }, error = function(e) FALSE)

    return(result)
  }

  # Fallback to container path
  fastspar_path <- "/opt/conda/bin/fastspar"
  if (file.exists(fastspar_path)) {
    return(TRUE)
  }

  return(FALSE)
}
