# =============================================================================
# Fast Correlation Methods for Network Analysis
# Supports: fastspar (C++ SparCC) and FlashWeave (Julia)
# =============================================================================

library(phyloseq)

#' Run fastspar for fast SparCC correlation
#' @param ps phyloseq object
#' @param output_dir Directory for temporary files
#' @param iterations Number of iterations (default: 50, faster than SpiecEasi's 20)
#' @param bootstraps Number of bootstrap iterations for p-values (default: 100)
#' @param threads Number of threads (default: 1)
#' @return List with correlation matrix and p-value matrix
run_fastspar <- function(ps, output_dir = tempdir(),
                         iterations = 50, bootstraps = 100, threads = 1) {

  # Check if fastspar is available
  fastspar_path <- Sys.which("fastspar")
  if (fastspar_path == "") {
    stop("fastspar not found. Install with: conda install -c bioconda fastspar")
  }

  message("Using fastspar for SparCC correlation...")
  message("  Iterations: ", iterations)
  message("  Bootstraps: ", bootstraps)
  message("  Threads: ", threads)

  # Create temporary directory
  temp_dir <- file.path(output_dir, paste0("fastspar_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  boot_dir <- file.path(temp_dir, "bootstrap")
  dir.create(boot_dir, showWarnings = FALSE)

  # Get OTU table (taxa as rows, samples as columns)
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Write OTU table in TSV format (fastspar format)
  otu_file <- file.path(temp_dir, "otu_table.tsv")
  write.table(otu_mat, otu_file, sep = "\t", quote = FALSE,
              col.names = NA, row.names = TRUE)

  # Output files
  cor_file <- file.path(temp_dir, "correlation.tsv")
  cov_file <- file.path(temp_dir, "covariance.tsv")
  pval_file <- file.path(temp_dir, "pvalues.tsv")

  # Run fastspar
  cmd <- sprintf(
    "fastspar --otu_table %s --correlation %s --covariance %s --iterations %d --threads %d",
    otu_file, cor_file, cov_file, iterations, threads
  )

  message("Running fastspar...")
  system(cmd, intern = FALSE)

  # Bootstrap for p-values
  message("Generating bootstrap samples...")
  boot_prefix <- file.path(boot_dir, "boot_")
  cmd_boot <- sprintf(
    "fastspar_bootstrap --otu_table %s --number %d --prefix %s --threads %d",
    otu_file, bootstraps, boot_prefix, threads
  )
  system(cmd_boot, intern = FALSE)

  # Run fastspar on bootstrap samples
  message("Computing bootstrap correlations...")
  for (i in 0:(bootstraps-1)) {
    boot_otu <- sprintf("%s%d.tsv", boot_prefix, i)
    boot_cor <- sprintf("%s%d_cor.tsv", boot_prefix, i)
    boot_cov <- sprintf("%s%d_cov.tsv", boot_prefix, i)

    cmd_boot_cor <- sprintf(
      "fastspar --otu_table %s --correlation %s --covariance %s --iterations %d --threads 1",
      boot_otu, boot_cor, boot_cov, iterations
    )
    system(cmd_boot_cor, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }

  # Calculate p-values
  message("Calculating p-values...")
  cmd_pval <- sprintf(
    "fastspar_pvalues --otu_table %s --correlation %s --prefix %s --permutations %d --outfile %s --threads %d",
    otu_file, cor_file, boot_prefix, bootstraps, pval_file, threads
  )
  system(cmd_pval, intern = FALSE)

  # Read results
  cor_matrix <- as.matrix(read.table(cor_file, sep = "\t", header = TRUE,
                                     row.names = 1, check.names = FALSE))
  p_matrix <- as.matrix(read.table(pval_file, sep = "\t", header = TRUE,
                                   row.names = 1, check.names = FALSE))

  # Clean up
  unlink(temp_dir, recursive = TRUE)

  message("fastspar completed successfully!")

  return(list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    method = "fastspar"
  ))
}

#' Run FlashWeave for fast network inference
#' @param ps phyloseq object
#' @param sensitive Boolean, use sensitive mode (slower but more accurate)
#' @param heterogeneous Boolean, account for heterogeneous data
#' @param max_k Maximum number of neighbors (default: 3)
#' @param alpha Significance level (default: 0.01)
#' @param conv Convergence threshold (default: 0.01)
#' @return List with correlation matrix (partial correlations)
run_flashweave <- function(ps, sensitive = TRUE, heterogeneous = TRUE,
                           max_k = 3, alpha = 0.01, conv = 0.01) {

  # Check if JuliaCall is available
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("JuliaCall package required. Install with: install.packages('JuliaCall')")
  }

  library(JuliaCall)

  message("Using FlashWeave for network inference...")
  message("  Sensitive mode: ", sensitive)
  message("  Heterogeneous: ", heterogeneous)
  message("  Max neighbors: ", max_k)

  # Initialize Julia
  julia_setup()

  # Load FlashWeave
  julia_library("FlashWeave")

  # Get OTU table (samples as rows, taxa as columns for FlashWeave)
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Convert to Julia array
  julia_assign("otu_data", otu_mat)

  # Run FlashWeave
  message("Running FlashWeave network inference...")
  julia_command(sprintf(
    "network_result = learn_network(otu_data, sensitive=%s, heterogeneous=%s, max_k=%d, alpha=%f, conv=%f)",
    ifelse(sensitive, "true", "false"),
    ifelse(heterogeneous, "true", "false"),
    max_k, alpha, conv
  ))

  # Get edge list
  edge_list <- julia_eval("collect(edges(network_result))")

  # Convert to correlation matrix format
  n_taxa <- ntaxa(ps)
  taxa_names <- taxa_names(ps)

  # Initialize matrix with zeros
  cor_matrix <- matrix(0, nrow = n_taxa, ncol = n_taxa)
  rownames(cor_matrix) <- colnames(cor_matrix) <- taxa_names

  # Fill in edges (FlashWeave returns binary edges, we set them to 1)
  if (nrow(edge_list) > 0) {
    for (i in 1:nrow(edge_list)) {
      from <- edge_list[i, 1]
      to <- edge_list[i, 2]
      cor_matrix[from, to] <- 1
      cor_matrix[to, from] <- 1
    }
  }

  message("FlashWeave completed successfully!")
  message("  Edges found: ", nrow(edge_list))

  return(list(
    cor_matrix = cor_matrix,
    p_matrix = NULL,  # FlashWeave doesn't provide p-values
    edge_list = edge_list,
    method = "flashweave",
    n_edges = nrow(edge_list)
  ))
}

#' Wrapper function to call appropriate fast correlation method
#' @param ps phyloseq object
#' @param method Method to use: "fastspar" or "flashweave"
#' @param ... Additional parameters passed to specific methods
#' @return List with correlation and p-value matrices
fast_correlation <- function(ps, method = "fastspar", ...) {

  method <- tolower(method)

  if (method == "fastspar") {
    return(run_fastspar(ps, ...))
  } else if (method == "flashweave") {
    return(run_flashweave(ps, ...))
  } else {
    stop("Unknown method: ", method, ". Use 'fastspar' or 'flashweave'")
  }
}
