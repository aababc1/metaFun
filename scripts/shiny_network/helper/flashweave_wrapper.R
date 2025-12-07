# ==============================================================================
# FlashWeave Wrapper Functions
# ==============================================================================

#' Run FlashWeave network inference via Julia
#'
#' @param ps Phyloseq object
#' @param sensitive Use sensitive mode (TRUE=fine-grained, FALSE=fast) (default: TRUE)
#' @param heterogeneous Account for heterogeneous data (default: TRUE)
#' @param max_k Maximum conditioning set size (default: 3, use 0 to disable)
#' @param alpha Significance level (default: 0.01)
#' @param conv Convergence threshold (default: 0.01)
#' @param n_threads Number of Julia threads to use (default: 1)
#' @param FDR Apply False Discovery Rate correction (default: TRUE)
#' @param feed_forward Enable feed-forward heuristic (default: TRUE)
#' @param julia_path Path to Julia binary (default: auto-detect)
#' @return List containing correlation matrix and adjacency matrix
run_flashweave <- function(ps, sensitive = TRUE, heterogeneous = TRUE,
                           max_k = 3, alpha = 0.01, conv = 0.01,
                           n_threads = 1, FDR = TRUE, feed_forward = TRUE,
                           julia_path = NULL) {

  # Check if microeco is available
  if (!requireNamespace("microeco", quietly = TRUE)) {
    stop("microeco package is required for FlashWeave. Install it first.")
  }

  # Auto-detect Julia if not provided
  if (is.null(julia_path)) {
    # Check common Julia paths in order of priority
    julia_candidates <- c(
      "/opt/julia/julia-1.10.9/bin/julia",  # Container path (primary)
      "/opt/julia/julia-1.10.7/bin/julia",  # Container path (alternative)
      "/usr/local/bin/julia",               # Common system location
      Sys.which("julia")                    # System PATH
    )

    julia_path <- NULL
    for (candidate in julia_candidates) {
      if (candidate != "" && file.exists(candidate)) {
        julia_path <- candidate
        break
      }
    }

    if (is.null(julia_path)) {
      stop("Julia not found. Please install Julia or provide path.")
    }
  }

  # Create temporary directory
  tmp_dir <- tempdir()
  otu_file <- file.path(tmp_dir, "flashweave_input.tsv")
  network_file <- file.path(tmp_dir, "flashweave_network.edgelist")
  weight_matrix_file <- file.path(tmp_dir, "flashweave_weights.tsv")
  taxa_names_file <- file.path(tmp_dir, "flashweave_taxa_names.txt")
  julia_script <- file.path(tmp_dir, "run_flashweave.jl")

  # Prepare OTU table (samples as rows, taxa as columns for FlashWeave)
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # FlashWeave with normalize=true handles relative abundance directly!
  # DO NOT convert to pseudo-counts - it breaks the analysis
  max_val <- max(otu_mat, na.rm = TRUE)
  if (max_val <= 1.0) {
    message("  Using relative abundance data directly (FlashWeave will normalize)...")
    # Keep as-is, FlashWeave's normalize=true will handle it
  } else {
    message("  Data appears to be in count format...")
  }

  # Write OTU table
  otu_df <- data.frame(
    Sample = rownames(otu_mat),
    otu_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(otu_df, otu_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)

  # Get taxa names for later
  taxa_names <- colnames(otu_mat)
  n_taxa <- length(taxa_names)

  # Create Julia script with threading support
  julia_code <- sprintf('
# Set number of threads
using Base.Threads
println("Julia threads available: ", Threads.nthreads())

using FlashWeave
using DelimitedFiles

# Read data
data, headers = readdlm("%s", \'\\t\', header=true)
otu_data = Float64.(data[:, 2:end])  # Skip sample names column
taxa_names = String.(headers[2:end])

# Run FlashWeave with parameters
println("Running FlashWeave...")
println("  sensitive: %s")
println("  heterogeneous: %s")
println("  max_k: %d")
println("  alpha: %f")
println("  conv: %f")
println("  FDR: %s")
println("  feed_forward: %s")
println("  normalize: true")
println("  n_obs_min: -1 (auto)")

network = learn_network(otu_data,
    sensitive=%s,
    heterogeneous=%s,
    max_k=%d,
    alpha=%f,
    conv=%f,
    FDR=%s,
    feed_forward=%s,
    normalize=true,
    n_obs_min=-1,
    header=taxa_names
)

# Save network edgelist
save_network("%s", network)

# Also save the weight matrix with signs
G = graph(network)
W = G.weights

# Get the actual variable names that FlashWeave kept (after filtering)
kept_var_names = names(network)

# Convert sparse matrix to dense for saving
W_dense = Matrix(W)

# Save weight matrix
writedlm("%s", W_dense, \'\\t\')

# Save the taxa names that were kept
writedlm("%s", kept_var_names)

println("FlashWeave completed!")
println("Input variables: ", length(taxa_names))
println("Kept variables: ", length(kept_var_names))
println("Weight matrix: ", size(W_dense))
println("Positive weights: ", sum(W_dense .> 0), ", Negative weights: ", sum(W_dense .< 0))
',
    otu_file,
    ifelse(sensitive, "true", "false"),
    ifelse(heterogeneous, "true", "false"),
    max_k,
    alpha,
    conv,
    ifelse(FDR, "true", "false"),
    ifelse(feed_forward, "true", "false"),
    ifelse(sensitive, "true", "false"),
    ifelse(heterogeneous, "true", "false"),
    max_k,
    alpha,
    conv,
    ifelse(FDR, "true", "false"),
    ifelse(feed_forward, "true", "false"),
    network_file,
    weight_matrix_file,
    taxa_names_file
  )

  writeLines(julia_code, julia_script)

  cat("Running FlashWeave via Julia...\n")
  cat(sprintf("  Threads: %d\n", n_threads))
  cat(sprintf("  Mode: %s\n", ifelse(sensitive, "Sensitive (fine-grained)", "Fast")))

  # Run Julia script with threading support
  # Set JULIA_NUM_THREADS environment variable
  old_threads <- Sys.getenv("JULIA_NUM_THREADS")
  Sys.setenv(JULIA_NUM_THREADS = n_threads)

  cmd <- sprintf("%s %s", julia_path, julia_script)
  result <- system(cmd, intern = TRUE)

  # Restore old thread setting
  if (old_threads != "") {
    Sys.setenv(JULIA_NUM_THREADS = old_threads)
  } else {
    Sys.unsetenv("JULIA_NUM_THREADS")
  }

  if (!file.exists(weight_matrix_file)) {
    stop("FlashWeave failed. Check Julia installation and FlashWeave package.")
  }

  # Read the taxa names that FlashWeave actually kept (after filtering)
  kept_taxa_names <- readLines(taxa_names_file)

  # Read weight matrix directly (includes signs!)
  cor_matrix <- as.matrix(read.table(weight_matrix_file, sep = "\t", header = FALSE))

  # Check dimensions match
  if (nrow(cor_matrix) != length(kept_taxa_names) || ncol(cor_matrix) != length(kept_taxa_names)) {
    stop(sprintf("Dimension mismatch: matrix is %dx%d but have %d taxa names",
                 nrow(cor_matrix), ncol(cor_matrix), length(kept_taxa_names)))
  }

  rownames(cor_matrix) <- colnames(cor_matrix) <- kept_taxa_names

  # If some taxa were filtered out, expand matrix to match original taxa
  if (length(kept_taxa_names) < length(taxa_names)) {
    message(sprintf("  FlashWeave filtered %d taxa (kept %d of %d)",
                    length(taxa_names) - length(kept_taxa_names),
                    length(kept_taxa_names),
                    length(taxa_names)))

    # Create full matrix with zeros for filtered taxa
    full_matrix <- matrix(0, nrow = n_taxa, ncol = n_taxa)
    rownames(full_matrix) <- colnames(full_matrix) <- taxa_names

    # Fill in the values for kept taxa
    idx <- match(kept_taxa_names, taxa_names)
    full_matrix[idx, idx] <- cor_matrix

    cor_matrix <- full_matrix
  }

  # Also read edge list for compatibility
  edges <- tryCatch({
    read.table(network_file, sep = "\t", header = FALSE,
               stringsAsFactors = FALSE, comment.char = "#")
  }, error = function(e) {
    data.frame(from = character(0), to = character(0), weight = numeric(0))
  })

  # Process edge list but now we can get real weights from the matrix
  if (nrow(edges) > 0 && ncol(edges) >= 2) {
    # Create edge dataframe with actual weights from correlation matrix
    edge_list <- data.frame(from = character(0), to = character(0), weight = numeric(0))

    for (i in seq_len(nrow(edges))) {
      from_val <- as.character(edges[i, 1])
      to_val <- as.character(edges[i, 2])

      # Get weight from correlation matrix
      if (from_val %in% taxa_names && to_val %in% taxa_names) {
        from_idx <- which(taxa_names == from_val)
        to_idx <- which(taxa_names == to_val)
        weight <- cor_matrix[from_idx, to_idx]

        edge_list <- rbind(edge_list, data.frame(
          from = from_val,
          to = to_val,
          weight = weight
        ))
      }
    }
    edges <- edge_list
  } else {
    edges <- data.frame(from = character(0), to = character(0), weight = numeric(0))
  }

  # Clean up
  unlink(otu_file)
  unlink(network_file)
  unlink(weight_matrix_file)
  unlink(taxa_names_file)
  unlink(julia_script)

  # Report edge statistics
  n_positive <- sum(cor_matrix > 0 & row(cor_matrix) < col(cor_matrix))
  n_negative <- sum(cor_matrix < 0 & row(cor_matrix) < col(cor_matrix))
  n_total <- n_positive + n_negative

  cat("FlashWeave completed!\n")
  cat(sprintf("  Total edges: %d\n", n_total))
  cat(sprintf("  Positive edges: %d (%.1f%%)\n", n_positive, 100 * n_positive / max(n_total, 1)))
  cat(sprintf("  Negative edges: %d (%.1f%%)\n", n_negative, 100 * n_negative / max(n_total, 1)))

  # Debug: Verify matrix has correct values
  cat("  Matrix value range: [", min(cor_matrix), ", ", max(cor_matrix), "]\n", sep = "")
  cat("  Total non-zero: ", sum(cor_matrix != 0), "\n", sep = "")

  return(list(
    cor_matrix = cor_matrix,
    adjacency_matrix = (cor_matrix != 0) * 1,  # Binary adjacency
    edges = edges,
    method = "flashweave"
  ))
}


#' Check if FlashWeave is available
#' @return Logical indicating if FlashWeave is installed
check_flashweave <- function() {
  # Check common Julia paths in order of priority
  julia_candidates <- c(
    "/opt/julia/julia-1.10.9/bin/julia",  # Container path (primary)
    "/opt/julia/julia-1.10.7/bin/julia",  # Container path (alternative)
    "/usr/local/bin/julia",               # Common system location
    Sys.which("julia")                    # System PATH
  )

  julia_path <- NULL
  for (candidate in julia_candidates) {
    if (candidate != "" && file.exists(candidate)) {
      julia_path <- candidate
      break
    }
  }

  if (is.null(julia_path)) {
    return(FALSE)
  }

  # Test if FlashWeave package is available
  result <- tryCatch({
    cmd <- sprintf('%s -e "using FlashWeave; println(\\"OK\\")"', julia_path)
    output <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    any(grepl("OK", output))
  }, error = function(e) FALSE)

  return(result)
}
