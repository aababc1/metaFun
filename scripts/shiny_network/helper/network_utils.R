# =============================================================================
# Network Analysis Helper Functions - SIF Container Version
# Based on ggClusterNet2 filtering and network computation strategies
# =============================================================================

library(phyloseq)
library(dplyr)
library(igraph)

# Load centrality metrics module - FIXED PATH FOR CONTAINER
source("/app/helper/centrality_metrics.R")

# Load fast correlation methods (fastspar, flashweave)
source("/app/helper/fast_correlation.R")

# Load working wrappers (override the basic versions in fast_correlation.R)
source("/app/helper/flashweave_wrapper.R")
source("/app/helper/fastspar_wrapper.R")  # Container-aware version

#' Aggregate phyloseq object to a specified taxonomic rank
#' @param ps phyloseq object
#' @param rank Taxonomic rank to aggregate to (e.g., "Genus", "Family", "Phylum")
#' @return Aggregated phyloseq object
aggregate_taxa <- function(ps, rank) {

  if (is.null(rank) || rank == "none") {
    message("No aggregation - using original data")
    return(ps)
  }

  message("Aggregating taxa to ", rank, " level...")

  # Get taxonomy table
  tax_mat <- as(tax_table(ps), "matrix")

  # Check if rank exists in taxonomy table
  if (!rank %in% colnames(tax_mat)) {
    # Try case-insensitive match
    rank_idx <- grep(paste0("^", rank, "$"), colnames(tax_mat), ignore.case = TRUE)
    if (length(rank_idx) == 0) {
      stop("Taxonomic rank '", rank, "' not found in taxonomy table. Available ranks: ",
           paste(colnames(tax_mat), collapse = ", "))
    }
    rank <- colnames(tax_mat)[rank_idx[1]]
  }

  # Get OTU table
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Get the taxonomic assignments at the specified rank
  rank_assignments <- tax_mat[, rank]

  # Remove NA or empty assignments
  valid_idx <- !is.na(rank_assignments) & rank_assignments != ""
  otu_mat <- otu_mat[valid_idx, , drop = FALSE]
  tax_mat <- tax_mat[valid_idx, , drop = FALSE]
  rank_assignments <- rank_assignments[valid_idx]

  message("  Original taxa: ", nrow(otu_mat))
  message("  Unique ", rank, " groups: ", length(unique(rank_assignments)))

  # Aggregate abundance data by summing all taxa within each taxonomic group
  unique_ranks <- unique(rank_assignments)
  aggregated_otu <- matrix(0, nrow = length(unique_ranks), ncol = ncol(otu_mat))
  rownames(aggregated_otu) <- unique_ranks
  colnames(aggregated_otu) <- colnames(otu_mat)

  # Create aggregated taxonomy table
  aggregated_tax <- matrix("", nrow = length(unique_ranks), ncol = ncol(tax_mat))
  rownames(aggregated_tax) <- unique_ranks
  colnames(aggregated_tax) <- colnames(tax_mat)

  for (r in unique_ranks) {
    # Sum abundances for all taxa in this group
    idx <- rank_assignments == r
    if (sum(idx) == 1) {
      aggregated_otu[r, ] <- otu_mat[idx, ]
      aggregated_tax[r, ] <- tax_mat[idx, ]
    } else {
      aggregated_otu[r, ] <- colSums(otu_mat[idx, , drop = FALSE])
      # Take the first occurrence for taxonomy (they should all be the same up to the selected rank)
      aggregated_tax[r, ] <- tax_mat[which(idx)[1], ]
    }
  }

  # Clear lower taxonomic ranks to avoid confusion
  # (e.g., if aggregating to Genus, clear Species; if aggregating to Family, clear Genus and Species)
  tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (rank %in% tax_ranks) {
    rank_idx <- which(tax_ranks == rank)
    if (rank_idx < length(tax_ranks)) {
      # Clear all ranks below the aggregation level
      lower_ranks <- tax_ranks[(rank_idx + 1):length(tax_ranks)]
      for (lower_rank in lower_ranks) {
        if (lower_rank %in% colnames(aggregated_tax)) {
          aggregated_tax[, lower_rank] <- NA
        }
      }
      message("  Cleared lower taxonomic ranks: ", paste(lower_ranks, collapse = ", "))
    }
  }

  # Create new phyloseq object with aggregated data
  ps_aggregated <- phyloseq(
    otu_table(aggregated_otu, taxa_are_rows = TRUE),
    tax_table(aggregated_tax),
    sample_data(ps)
  )

  message("  Aggregated to ", ntaxa(ps_aggregated), " ", rank, "-level taxa")

  return(ps_aggregated)
}

#' Apply filtering strategy to phyloseq object
#' @param ps phyloseq object
#' @param strategy Filtering strategy: "topn", "prev_topn", "ra_threshold"
#' @param top_n Number of top features to keep
#' @param prevalence Prevalence threshold (0-1)
#' @param ra_threshold Relative abundance threshold (0-1)
#' @param group Group column name for group-mean calculation
#' @return Filtered phyloseq object
apply_filtering <- function(ps, strategy = "topn", top_n = 200,
                          prevalence = 0.1, ra_threshold = 0.01,
                          group = NULL) {

  message("Applying filtering strategy: ", strategy)
  message("Original taxa: ", ntaxa(ps))

  if (strategy == "prev_topn") {
    # Step 1: Prevalence filtering
    ps <- filter_by_prevalence(ps, prevalence)
    message("After prevalence filter (", prevalence*100, "%): ", ntaxa(ps))

    # Step 2: Top-N by mean abundance
    ps <- filter_top_n(ps, top_n, group)
    message("After top-N filter (N=", top_n, "): ", ntaxa(ps))

  } else if (strategy == "topn") {
    # Top-N by mean abundance only
    ps <- filter_top_n(ps, top_n, group)
    message("After top-N filter (N=", top_n, "): ", ntaxa(ps))

  } else if (strategy == "ra_threshold") {
    # Filter by total relative abundance threshold
    ps <- filter_by_ra_threshold(ps, ra_threshold)
    message("After RA threshold filter (", ra_threshold*100, "%): ", ntaxa(ps))

  } else if (strategy == "prev_ra") {
    # Prevalence + Relative abundance filtering
    ps <- filter_by_prevalence(ps, prevalence)
    message("After prevalence filter (", prevalence*100, "%): ", ntaxa(ps))

    ps <- filter_by_ra_threshold(ps, ra_threshold)
    message("After RA threshold filter (", ra_threshold*100, "%): ", ntaxa(ps))
  }

  if (ntaxa(ps) == 0) {
    stop("No taxa remaining after filtering!")
  }

  return(ps)
}

#' Filter taxa by prevalence
#' @param ps phyloseq object
#' @param prevalence Minimum prevalence (proportion of samples, 0-1)
#' @return Filtered phyloseq object
filter_by_prevalence <- function(ps, prevalence = 0.1) {

  # Filter taxa present in at least X% of samples
  ps_filtered <- filter_taxa(ps, function(x) {
    sum(x > 0) / length(x) >= prevalence
  }, prune = TRUE)

  return(ps_filtered)
}

#' Filter top-N taxa by mean abundance
#' @param ps phyloseq object
#' @param top_n Number of top taxa to keep
#' @param group Group column for group-mean calculation (optional)
#' @return Filtered phyloseq object
filter_top_n <- function(ps, top_n = 200, group = NULL) {

  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  if (!is.null(group) && group %in% colnames(sample_data(ps))) {
    # Calculate group means for fair comparison
    # This ensures common feature set across groups
    group_means <- calculate_group_means(ps, group)
    overall_mean <- rowMeans(group_means, na.rm = TRUE)

  } else {
    # Simple overall mean
    overall_mean <- rowMeans(otu_mat, na.rm = TRUE)
  }

  # Sort by mean abundance
  top_taxa <- names(sort(overall_mean, decreasing = TRUE)[1:min(top_n, length(overall_mean))])

  # Prune to top taxa
  ps_filtered <- prune_taxa(top_taxa, ps)

  return(ps_filtered)
}

#' Calculate group-wise mean abundances
#' @param ps phyloseq object
#' @param group Group column name
#' @return Matrix of group means (taxa x groups)
calculate_group_means <- function(ps, group) {

  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  metadata <- as(sample_data(ps), "data.frame")
  groups <- unique(metadata[[group]])

  group_means <- matrix(0, nrow = nrow(otu_mat), ncol = length(groups))
  rownames(group_means) <- rownames(otu_mat)
  colnames(group_means) <- groups

  for (g in groups) {
    samples_in_group <- rownames(metadata)[metadata[[group]] == g]
    group_means[, g] <- rowMeans(otu_mat[, samples_in_group, drop = FALSE], na.rm = TRUE)
  }

  return(group_means)
}

#' Filter by total relative abundance threshold
#' @param ps phyloseq object
#' @param threshold Minimum total relative abundance across all samples
#' @return Filtered phyloseq object
filter_by_ra_threshold <- function(ps, threshold = 0.01) {

  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Calculate total relative abundance for each taxon
  total_abundance <- rowSums(otu_mat) / sum(otu_mat)

  # Keep taxa above threshold
  keep_taxa <- names(total_abundance[total_abundance >= threshold])

  ps_filtered <- prune_taxa(keep_taxa, ps)

  return(ps_filtered)
}

#' Run group-wise network analysis using ggClusterNet2 approach
#' @param ps Filtered phyloseq object
#' @param group Group column name
#' @param r.threshold Correlation threshold
#' @param p.threshold P-value threshold
#' @param method Correlation method ("spearman", "pearson", "kendall")
#' @return List of network results per group
run_network_analysis <- function(ps, group, r.threshold = 0.2,
                                p.threshold = 0.05, method = "spearman",
                                method_params = list()) {

  message("Running group-wise network analysis...")

  # Check if phyloseq has sample_data
  if (is.null(sample_data(ps, errorIfNULL = FALSE))) {
    stop("Phyloseq object has no sample_data. Cannot perform group-wise analysis.")
  }

  # Get sample data as data frame
  samp_data <- as(sample_data(ps), "data.frame")

  # Check if group column exists
  if (!group %in% colnames(samp_data)) {
    stop("Group column '", group, "' not found in sample data. Available columns: ",
         paste(colnames(samp_data), collapse = ", "))
  }

  # Get unique groups
  groups <- unique(as.character(samp_data[[group]]))

  message("Groups: ", paste(groups, collapse = ", "))

  if (length(groups) == 0) {
    stop("No groups found in column '", group, "'")
  }

  # Initialize results list
  network_results <- list()

  # Run network analysis for each group separately
  for (g in groups) {
    message("Processing group: ", g)

    tryCatch({
      # Subset phyloseq to current group using subset_samples
      # Use get_variable to safely access the column
      sample_filter <- get_variable(ps, group) == g
      ps_group <- prune_samples(sample_filter, ps)

      # Remove taxa with zero abundance in this group
      ps_group <- prune_taxa(taxa_sums(ps_group) > 0, ps_group)

      message("  Samples: ", nsamples(ps_group), ", Taxa: ", ntaxa(ps_group))

      # Calculate correlation network with method-specific parameters
      network <- calculate_correlation_network(
        ps_group,
        r.threshold = r.threshold,
        p.threshold = p.threshold,
        method = method,
        method_params = method_params
      )

      # Calculate basic network metrics
      metrics <- calculate_network_metrics(network$igraph)

      # Calculate comprehensive centrality metrics (Betweenness, Degree, Closeness, Zi-Pi, IVI)
      comprehensive_metrics <- get_comprehensive_network_metrics(network$igraph)

      # Store results
      network_results[[g]] <- list(
        phyloseq = ps_group,
        correlation_matrix = network$cor_matrix,
        p_matrix = network$p_matrix,
        igraph = network$igraph,
        metrics = metrics,
        node_metrics = comprehensive_metrics$node_metrics,
        zi_pi = comprehensive_metrics$zi_pi,
        ivi = comprehensive_metrics$ivi,
        robustness = comprehensive_metrics$robustness,
        communities = comprehensive_metrics$communities,
        n_samples = nsamples(ps_group),
        n_taxa = ntaxa(ps_group),
        method = method  # Store which method was used
      )

    }, error = function(e) {
      message("ERROR processing group ", g, ": ", e$message)
      message("Full error: ", toString(e))
      traceback()
      network_results[[g]] <- NULL
    })
  }

  # Add comparison metrics
  network_results$comparison <- compare_network_metrics(network_results)

  # brainGraph robustness analysis is now on-demand via "Run Stability Test" button
  # in Network Robustness tab (not automatic)

  return(network_results)
}

#' Calculate correlation network for a phyloseq object
#' @param ps phyloseq object
#' @param r.threshold Correlation threshold
#' @param p.threshold P-value threshold
#' @param method Correlation method ("spearman", "pearson", "fastspar", "flashweave")
#' @return List with correlation matrix, p-values, and igraph object
calculate_correlation_network <- function(ps, r.threshold = 0.2,
                                        p.threshold = 0.05,
                                        method = "spearman",
                                        method_params = list()) {

  # Get OTU table
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  n_taxa <- nrow(otu_mat)

  # Fast methods: fastspar and flashweave
  if (method == "fastspar") {
    message("  Using fastspar (fast C++ implementation)...")

    # Use user-provided parameters or defaults (proper NULL handling)
    iterations <- if (!is.null(method_params$iterations)) method_params$iterations else 50
    bootstraps <- if (!is.null(method_params$bootstraps)) method_params$bootstraps else 100
    threads <- if (!is.null(method_params$threads)) method_params$threads else 4

    message(sprintf("  Parameters: iterations=%d, bootstraps=%d, threads=%d",
                   iterations, bootstraps, threads))

    result <- run_fastspar(ps, iterations = iterations, bootstraps = bootstraps, threads = threads)
    cor_matrix <- result$cor_matrix
    p_matrix <- result$p_matrix

    # Count significant correlations before filtering
    n_all_edges <- sum(cor_matrix != 0 & !is.na(cor_matrix)) / 2
    n_significant <- sum(p_matrix < 0.05 & !is.na(p_matrix)) / 2

    # Convert to integer to avoid sprintf format issues
    n_all_edges <- as.integer(n_all_edges)
    n_significant <- as.integer(n_significant)

    message(sprintf("  FastSpar found %d correlations (%d significant at p<0.05)", n_all_edges, n_significant))

    message("  FastSpar completed successfully!")

  } else if (method == "flashweave") {
    message("  Using FlashWeave (Julia-based network inference)...")

    # Use user-provided parameters or defaults (proper NULL handling)
    sensitive <- if (!is.null(method_params$sensitive)) method_params$sensitive else TRUE
    heterogeneous <- if (!is.null(method_params$heterogeneous)) method_params$heterogeneous else TRUE
    max_k <- if (!is.null(method_params$max_k)) method_params$max_k else 3
    alpha <- if (!is.null(method_params$alpha)) method_params$alpha else 0.01
    n_threads <- if (!is.null(method_params$n_threads)) method_params$n_threads else 4
    FDR <- if (!is.null(method_params$FDR)) method_params$FDR else TRUE
    feed_forward <- if (!is.null(method_params$feed_forward)) method_params$feed_forward else TRUE

    message(sprintf("  Parameters: sensitive=%s, max_k=%d, alpha=%.3f, threads=%d, heterogeneous=%s, FDR=%s",
                   sensitive, max_k, alpha, n_threads, heterogeneous, FDR))

    result <- run_flashweave(ps,
                            sensitive = sensitive,
                            heterogeneous = heterogeneous,
                            max_k = max_k,
                            alpha = alpha,
                            conv = 0.01,
                            n_threads = n_threads,
                            FDR = FDR,
                            feed_forward = feed_forward)
    cor_matrix <- result$cor_matrix

    # Count detected edges
    n_edges <- sum(cor_matrix != 0) / 2  # Divide by 2 because matrix is symmetric
    n_edges <- as.integer(n_edges)  # Convert to integer for sprintf
    message(sprintf("  FlashWeave detected %d edges", n_edges))

    # FlashWeave returns SIGNED correlation weights (partial correlations)
    # These already represent correlation strength with sign (-1 to +1)
    # DO NOT convert to binary - preserve the signs!
    # Just use cor_matrix directly (it already has signed values)

    # FlashWeave already did significance testing internally
    # Set p-values to 0 for edges (significant), 1 for non-edges
    p_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
    p_matrix[cor_matrix != 0] <- 0
    rownames(p_matrix) <- rownames(cor_matrix)
    colnames(p_matrix) <- colnames(cor_matrix)

    # cor_matrix already has signed values - use it directly!

    if (n_edges == 0) {
      warning("  FlashWeave found 0 edges! Network will be empty. Try:\n",
              "    - Increase alpha (e.g., 0.05 or 0.1)\n",
              "    - Decrease max_k (e.g., 0 or 1)\n",
              "    - Set FDR=FALSE (not recommended for publication)\n",
              "    - Use 'Fast' mode instead of 'Sensitive'")
    }

    message("  FlashWeave completed successfully!")

  } else {
    # Standard correlation methods (Spearman, Pearson, Kendall)
    cor_matrix <- matrix(0, nrow = n_taxa, ncol = n_taxa)
    p_matrix <- matrix(1, nrow = n_taxa, ncol = n_taxa)
    rownames(cor_matrix) <- colnames(cor_matrix) <- rownames(otu_mat)
    rownames(p_matrix) <- colnames(p_matrix) <- rownames(otu_mat)

    # Calculate pairwise correlations
    for (i in 1:(n_taxa-1)) {
      for (j in (i+1):n_taxa) {
        test_result <- cor.test(otu_mat[i, ], otu_mat[j, ], method = method)
        cor_matrix[i, j] <- cor_matrix[j, i] <- test_result$estimate
        p_matrix[i, j] <- p_matrix[j, i] <- test_result$p.value
      }
    }
  }

  # Filter edges by threshold
  # For FlashWeave: skip filtering (it already did its own with FDR)
  # For fastspar and other correlation methods: apply r.threshold and p.threshold
  if (method == "flashweave") {
    # FlashWeave already filtered internally with FDR
    # Just use correlation matrix directly (already has only significant edges)
    adj_matrix <- cor_matrix
    n_edges_after_filter <- sum(adj_matrix != 0 & !is.na(adj_matrix)) / 2
    n_edges_after_filter <- as.integer(n_edges_after_filter)
    message(sprintf("  Using %s results directly (already filtered): %d edges",
                   method, n_edges_after_filter))
  } else {
    # Standard correlation methods: apply thresholds
    adj_matrix <- cor_matrix

    # Apply thresholds uniformly for all correlation methods
    # Note: SparCC/FastSpar produces valid correlation coefficients for compositional data
    # No special threshold adjustment needed - use same thresholds as other methods
    adj_matrix[abs(cor_matrix) < r.threshold | p_matrix > p.threshold] <- 0

    # Count edges after filtering
    n_edges_after_filter <- sum(adj_matrix != 0 & !is.na(adj_matrix)) / 2
    n_edges_after_filter <- as.integer(n_edges_after_filter)
    message(sprintf("  After filtering (|r|>%.2f, p<%.3f): %d edges remain",
                   r.threshold, p.threshold, n_edges_after_filter))

    if (n_edges_after_filter == 0) {
      warning("  No edges remain after filtering! Try:\n",
              "    - Decrease r.threshold (currently: ", r.threshold, ")\n",
              "    - Increase p.threshold (currently: ", p.threshold, ")\n",
              "    - Use less stringent method parameters")
    }
  }

  # Create igraph object
  # Use absolute values for graph structure (igraph needs non-negative for some algorithms)
  # BUT store original signed correlations as edge attribute for visualization
  adj_matrix_abs <- abs(adj_matrix)
  g <- graph_from_adjacency_matrix(adj_matrix_abs, mode = "undirected",
                                  weighted = TRUE, diag = FALSE)

  # IMPORTANT: Store original SIGNED correlations as separate edge attribute
  # This preserves positive/negative information for visualization
  if (ecount(g) > 0) {
    edge_list <- get.edgelist(g, names = TRUE)
    signed_weights <- sapply(1:nrow(edge_list), function(i) {
      adj_matrix[edge_list[i, 1], edge_list[i, 2]]
    })
    E(g)$signed_weight <- signed_weights  # Original signed correlation values

    # Debug: Show weight distribution for verification
    abs_weights <- E(g)$weight
    message(sprintf("  Edge weights: %d positive, %d negative",
                   sum(signed_weights > 0), sum(signed_weights < 0)))
    message(sprintf("  Weight range (|r|): %.3f ~ %.3f (mean: %.3f)",
                   min(abs_weights), max(abs_weights), mean(abs_weights)))
    message(sprintf("  |r| >= 0.3: %d edges, |r| >= 0.5: %d edges",
                   sum(abs_weights >= 0.3), sum(abs_weights >= 0.5)))
  }

  # Add taxonomy information to vertices if available
  if (!is.null(tax_table(ps, errorIfNULL = FALSE))) {
    tax <- as(tax_table(ps), "matrix")
    # Match node names to taxonomy
    node_names <- V(g)$name
    tax_matched <- tax[node_names, , drop = FALSE]

    # Create readable taxonomic labels (Genus species or Family if Genus missing)
    tax_labels <- apply(tax_matched, 1, function(x) {
      if (!is.na(x["Genus"]) && x["Genus"] != "") {
        genus <- x["Genus"]
        species <- if (!is.na(x["Species"]) && x["Species"] != "") {
          paste0(" ", gsub("^\\[|\\]$", "", x["Species"]))  # Remove brackets if present
        } else ""
        paste0(genus, species)
      } else if (!is.na(x["Family"]) && x["Family"] != "") {
        paste0(x["Family"], " (Family)")
      } else if (!is.na(x["Order"]) && x["Order"] != "") {
        paste0(x["Order"], " (Order)")
      } else {
        x["Phylum"]  # Fallback to Phylum
      }
    })

    V(g)$taxonomy <- tax_labels
    V(g)$genus <- tax_matched[, "Genus"]
    V(g)$species <- tax_matched[, "Species"]
    V(g)$family <- tax_matched[, "Family"]
  }

  # Remove isolated nodes
  isolated <- which(degree(g) == 0)
  if (length(isolated) > 0) {
    message(sprintf("  Removing %d isolated nodes", length(isolated)))
    g <- delete_vertices(g, isolated)
  }

  message(sprintf("  Final network: %d nodes, %d edges", as.integer(vcount(g)), as.integer(ecount(g))))

  return(list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    adj_matrix = adj_matrix,
    igraph = g
  ))
}

#' Calculate network topology metrics
#' @param g igraph object
#' @return Data frame of network metrics
calculate_network_metrics <- function(g) {

  if (vcount(g) == 0) {
    return(data.frame(
      n_nodes = 0,
      n_edges = 0,
      density = 0,
      avg_degree = 0,
      transitivity = 0,
      avg_path_length = NA,
      diameter = NA,
      modularity = NA
    ))
  }

  # Basic metrics
  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  density <- edge_density(g)
  avg_degree <- mean(degree(g))

  # Transitivity (clustering coefficient)
  transitivity_global <- transitivity(g, type = "global")

  # Path-based metrics (only for connected graphs)
  if (is_connected(g)) {
    avg_path <- mean_distance(g)
    diam <- diameter(g)
  } else {
    avg_path <- NA
    diam <- NA
  }

  # Modularity
  tryCatch({
    communities <- cluster_fast_greedy(g)
    mod <- modularity(communities)
  }, error = function(e) {
    mod <- NA
  })

  metrics <- data.frame(
    n_nodes = n_nodes,
    n_edges = n_edges,
    density = density,
    avg_degree = avg_degree,
    transitivity = transitivity_global,
    avg_path_length = avg_path,
    diameter = diam,
    modularity = mod
  )

  return(metrics)
}

#' Compare network metrics across groups
#' @param network_results List of network results
#' @return Data frame of comparison metrics
compare_network_metrics <- function(network_results) {

  message("=== COMPARE_NETWORK_METRICS DEBUG ===")

  # Remove non-group entries
  groups <- names(network_results)
  message("All names in network_results: ", paste(groups, collapse = ", "))
  groups <- groups[groups != "comparison"]
  message("Groups (without 'comparison'): ", paste(groups, collapse = ", "))

  # Extract metrics for each group
  metrics_list <- lapply(groups, function(g) {
    if (is.null(network_results[[g]])) {
      message("Group ", g, " is NULL, skipping")
      return(NULL)
    }

    message("Processing group: ", g)
    metrics <- network_results[[g]]$metrics
    message("  Metrics for ", g, ": ", paste(names(metrics), collapse = ", "))
    message("  Values: n_nodes=", metrics$n_nodes, ", n_edges=", metrics$n_edges)

    metrics$group <- g
    metrics$n_samples <- network_results[[g]]$n_samples
    metrics$n_taxa <- network_results[[g]]$n_taxa
    return(metrics)
  })

  # Combine into data frame
  message("Combining metrics from ", length(metrics_list), " groups")
  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- NULL

  message("Final comparison dataframe:")
  message("  Dimensions: ", nrow(metrics_df), " rows x ", ncol(metrics_df), " cols")
  message("  Columns: ", paste(colnames(metrics_df), collapse = ", "))
  if (nrow(metrics_df) > 0) {
    message("  First row: ", paste(metrics_df[1, ], collapse = ", "))
  }
  message("=====================================")

  return(metrics_df)
}

