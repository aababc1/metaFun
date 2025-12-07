# =============================================================================
# Network Centrality and Robustness Metrics
# Includes: Betweenness, Degree, Closeness, Zi-Pi, IVI, and Robustness tests
# =============================================================================

library(igraph)
library(dplyr)

#' Calculate comprehensive node centrality metrics
#' @param g igraph object
#' @return Data frame with node-level centrality metrics
calculate_node_centralities <- function(g) {

  # Handle empty graph
  if (vcount(g) == 0) {
    return(data.frame(
      node = character(0),
      degree = numeric(0),
      betweenness = numeric(0),
      closeness = numeric(0),
      eigenvector = numeric(0),
      pagerank = numeric(0),
      clustering_coefficient = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  node_names <- V(g)$name

  # Handle graph with nodes but no edges
  if (ecount(g) == 0) {
    return(data.frame(
      node = node_names,
      degree = rep(0, length(node_names)),
      betweenness = rep(0, length(node_names)),
      closeness = rep(NA, length(node_names)),
      eigenvector = rep(NA, length(node_names)),
      pagerank = rep(1/length(node_names), length(node_names)),
      clustering_coefficient = rep(0, length(node_names)),
      stringsAsFactors = FALSE
    ))
  }

  # Basic centrality metrics
  degree_cent <- degree(g)
  betweenness_cent <- betweenness(g, normalized = TRUE)

  # Closeness centrality (only for connected graphs)
  if (is_connected(g)) {
    closeness_cent <- closeness(g, normalized = TRUE)
  } else {
    # For disconnected graphs, calculate per component
    closeness_cent <- rep(NA, vcount(g))
    comps <- components(g)
    for (i in unique(comps$membership)) {
      comp_nodes <- which(comps$membership == i)
      if (length(comp_nodes) > 1) {
        subg <- induced_subgraph(g, comp_nodes)
        closeness_cent[comp_nodes] <- closeness(subg, normalized = TRUE)
      }
    }
  }

  # Eigenvector centrality
  eigen_cent <- tryCatch({
    eigen_centrality(g)$vector
  }, error = function(e) {
    rep(NA, vcount(g))
  })

  # PageRank
  pagerank_cent <- page_rank(g)$vector

  # Local clustering coefficient
  clustering_coef <- transitivity(g, type = "local")
  clustering_coef[is.nan(clustering_coef)] <- 0

  # Create data frame
  centralities <- data.frame(
    node = node_names,
    degree = degree_cent,
    betweenness = betweenness_cent,
    closeness = closeness_cent,
    eigenvector = eigen_cent,
    pagerank = pagerank_cent,
    clustering_coefficient = clustering_coef,
    stringsAsFactors = FALSE
  )

  return(centralities)
}

#' Calculate Zi-Pi roles (within-module connectivity and among-module connectivity)
#' @param g igraph object
#' @param communities Community structure object from igraph
#' @return Data frame with Zi and Pi values
calculate_zi_pi <- function(g, communities = NULL) {

  # Handle empty graph
  if (vcount(g) == 0) {
    return(data.frame(
      node = character(0),
      module = integer(0),
      Zi = numeric(0),
      Pi = numeric(0),
      role = character(0),
      stringsAsFactors = FALSE
    ))
  }

  node_names <- V(g)$name

  # Handle graph with nodes but no edges
  if (ecount(g) == 0) {
    return(data.frame(
      node = node_names,
      module = rep(1, length(node_names)),  # All in same "module"
      Zi = rep(0, length(node_names)),
      Pi = rep(0, length(node_names)),
      role = rep("Peripheral (no edges)", length(node_names)),
      stringsAsFactors = FALSE
    ))
  }

  # Detect communities if not provided
  if (is.null(communities)) {
    communities <- tryCatch({
      cluster_fast_greedy(g)
    }, error = function(e) {
      # Try alternative clustering for disconnected graphs
      tryCatch({
        cluster_louvain(g)
      }, error = function(e2) {
        return(NULL)
      })
    })
  }

  if (is.null(communities)) {
    return(data.frame(
      node = node_names,
      module = rep(1, length(node_names)),
      Zi = rep(NA, length(node_names)),
      Pi = rep(NA, length(node_names)),
      role = rep("Unable to detect communities", length(node_names)),
      stringsAsFactors = FALSE
    ))
  }

  membership <- as.numeric(membership(communities))  # Convert to numeric vector
  n_modules <- max(membership)

  # Calculate Zi (within-module degree z-score)
  Zi <- rep(0, vcount(g))

  for (m in 1:n_modules) {
    module_nodes <- which(membership == m)
    if (length(module_nodes) > 1) {
      # Get subgraph for this module
      subg <- induced_subgraph(g, module_nodes)
      module_degrees <- degree(subg)

      # Calculate z-score
      mean_deg <- mean(module_degrees)
      sd_deg <- sd(module_degrees)

      if (sd_deg > 0) {
        z_scores <- (module_degrees - mean_deg) / sd_deg
        Zi[module_nodes] <- z_scores
      }
    }
  }

  # Calculate Pi (among-module connectivity)
  Pi <- rep(0, vcount(g))

  for (i in 1:vcount(g)) {
    neighbors <- neighbors(g, i)
    if (length(neighbors) > 0) {
      neighbor_modules <- membership[neighbors]
      k_total <- length(neighbors)

      # Count connections to each module
      module_counts <- table(neighbor_modules)

      # Calculate Pi = 1 - sum((k_i / k)^2)
      pi_val <- 1 - sum((module_counts / k_total)^2)
      Pi[i] <- pi_val
    }
  }

  # Classify nodes into roles based on Zi and Pi
  roles <- rep("Peripheral", vcount(g))
  roles[Zi > 2.5] <- "Module hub"
  roles[Pi > 0.62] <- "Connector"
  roles[Zi > 2.5 & Pi > 0.62] <- "Network hub"

  zi_pi_df <- data.frame(
    node = node_names,
    module = membership,
    Zi = Zi,
    Pi = Pi,
    role = roles,
    stringsAsFactors = FALSE
  )

  return(zi_pi_df)
}

#' Calculate Integrated Value of Influence (IVI)
#' Uses the influential package for correct IVI calculation (Salavaty et al., 2020)
#' IVI = (DC' + LH') × ((NC' + CR') × (BC' + CI'))
#' @param g igraph object
#' @return Data frame with IVI scores and all component metrics
calculate_ivi <- function(g) {

  # Handle empty graph
  if (vcount(g) == 0) {
    return(data.frame(
      node = character(0),
      hubness_score = numeric(0),
      spreading_score = numeric(0),
      IVI = numeric(0),
      rank = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  node_names <- V(g)$name

  # Handle graph with nodes but no edges
  if (ecount(g) == 0) {
    return(data.frame(
      node = node_names,
      hubness_score = rep(0, length(node_names)),
      spreading_score = rep(0, length(node_names)),
      IVI = rep(0, length(node_names)),
      rank = seq_along(node_names),
      stringsAsFactors = FALSE
    ))
  }

  # REQUIRE influential package - no fallback
  if (!requireNamespace("influential", quietly = TRUE)) {
    stop("influential package is required for IVI calculation. Please install it: install.packages('influential')")
  }

  # Use influential package for correct IVI calculation
  # d=3: Collective influence at distance 3 (as per original paper)
  # scale="range": Min-Max normalization to [0,1]

  # Calculate IVI scores using influential package
  ivi_scores <- influential::ivi(
    graph = g,
    vertices = V(g),
    directed = FALSE,
    mode = "all",
    d = 3,
    scale = "range"
  )

  # Get individual component scores for detailed output
  # Hubness score = DC' + LH'
  hubness <- influential::hubness.score(g, scale = "range")

  # Spreading score = (NC' + CR') × (BC' + CI')
  spreading <- influential::spreading.score(g, d = 3, scale = "range")

  ivi_df <- data.frame(
    node = V(g)$name,
    hubness_score = hubness,
    spreading_score = spreading,
    IVI = ivi_scores,
    stringsAsFactors = FALSE
  )

  # Rank by IVI
  ivi_df <- ivi_df[order(ivi_df$IVI, decreasing = TRUE), ]
  ivi_df$rank <- 1:nrow(ivi_df)

  return(ivi_df)
}

#' Network robustness test - targeted vs random removal
#' @param g igraph object
#' @param metric Metric to use for targeted removal ("degree", "betweenness", "IVI")
#' @param steps Number of removal steps (default 10)
#' @return Data frame with robustness curves
calculate_robustness <- function(g, metric = "degree", steps = 10) {

  if (vcount(g) == 0) {
    return(data.frame())
  }

  n_nodes <- vcount(g)
  removal_fractions <- seq(0, 1, length.out = steps + 1)

  # Calculate node importance for targeted removal
  if (metric == "degree") {
    importance <- degree(g)
  } else if (metric == "betweenness") {
    importance <- betweenness(g, normalized = TRUE)
  } else if (metric == "IVI") {
    ivi_scores <- calculate_ivi(g)
    importance <- ivi_scores$IVI[match(V(g)$name, ivi_scores$node)]
  } else {
    stop("Unknown metric: ", metric)
  }

  # Order nodes by importance (descending)
  targeted_order <- order(importance, decreasing = TRUE)

  # Initialize results with proper column names
  results_targeted <- data.frame(
    fraction_removed = removal_fractions,
    largest_component = 0,
    avg_path_length = 0,
    efficiency = 0,
    strategy = "Targeted"
  )

  results_random <- data.frame(
    fraction_removed = removal_fractions,
    largest_component = 0,
    avg_path_length = 0,
    efficiency = 0,
    strategy = "Random"
  )

  # Targeted removal
  g_temp <- g
  for (i in 1:length(removal_fractions)) {
    frac <- removal_fractions[i]
    n_remove <- floor(frac * n_nodes)

    if (n_remove > 0 && n_remove < vcount(g)) {
      nodes_to_remove <- targeted_order[1:n_remove]
      g_temp <- delete_vertices(g, nodes_to_remove)
    } else {
      g_temp <- g
    }

    if (vcount(g_temp) > 0) {
      comps <- components(g_temp)
      results_targeted$largest_component[i] <- max(comps$csize) / n_nodes

      if (is_connected(g_temp)) {
        results_targeted$avg_path_length[i] <- average.path.length(g_temp)
        results_targeted$efficiency[i] <- 1 / results_targeted$avg_path_length[i]
      } else {
        results_targeted$avg_path_length[i] <- NA
        results_targeted$efficiency[i] <- 0
      }
    } else {
      results_targeted$largest_component[i] <- 0
      results_targeted$avg_path_length[i] <- NA
      results_targeted$efficiency[i] <- 0
    }
  }

  # Random removal (average over 10 runs)
  n_runs <- 10
  for (i in 1:length(removal_fractions)) {
    frac <- removal_fractions[i]
    n_remove <- floor(frac * n_nodes)

    lcc_sizes <- numeric(n_runs)
    apl_values <- numeric(n_runs)

    for (run in 1:n_runs) {
      g_temp <- g

      if (n_remove > 0 && n_remove < vcount(g)) {
        nodes_to_remove <- sample(1:vcount(g), n_remove)
        g_temp <- delete_vertices(g, nodes_to_remove)
      }

      if (vcount(g_temp) > 0) {
        comps <- components(g_temp)
        lcc_sizes[run] <- max(comps$csize) / n_nodes

        if (is_connected(g_temp)) {
          apl_values[run] <- average.path.length(g_temp)
        } else {
          apl_values[run] <- NA
        }
      } else {
        lcc_sizes[run] <- 0
        apl_values[run] <- NA
      }
    }

    results_random$largest_component[i] <- mean(lcc_sizes)
    results_random$avg_path_length[i] <- mean(apl_values, na.rm = TRUE)
    if (!is.na(results_random$avg_path_length[i]) && results_random$avg_path_length[i] > 0) {
      results_random$efficiency[i] <- 1 / results_random$avg_path_length[i]
    } else {
      results_random$efficiency[i] <- 0
    }
  }

  # Combine results
  robustness_df <- rbind(results_targeted, results_random)

  return(robustness_df)
}

#' Get comprehensive network analysis with all metrics
#' @param g igraph object
#' @return List with all centrality metrics, Zi-Pi, IVI, and robustness
get_comprehensive_network_metrics <- function(g) {

  if (vcount(g) == 0) {
    return(list(
      centralities = data.frame(),
      zi_pi = data.frame(),
      ivi = data.frame(),
      robustness = data.frame()
    ))
  }

  # Calculate communities for Zi-Pi
  communities <- tryCatch({
    cluster_fast_greedy(g)
  }, error = function(e) {
    NULL
  })

  # Calculate all metrics
  centralities <- calculate_node_centralities(g)
  zi_pi <- calculate_zi_pi(g, communities)
  ivi <- calculate_ivi(g)

  # Merge centralities with Zi-Pi and IVI
  node_metrics <- centralities %>%
    left_join(zi_pi, by = "node") %>%
    left_join(ivi[, c("node", "IVI", "rank")], by = "node")

  # Calculate robustness curves
  robustness_degree <- calculate_robustness(g, metric = "degree", steps = 20)
  robustness_betweenness <- calculate_robustness(g, metric = "betweenness", steps = 20)
  robustness_ivi <- calculate_robustness(g, metric = "IVI", steps = 20)

  # Store as named list for compatibility with module expectations
  robustness_list <- list(
    degree = robustness_degree,
    betweenness = robustness_betweenness,
    ivi = robustness_ivi,
    random = robustness_degree[robustness_degree$strategy == "Random", ]
  )

  return(list(
    node_metrics = node_metrics,
    zi_pi = zi_pi,
    ivi = ivi,
    robustness = robustness_list,
    communities = communities
  ))
}

# =============================================================================
# brainGraph-based Network Robustness Analysis
# =============================================================================

#' Run comprehensive network robustness analysis using brainGraph
#' @param igraph_list Named list of igraph objects (one per group)
#' @param n_iterations Number of random attack simulations (default = 100)
#' @param centrality_methods Character vector of attack strategies
#' @param optional_ivi_scores Named list of IVI score vectors (optional)
#' @param edge_removal Include edge-based attacks (default = TRUE)
#' @param output_dir Directory for saving plots (optional)
#' @param group_names Custom group names (optional)
#' @return List with robustness results and plots for each group
run_network_robustness_analysis <- function(
    igraph_list,
    n_iterations = 100,
    centrality_methods = c("random", "degree", "betweenness"),
    optional_ivi_scores = NULL,
    edge_removal = TRUE,
    output_dir = NULL,
    group_names = NULL
) {

  # Load brainGraph
  if (!requireNamespace("brainGraph", quietly = TRUE)) {
    stop("brainGraph package is required. Install with: install.packages('brainGraph')")
  }

  # Validate inputs
  if (!is.list(igraph_list) || length(igraph_list) == 0) {
    stop("igraph_list must be a non-empty named list of igraph objects")
  }

  # Use names from list if group_names not provided
  if (is.null(group_names)) {
    group_names <- names(igraph_list)
    if (is.null(group_names)) {
      group_names <- paste0("Group", seq_along(igraph_list))
      names(igraph_list) <- group_names
    }
  }

  cat("=============================================================================\n")
  cat("NETWORK ROBUSTNESS ANALYSIS WITH BRAINGRAPH\n")
  cat("=============================================================================\n\n")
  cat("Groups:", paste(group_names, collapse = ", "), "\n")
  cat("Iterations:", n_iterations, "\n")
  cat("Attack strategies:", paste(centrality_methods, collapse = ", "), "\n")
  if (!is.null(optional_ivi_scores)) cat("IVI-based attack: Enabled\n")
  if (edge_removal) cat("Edge removal: Enabled\n")
  cat("\n")

  # Initialize results storage
  results <- list()

  # Process each group
  for (group_name in group_names) {
    cat("Processing group:", group_name, "\n")

    g <- igraph_list[[group_name]]

    if (is.null(g) || vcount(g) == 0) {
      warning("  Skipping ", group_name, " - empty or NULL graph")
      next
    }

    # Ensure graph has vertex names (required by brainGraph)
    # brainGraph's random attack is very strict about vertex names
    original_names <- V(g)$name
    if (is.null(original_names)) {
      original_names <- as.character(1:vcount(g))
    }

    # Store original names for mapping back later
    V(g)$original_name <- original_names

    # Create ultra-simple vertex names for brainGraph compatibility
    # Use format: node_001, node_002, etc. (brainGraph accepts these)
    simple_names <- sprintf("node_%03d", 1:vcount(g))
    V(g)$name <- simple_names

    # Simplify graph (remove multiple edges and loops - required by brainGraph)
    if (!is.simple(g)) {
      cat("  Simplifying graph (removing multiple edges and loops)...\n")
      g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    }

    cat("  Nodes:", vcount(g), ", Edges:", ecount(g), "\n")
    cat("  Vertex names: ", head(V(g)$name, 3), "... (simplified for brainGraph)\n")

    # Verify degree distribution for targeted attacks
    degrees <- degree(g)
    cat("  Degree range:", min(degrees), "-", max(degrees), "\n")
    cat("  Top 10% degree threshold:", quantile(degrees, 0.9), "\n")

    # Initialize group results
    group_results <- list()
    all_data <- data.frame()

    # --- VERTEX-BASED ATTACKS ---

    # Random attack - CUSTOM IMPLEMENTATION
    # brainGraph's random attack has a known bug (Issue #25) that causes
    # "Invalid vertex names" errors due to parallel processing issues.
    # We implement our own random attack to avoid this bug.
    if ("random" %in% centrality_methods) {
      cat("  Running random attack (", n_iterations, " iterations) - Custom implementation\n", sep = "")
      cat("  [Info] Using custom implementation due to brainGraph bug #25\n")

      tryCatch({
        n_nodes <- vcount(g)
        orig_max <- max(components(g)$csize)

        # We'll test at each node removal (0%, 0.5%, 1%, ..., 100%)
        removal_fractions <- seq(0, 1, length.out = n_nodes + 1)
        n_fractions <- length(removal_fractions)

        # Storage for results from all iterations
        all_comp_sizes <- matrix(0, nrow = n_fractions, ncol = n_iterations)

        # Run N iterations of random removal
        for (iter in 1:n_iterations) {
          if (iter %% 20 == 0) {
            cat("    Iteration", iter, "/", n_iterations, "\n")
          }

          # Random ordering of vertices
          random_order <- sample(1:n_nodes)

          # Test each removal fraction
          for (i in 1:n_fractions) {
            frac <- removal_fractions[i]
            n_remove <- floor(frac * n_nodes)

            if (n_remove == 0) {
              # No removal - original component size
              all_comp_sizes[i, iter] <- orig_max
            } else if (n_remove >= n_nodes) {
              # All removed
              all_comp_sizes[i, iter] <- 0
            } else {
              # Remove vertices in random order
              g_temp <- delete_vertices(g, random_order[1:n_remove])
              if (vcount(g_temp) > 0) {
                comps <- components(g_temp)
                all_comp_sizes[i, iter] <- max(comps$csize)
              } else {
                all_comp_sizes[i, iter] <- 0
              }
            }
          }
        }

        # Average across all iterations
        avg_comp_sizes <- rowMeans(all_comp_sizes)

        # Normalize by original max component
        comp_pct <- avg_comp_sizes / orig_max

        random_df <- data.frame(
          fraction_removed = removal_fractions,
          largest_component = comp_pct,
          strategy = "Random",
          attack_type = "Vertex",
          group = group_name,
          stringsAsFactors = FALSE
        )

        cat("  ✓ Random attack SUCCESS - Generated", nrow(random_df), "data points\n")
        cat("    First point: removed=", random_df$fraction_removed[1],
            ", comp=", random_df$largest_component[1], "\n")
        cat("    Last point: removed=", tail(random_df$fraction_removed, 1),
            ", comp=", tail(random_df$largest_component, 1), "\n")

        all_data <- rbind(all_data, random_df)
        group_results$random_vertex <- random_df

      }, error = function(e) {
        cat("  ✗ Custom random attack FAILED:", conditionMessage(e), "\n")
        warning("Random attack failed for ", group_name, ": ", conditionMessage(e))
      })
    }

    # Degree-based attack
    if ("degree" %in% centrality_methods) {
      cat("  Running degree-based attack...\n")

      tryCatch({
        rob_degree <- brainGraph::robustness(
          g,
          type = "vertex",
          measure = "degree",
          N = 1  # N not used for targeted attacks
        )

        degree_df <- data.frame(
          fraction_removed = rob_degree$removed.pct,
          largest_component = rob_degree$comp.pct,
          strategy = "Degree",
          attack_type = "Vertex",
          group = group_name,
          stringsAsFactors = FALSE
        )

        cat("  Degree attack - Generated", nrow(degree_df), "data points\n")
        cat("  Degree attack - First point: removed=", degree_df$fraction_removed[1],
            ", comp=", degree_df$largest_component[1], "\n")
        cat("  Degree attack - At 10% removed: comp=",
            degree_df$largest_component[which.min(abs(degree_df$fraction_removed - 0.1))], "\n")
        cat("  Degree attack - Last point: removed=", tail(degree_df$fraction_removed, 1),
            ", comp=", tail(degree_df$largest_component, 1), "\n")

        # CRITICAL CHECK: Verify node removal order WITH ACTUAL NODE NAMES
        node_degrees <- degree(g)
        sorted_indices <- order(node_degrees, decreasing = TRUE)
        top5_indices <- sorted_indices[1:min(5, length(sorted_indices))]

        cat("  ============================================================\n")
        cat("  VERIFICATION: Top 5 nodes that will be removed FIRST\n")
        cat("  ============================================================\n")
        for (i in seq_along(top5_indices)) {
          idx <- top5_indices[i]
          simple_name <- V(g)$name[idx]
          original_name <- V(g)$original_name[idx]
          deg_value <- node_degrees[idx]
          cat(sprintf("  Rank %d: %s (original: %s) - Degree: %d\n",
                      i, simple_name, original_name, deg_value))
        }
        cat("  ============================================================\n")

        sorted_degrees <- sort(node_degrees, decreasing = TRUE)
        cat("  Degree value range - Top 5:", paste(head(sorted_degrees, 5), collapse=", "), "\n")
        cat("  Degree value range - Bottom 5:", paste(tail(sorted_degrees, 5), collapse=", "), "\n")

        # Check if brainGraph actually uses high-degree nodes first
        # brainGraph should remove nodes in order of vertex attribute 'degree'
        # Let's verify the decline rate
        early_comp <- degree_df$largest_component[which.min(abs(degree_df$fraction_removed - 0.1))]
        mid_comp <- degree_df$largest_component[which.min(abs(degree_df$fraction_removed - 0.5))]
        late_comp <- degree_df$largest_component[which.min(abs(degree_df$fraction_removed - 0.9))]

        early_decline <- 1.0 - early_comp  # Decline in first 10%
        late_decline <- mid_comp - late_comp  # Decline in last 40%

        cat("  Decline analysis:\n")
        cat("    - First 10% removed: ", round(early_decline * 100, 1), "% component lost\n", sep="")
        cat("    - At 50% removed: ", round((1.0 - mid_comp) * 100, 1), "% total lost\n", sep="")
        cat("    - Last 40% removed: ", round(late_decline * 100, 1), "% component lost\n", sep="")

        if (early_decline < 0.15) {
          cat("  WARNING: Only ", round(early_decline * 100, 1), "% decline after removing top 10% nodes!\n", sep="")
          cat("  Expected: >15% decline for scale-free networks with hubs\n")
          cat("  Possible causes:\n")
          cat("    1. Network has uniform degree distribution (no hubs)\n")
          cat("    2. High redundancy/mesh-like structure\n")
          cat("    3. brainGraph may not be removing highest-degree nodes first\n")
        } else {
          cat("  OK: Significant early decline indicates targeted attack working correctly\n")
        }

        all_data <- rbind(all_data, degree_df)
        group_results$degree_vertex <- degree_df

      }, error = function(e) {
        warning("  Failed to run degree attack: ", e$message)
      })
    }

    # Betweenness-based attack
    if ("betweenness" %in% centrality_methods) {
      cat("  Running betweenness-based attack...\n")

      tryCatch({
        rob_btwn <- brainGraph::robustness(
          g,
          type = "vertex",
          measure = "btwn.cent",
          N = 1
        )

        btwn_df <- data.frame(
          fraction_removed = rob_btwn$removed.pct,
          largest_component = rob_btwn$comp.pct,
          strategy = "Betweenness",
          attack_type = "Vertex",
          group = group_name,
          stringsAsFactors = FALSE
        )

        # VERIFICATION: Show top 5 nodes by betweenness
        node_btwn <- betweenness(g, normalized = TRUE)
        sorted_indices <- order(node_btwn, decreasing = TRUE)
        top5_indices <- sorted_indices[1:min(5, length(sorted_indices))]

        cat("  ============================================================\n")
        cat("  VERIFICATION: Top 5 nodes by Betweenness (will be removed FIRST)\n")
        cat("  ============================================================\n")
        for (i in seq_along(top5_indices)) {
          idx <- top5_indices[i]
          simple_name <- V(g)$name[idx]
          original_name <- V(g)$original_name[idx]
          btwn_value <- node_btwn[idx]
          cat(sprintf("  Rank %d: %s (original: %s) - Betweenness: %.4f\n",
                      i, simple_name, original_name, btwn_value))
        }
        cat("  ============================================================\n")

        all_data <- rbind(all_data, btwn_df)
        group_results$betweenness_vertex <- btwn_df

      }, error = function(e) {
        warning("  Failed to run betweenness attack: ", e$message)
      })
    }

    # IVI-based attack (custom implementation)
    if (!is.null(optional_ivi_scores) && group_name %in% names(optional_ivi_scores)) {
      cat("  Running IVI-based attack...\n")

      tryCatch({
        ivi_scores <- optional_ivi_scores[[group_name]]

        # Sort nodes by IVI (descending)
        node_order <- order(ivi_scores, decreasing = TRUE)
        n_nodes <- vcount(g)
        removal_fractions <- seq(0, 1, length.out = 50)

        # VERIFICATION: Show top 5 nodes by IVI
        top5_indices <- node_order[1:min(5, length(node_order))]

        cat("  ============================================================\n")
        cat("  VERIFICATION: Top 5 nodes by IVI (will be removed FIRST)\n")
        cat("  ============================================================\n")
        for (i in seq_along(top5_indices)) {
          idx <- top5_indices[i]
          simple_name <- V(g)$name[idx]
          original_name <- V(g)$original_name[idx]
          ivi_value <- ivi_scores[idx]
          cat(sprintf("  Rank %d: %s (original: %s) - IVI: %.2f\n",
                      i, simple_name, original_name, ivi_value))
        }
        cat("  ============================================================\n")

        lcc_sizes <- numeric(length(removal_fractions))

        for (i in seq_along(removal_fractions)) {
          frac <- removal_fractions[i]
          n_remove <- floor(frac * n_nodes)

          if (n_remove == 0) {
            lcc_sizes[i] <- 1.0
          } else if (n_remove >= n_nodes) {
            lcc_sizes[i] <- 0
          } else {
            g_temp <- delete_vertices(g, node_order[1:n_remove])
            if (vcount(g_temp) > 0) {
              comps <- components(g_temp)
              lcc_sizes[i] <- max(comps$csize) / n_nodes
            } else {
              lcc_sizes[i] <- 0
            }
          }
        }

        ivi_df <- data.frame(
          fraction_removed = removal_fractions,
          largest_component = lcc_sizes,
          strategy = "IVI",
          attack_type = "Vertex",
          group = group_name,
          stringsAsFactors = FALSE
        )

        all_data <- rbind(all_data, ivi_df)
        group_results$ivi_vertex <- ivi_df

      }, error = function(e) {
        warning("  Failed to run IVI attack: ", e$message)
      })
    }

    # --- EDGE-BASED ATTACKS ---

    if (edge_removal) {
      cat("  Running edge-based attacks...\n")

      # Random edge removal
      if ("random" %in% centrality_methods) {
        tryCatch({
          rob_edge_random <- brainGraph::robustness(
            g,
            type = "edge",
            measure = "random",
            N = n_iterations
          )

          edge_random_df <- data.frame(
            fraction_removed = rob_edge_random$removed.pct,
            largest_component = rob_edge_random$comp.pct,
            strategy = "Random",
            attack_type = "Edge",
            group = group_name,
            stringsAsFactors = FALSE
          )

          all_data <- rbind(all_data, edge_random_df)
          group_results$random_edge <- edge_random_df

        }, error = function(e) {
          warning("  Failed to run random edge removal: ", e$message)
        })
      }

      # Betweenness-based edge removal
      if ("betweenness" %in% centrality_methods) {
        tryCatch({
          rob_edge_btwn <- brainGraph::robustness(
            g,
            type = "edge",
            measure = "btwn.cent",
            N = 1
          )

          edge_btwn_df <- data.frame(
            fraction_removed = rob_edge_btwn$removed.pct,
            largest_component = rob_edge_btwn$comp.pct,
            strategy = "Betweenness",
            attack_type = "Edge",
            group = group_name,
            stringsAsFactors = FALSE
          )

          all_data <- rbind(all_data, edge_btwn_df)
          group_results$betweenness_edge <- edge_btwn_df

        }, error = function(e) {
          warning("  Failed to run betweenness edge removal: ", e$message)
        })
      }
    }

    # Store group results
    results[[group_name]] <- list(
      data = all_data,
      individual_attacks = group_results
    )

    cat("  Completed ", nrow(all_data), " attack scenarios\n\n", sep = "")
  }

  # Calculate summary metrics
  cat("Calculating robustness metrics...\n")

  summary_table <- data.frame()

  for (group_name in names(results)) {
    group_data <- results[[group_name]]$data

    if (nrow(group_data) == 0) next

    # Calculate area under curve (AUC) for each strategy
    strategies <- unique(group_data$strategy)
    attack_types <- unique(group_data$attack_type)

    for (strat in strategies) {
      for (atype in attack_types) {
        subset_data <- group_data[group_data$strategy == strat &
                                   group_data$attack_type == atype, ]

        if (nrow(subset_data) > 1) {
          # Calculate AUC using trapezoidal rule
          auc <- sum(diff(subset_data$fraction_removed) *
                    (head(subset_data$largest_component, -1) +
                     tail(subset_data$largest_component, -1))) / 2

          # R50: fraction removed to reach 50% of largest component
          r50_idx <- which(subset_data$largest_component <= 0.5)[1]
          r50 <- if (!is.na(r50_idx)) subset_data$fraction_removed[r50_idx] else 1.0

          summary_table <- rbind(summary_table, data.frame(
            group = group_name,
            strategy = strat,
            attack_type = atype,
            AUC = round(auc, 3),
            R50 = round(r50, 3),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  cat("Analysis complete!\n\n")

  # Return comprehensive results
  return(list(
    group_results = results,
    summary_table = summary_table,
    parameters = list(
      n_iterations = n_iterations,
      centrality_methods = centrality_methods,
      edge_removal = edge_removal,
      group_names = group_names
    )
  ))
}
