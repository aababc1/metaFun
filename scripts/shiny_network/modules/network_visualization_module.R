# =============================================================================
# Network Visualization Module
# Interactive network graph visualization using igraph and plotly
# =============================================================================

library(shiny)
library(shinyjs)
library(bslib)
library(plotly)
library(igraph)
library(ggplot2)
library(dplyr)

#' Network Visualization UI
network_visualization_UI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),

    card(
      full_screen = TRUE,
      class = "main-card",
      card_header(
        style = "background-color: #2FA4E7; color: white;",
        "Network Visualization"
      ),
      card_body(
        layout_sidebar(
          sidebar = sidebar(
            title = "Visualization Settings",
            width = 300,

            # View mode selection
            radioButtons(
              ns("view_mode"),
              "View Mode:",
              choices = c(
                "Compare Side-by-Side" = "compare",
                "Single Group" = "single"
              ),
              selected = "compare"
            ),

            # Group selection (only shown in single mode)
            conditionalPanel(
              condition = "input.view_mode == 'single'",
              ns = ns,
              selectInput(
                ns("selected_group"),
                "Select Group to Visualize:",
                choices = NULL
              )
            ),

            hr(),

            # Layout algorithm
            selectInput(
              ns("layout_algorithm"),
              "Network Layout:",
              choices = c(
                "Fruchterman-Reingold" = "fr",
                "Kamada-Kawai" = "kk",
                "Circle" = "circle",
                "Grid" = "grid",
                "Sphere" = "sphere",
                "Random" = "random"
              ),
              selected = "fr"
            ),

            # Node coloring
            selectInput(
              ns("node_color_by"),
              "Color Nodes By:",
              choices = c(
                "Degree" = "degree",
                "Betweenness" = "betweenness",
                "Module" = "module",
                "Taxonomy" = "taxonomy"
              ),
              selected = "degree"
            ),

            # Node size
            selectInput(
              ns("node_size_by"),
              "Size Nodes By:",
              choices = c(
                "Degree" = "degree",
                "Betweenness" = "betweenness",
                "IVI" = "ivi",
                "Uniform" = "uniform"
              ),
              selected = "degree"
            ),

            # Edge opacity control
            sliderInput(
              ns("edge_opacity"),
              "Edge Opacity (transparency):",
              min = 0.1,
              max = 1.0,
              value = 0.3,
              step = 0.1
            ) |>
              tooltip("Lower values make edges more transparent (recommended: 0.2-0.4 for dense networks)"),

            # Edge filtering by correlation strength
            sliderInput(
              ns("min_edge_weight"),
              "Filter by Correlation Strength (|r|):",
              min = 0,
              max = 1,
              value = 0,
              step = 0.05
            ) |>
              tooltip("Filter edges by absolute correlation value. Set higher to show only strong correlations (e.g., |r| ≥ 0.5). This filters the existing network without redrawing it."),

            hr(),

            # [COMMENTED OUT] Differential Network Analysis - To be implemented later
            # conditionalPanel(
            #   condition = "input.view_mode == 'compare'",
            #   ns = ns,
            #   h5("Differential Network Analysis", style = "margin-top: 15px; color: #28a745;"),
            #   p("Compare network differences between groups using NetCoMi",
            #     style = "font-size: 12px; color: #666;"),
            #
            #   actionButton(
            #     ns("run_diffnet"),
            #     "Run Differential Analysis",
            #     icon = icon("chart-line"),
            #     class = "btn-success",
            #     style = "width: 100%; margin-bottom: 10px;"
            #   ),
            #
            #   selectInput(
            #     ns("diffnet_method"),
            #     "Test Method:",
            #     choices = c(
            #       "Fisher's Exact Test" = "fisherTest",
            #       "Permutation Test" = "permute"
            #     ),
            #     selected = "fisherTest"
            #   ) |>
            #     tooltip("Fisher's test: Fast, assumes independence. Permutation: More accurate but slower."),
            #
            #   selectInput(
            #     ns("diffnet_adjust"),
            #     "Multiple Testing Adjustment:",
            #     choices = c(
            #       "Local FDR" = "lfdr",
            #       "Benjamini-Hochberg (FDR)" = "BH",
            #       "Benjamini-Yekutieli (BY)" = "BY",
            #       "Bonferroni" = "bonferroni",
            #       "None" = "none"
            #     ),
            #     selected = "lfdr"
            #   ) |>
            #     tooltip("LFDR: Recommended for network data. BH/BY: Standard FDR control."),
            #
            #   hr()
            # ),

            # Download
            downloadButton(ns("download_network"), "Download Network Plot"),
            br(), br(),
            downloadButton(ns("download_graphml"), "Download GraphML")
          ),

          # Main panel - network plot
          div(
            id = ns("analysis_message"),
            style = "color: red; font-weight: bold; margin: 20px;",
            icon("exclamation-triangle"),
            " Run network analysis first"
          ),

          # Comparison view (side-by-side)
          conditionalPanel(
            condition = "input.view_mode == 'compare'",
            ns = ns,
            layout_columns(
              col_widths = c(6, 6),
              card(
                full_screen = TRUE,
                fill = FALSE,
                card_header(textOutput(ns("group1_title"))),
                plotlyOutput(ns("network_plot_group1"), height = "900px")
              ),
              card(
                full_screen = TRUE,
                fill = FALSE,
                card_header(textOutput(ns("group2_title"))),
                plotlyOutput(ns("network_plot_group2"), height = "900px")
              )
            )
          ),

          # Single view
          conditionalPanel(
            condition = "input.view_mode == 'single'",
            ns = ns,
            plotlyOutput(ns("network_plot"), height = "1000px")
          ),

          br(),

          # [COMMENTED OUT] Differential Network Analysis Section - To be implemented later
          # conditionalPanel(
          #   condition = "input.view_mode == 'compare'",
          #   ns = ns,
          #
          #   h4("Differential Network Analysis"),
          #
          #   # Warning message - run analysis
          #   div(
          #     id = ns("diffnet_warning"),
          #     style = "color: #856404; background-color: #fff3cd; border: 1px solid #ffeeba; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
          #     icon("exclamation-triangle"),
          #     " Click 'Run Differential Analysis' button in the left sidebar to compare networks between groups"
          #   ),
          #
          #   # Summary - simple layout like Network Statistics
          #   layout_column_wrap(
          #     width = "200px",
          #     fill = FALSE,
          #     value_box(
          #       title = "Test Method",
          #       value = textOutput(ns("diffnet_method_display")),
          #       theme = "success",
          #       height = "100px"
          #     ),
          #     value_box(
          #       title = "Adjustment",
          #       value = textOutput(ns("diffnet_adjust_display")),
          #       theme = "info",
          #       height = "100px"
          #     ),
          #     value_box(
          #       title = "Diff. Nodes (p<0.05)",
          #       value = textOutput(ns("diffnet_n_diff_nodes")),
          #       theme = "warning",
          #       height = "100px"
          #     ),
          #     value_box(
          #       title = "Diff. Edges (p<0.05)",
          #       value = textOutput(ns("diffnet_n_diff_edges")),
          #       theme = "danger",
          #       height = "100px"
          #     )
          #   ),
          #
          #   br(),
          #
          #   # Differential plot - same structure as Network Analysis tab
          #   h5("Differential Components by P-value Threshold"),
          #   div(
          #     style = "margin-bottom: 30px; padding: 20px; background-color: #ffffff; border: 1px solid #dee2e6; border-radius: 5px;",
          #     plotlyOutput(ns("diffnet_plot"), height = "400px")
          #   ),
          #
          #   br(),
          #
          #   # Differential Nodes Table - nodes involved in significant differential edges
          #   h5("Nodes with Differential Connectivity"),
          #   p(style = "font-size: 12px; color: #666; margin-bottom: 10px;",
          #     "Nodes ranked by number of significant differential edges they are involved in (p < 0.05)"),
          #   DTOutput(ns("diffnet_nodes_table")),
          #
          #   br(),
          #
          #   # Differential Edges Table - significant differential edges
          #   h5("Significant Differential Edges"),
          #   p(style = "font-size: 12px; color: #666; margin-bottom: 10px;",
          #     "Edges with significantly different associations between groups (adjusted p < 0.05)"),
          #   DTOutput(ns("diffnet_edges_table"))
          # ),

          br(),

          # Node Properties Table - wrapped in card for visibility
          card(
            full_screen = TRUE,
            fill = FALSE,
            card_header(
              style = "background-color: #2FA4E7; color: white;",
              "Node Properties"
            ),
            card_body(
              DTOutput(ns("node_table"))
            )
          ),

          br(),

          # Module Assignment Table - wrapped in card for visibility
          card(
            full_screen = TRUE,
            fill = FALSE,
            card_header(
              style = "background-color: #2FA4E7; color: white;",
              "Module Assignment (Community Detection)"
            ),
            card_body(
              p(
                style = "font-size: 12px; color: #666;",
                "Method: Louvain algorithm (igraph::cluster_louvain)"
              ),
              div(
                style = "margin-bottom: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; font-size: 13px;",
                tags$p(
                  tags$strong("How modules are calculated: "),
                  "The Louvain algorithm detects communities by optimizing modularity - a measure of how densely connected nodes are within communities compared to between communities. ",
                  "Higher modularity indicates clearer community structure."
                ),
                tags$p(
                  tags$strong("Edge visualization: "),
                  tags$span(style = "color: forestgreen; font-weight: bold;", "Green edges = positive correlations"),
                  " | ",
                  tags$span(style = "color: crimson; font-weight: bold;", "Red edges = negative correlations"),
                  " | Edge width ∝ |correlation|"
                )
              ),
              DTOutput(ns("module_table"))
            )
          )
        )
      )
    )
  )
}

#' Network Visualization Server
network_visualization_Server <- function(id, network_results_reactive,
                                        group_column_reactive,
                                        analysis_completed_reactive,
                                        tax_rank_reactive = reactive("none")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Show/hide analysis message
    observe({
      if (analysis_completed_reactive()) {
        shinyjs::hide("analysis_message")
      } else {
        shinyjs::show("analysis_message")
      }
    })

    # Update group choices
    observe({
      req(network_results_reactive())

      results <- network_results_reactive()
      groups <- names(results)[names(results) != "comparison"]

      updateSelectInput(session, "selected_group",
                       choices = groups,
                       selected = groups[1])
    })

    # Get selected group's network
    selected_network <- reactive({
      req(input$selected_group, network_results_reactive())
      network_results_reactive()[[input$selected_group]]
    })

    # Get group names for comparison view
    comparison_groups <- reactive({
      req(network_results_reactive())
      results <- network_results_reactive()
      groups <- names(results)[names(results) != "comparison"]
      groups
    })

    # Group titles for comparison view
    output$group1_title <- renderText({
      groups <- comparison_groups()
      if (length(groups) >= 1) groups[1] else ""
    })

    output$group2_title <- renderText({
      groups <- comparison_groups()
      if (length(groups) >= 2) groups[2] else ""
    })

    # Network plot
    output$network_plot <- renderPlotly({
      req(selected_network())

      net <- selected_network()
      g <- net$igraph

      if (vcount(g) == 0) {
        return(plotly_empty() %>%
                layout(title = "No network to display (no edges met threshold)"))
      }

      # Filter edges by weight if needed
      if (input$min_edge_weight > 0) {
        edges_to_keep <- which(abs(E(g)$weight) >= input$min_edge_weight)
        if (length(edges_to_keep) == 0) {
          return(plotly_empty() %>%
                  layout(title = "No edges meet the minimum weight threshold"))
        }
        g <- subgraph.edges(g, edges_to_keep, delete.vertices = TRUE)
      }

      # Calculate layout
      layout_coords <- switch(
        input$layout_algorithm,
        "fr" = layout_with_fr(g),
        "kk" = layout_with_kk(g),
        "circle" = layout_in_circle(g),
        "grid" = layout_on_grid(g),
        "sphere" = layout_on_sphere(g),
        "random" = layout_randomly(g)
      )

      # Calculate node properties
      node_degree <- degree(g)
      node_betweenness <- betweenness(g)

      # Get IVI values (from network results, matched to current graph nodes)
      node_ivi <- rep(0, vcount(g))
      if (!is.null(net$ivi) && nrow(net$ivi) > 0) {
        ivi_match <- match(V(g)$name, net$ivi$node)
        node_ivi <- ifelse(is.na(ivi_match), 0, net$ivi$IVI[ivi_match])
      }

      # Detect communities for coloring
      communities <- cluster_fast_greedy(g)

      # Determine node colors and prepare discrete coloring info
      color_is_discrete <- FALSE
      color_categories <- NULL
      color_palette <- NULL

      if (input$node_color_by == "taxonomy") {
        # Get taxonomy from vertex attributes
        if (!is.null(V(g)$taxonomy)) {
          taxonomy_names <- V(g)$taxonomy
          # Convert to factor and then to numeric for coloring
          taxonomy_factor <- factor(taxonomy_names)
          node_colors <- as.numeric(taxonomy_factor)
          color_is_discrete <- TRUE
          color_categories <- levels(taxonomy_factor)
          # Create discrete color palette (using Set3 for many categories)
          n_colors <- length(color_categories)
          if (n_colors <= 12) {
            color_palette <- RColorBrewer::brewer.pal(max(3, n_colors), "Set3")[1:n_colors]
          } else {
            # For more categories, use rainbow
            color_palette <- rainbow(n_colors)
          }
        } else {
          # Fallback if no taxonomy data
          node_colors <- rep(1, vcount(g))
        }
      } else {
        # Continuous coloring for other metrics
        node_colors <- switch(
          input$node_color_by,
          "degree" = node_degree,
          "betweenness" = node_betweenness,
          "module" = membership(communities)
        )
      }

      # Determine node sizes
      node_sizes <- switch(
        input$node_size_by,
        "degree" = node_degree,
        "betweenness" = node_betweenness,
        "ivi" = node_ivi,
        "uniform" = rep(10, vcount(g))
      )

      # Normalize sizes
      if (input$node_size_by != "uniform") {
        node_sizes <- scales::rescale(node_sizes, to = c(5, 30))
      }

      # Create edge list (get.edgelist returns character vertex names)
      edge_list <- get.edgelist(g, names = FALSE)  # Get numeric indices

      # Use signed_weight for color (positive=green, negative=red)
      # Use absolute weight for width (correlation strength)
      edge_signed_weights <- if (!is.null(E(g)$signed_weight)) E(g)$signed_weight else E(g)$weight
      edge_abs_weights <- E(g)$weight  # Already absolute values

      # Prepare data for plotly with transparent edges
      edge_shapes <- list()
      for (i in 1:nrow(edge_list)) {
        v1 <- as.integer(edge_list[i, 1])
        v2 <- as.integer(edge_list[i, 2])

        # Edge color: Green for positive correlation, Red for negative correlation
        rgba_color <- ifelse(edge_signed_weights[i] > 0,
                            sprintf("rgba(34, 139, 34, %s)", input$edge_opacity),  # forest green (positive)
                            sprintf("rgba(220, 20, 60, %s)", input$edge_opacity))  # crimson red (negative)

        # Edge width: Scale based on correlation magnitude (absolute value)
        # Range from 0.5 (weak) to 6 (strong) pixels
        abs_weight <- edge_abs_weights[i]
        edge_width <- 0.5 + (abs_weight * 5.5)  # 0.5 to 6 pixels

        edge_shapes[[i]] <- list(
          type = "line",
          x0 = layout_coords[v1, 1],
          y0 = layout_coords[v1, 2],
          x1 = layout_coords[v2, 1],
          y1 = layout_coords[v2, 2],
          line = list(
            color = rgba_color,
            width = edge_width
          ),
          layer = "below"  # Draw edges below nodes
        )
      }

      # Node labels (only shown on plot when checkbox is checked)
      node_labels <- V(g)$name

      # Create hover labels with taxonomy information (shown on mouse hover)
      if (!is.null(V(g)$taxonomy)) {
        hover_labels <- paste0(
          "<b style='font-size:14px'>", V(g)$taxonomy, "</b><br>",
          "<span style='color:#666'>ID: ", V(g)$name, "</span><br>",
          "<br>",
          "Degree: ", node_degree, "<br>",
          "Betweenness: ", round(node_betweenness, 2)
        )
      } else {
        hover_labels <- paste0(
          "<b style='font-size:14px'>", V(g)$name, "</b><br>",
          "<br>",
          "Degree: ", node_degree, "<br>",
          "Betweenness: ", round(node_betweenness, 2)
        )
      }

      # Configure marker appearance based on color type
      marker_config <- list(
        size = node_sizes,
        line = list(width = 1, color = "white")
      )

      if (color_is_discrete && !is.null(color_palette)) {
        # Discrete coloring for taxonomy
        # Map numeric indices to colors
        marker_config$color <- sapply(node_colors, function(i) color_palette[i])
        marker_config$showscale = FALSE
      } else {
        # Continuous coloring for metrics
        marker_config$color <- node_colors
        marker_config$colorscale = "Viridis"
        marker_config$showscale = TRUE
        marker_config$colorbar = list(title = tools::toTitleCase(input$node_color_by))
      }

      node_trace <- list(
        x = layout_coords[, 1],
        y = layout_coords[, 2],
        mode = "markers",  # Always markers only, labels via hover
        type = "scatter",
        text = hover_labels,  # Use hover_labels for hover text
        marker = marker_config,
        hoverinfo = "text"
      )

      # Create plot
      plot_ly() %>%
        add_trace(
          x = node_trace$x,
          y = node_trace$y,
          mode = node_trace$mode,
          type = node_trace$type,
          text = node_trace$text,
          marker = node_trace$marker,
          hoverinfo = node_trace$hoverinfo
        ) %>%
        layout(
          title = paste("Network Graph -", input$selected_group,
                       ifelse(!is.null(tax_rank_reactive()) && tax_rank_reactive() != "none",
                              paste0(" (", tax_rank_reactive(), " level)"),
                              "")),
          showlegend = FALSE,
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          shapes = edge_shapes,
          hoverlabel = list(bgcolor = "white", font = list(size = 12))
        )
    })

    # Helper function to create network plot for a given group
    create_network_plot <- function(group_name) {
      req(network_results_reactive())

      results <- network_results_reactive()
      net <- results[[group_name]]

      if (is.null(net)) {
        return(plotly_empty() %>% layout(title = paste("No data for", group_name)))
      }

      g <- net$igraph

      if (vcount(g) == 0) {
        return(plotly_empty() %>%
                layout(title = paste("No network for", group_name)))
      }

      # Filter edges by weight if needed
      if (input$min_edge_weight > 0) {
        edges_to_keep <- which(abs(E(g)$weight) >= input$min_edge_weight)
        if (length(edges_to_keep) == 0) {
          return(plotly_empty() %>%
                  layout(title = paste(group_name, "- No edges meet threshold")))
        }
        g <- subgraph.edges(g, edges_to_keep, delete.vertices = TRUE)
      }

      # Calculate layout
      layout_coords <- switch(
        input$layout_algorithm,
        "fr" = layout_with_fr(g),
        "kk" = layout_with_kk(g),
        "circle" = layout_in_circle(g),
        "grid" = layout_on_grid(g),
        "sphere" = layout_on_sphere(g),
        "random" = layout_randomly(g)
      )

      # Calculate node properties
      node_degree <- degree(g)
      node_betweenness <- betweenness(g)
      communities <- cluster_fast_greedy(g)

      # Get IVI values (from network results, matched to current graph nodes)
      node_ivi <- rep(0, vcount(g))
      if (!is.null(net$ivi) && nrow(net$ivi) > 0) {
        ivi_match <- match(V(g)$name, net$ivi$node)
        node_ivi <- ifelse(is.na(ivi_match), 0, net$ivi$IVI[ivi_match])
      }

      # Determine node colors and prepare discrete coloring info
      color_is_discrete <- FALSE
      color_categories <- NULL
      color_palette <- NULL

      if (input$node_color_by == "taxonomy") {
        # Get taxonomy from vertex attributes
        if (!is.null(V(g)$taxonomy)) {
          taxonomy_names <- V(g)$taxonomy
          # Convert to factor and then to numeric for coloring
          taxonomy_factor <- factor(taxonomy_names)
          node_colors <- as.numeric(taxonomy_factor)
          color_is_discrete <- TRUE
          color_categories <- levels(taxonomy_factor)
          # Create discrete color palette (using Set3 for many categories)
          n_colors <- length(color_categories)
          if (n_colors <= 12) {
            color_palette <- RColorBrewer::brewer.pal(max(3, n_colors), "Set3")[1:n_colors]
          } else {
            # For more categories, use rainbow
            color_palette <- rainbow(n_colors)
          }
        } else {
          # Fallback if no taxonomy data
          node_colors <- rep(1, vcount(g))
        }
      } else {
        # Continuous coloring for other metrics
        node_colors <- switch(
          input$node_color_by,
          "degree" = node_degree,
          "betweenness" = node_betweenness,
          "module" = membership(communities)
        )
      }

      # Determine node sizes
      node_sizes <- switch(
        input$node_size_by,
        "degree" = node_degree,
        "betweenness" = node_betweenness,
        "ivi" = node_ivi,
        "uniform" = rep(10, vcount(g))
      )

      if (input$node_size_by != "uniform") {
        node_sizes <- scales::rescale(node_sizes, to = c(5, 30))
      }

      # Node labels (only shown on plot when checkbox is checked)
      node_labels <- V(g)$name

      # Create hover labels with taxonomy information (shown on mouse hover)
      if (!is.null(V(g)$taxonomy)) {
        hover_labels <- paste0(
          "<b style='font-size:14px'>", V(g)$taxonomy, "</b><br>",
          "<span style='color:#666'>ID: ", V(g)$name, "</span><br>",
          "<br>",
          "Degree: ", node_degree, "<br>",
          "Betweenness: ", round(node_betweenness, 2)
        )
      } else {
        hover_labels <- paste0(
          "<b style='font-size:14px'>", V(g)$name, "</b><br>",
          "<br>",
          "Degree: ", node_degree, "<br>",
          "Betweenness: ", round(node_betweenness, 2)
        )
      }

      # Configure marker appearance based on color type
      marker_config <- list(
        size = node_sizes,
        line = list(width = 1, color = "white")
      )

      if (color_is_discrete && !is.null(color_palette)) {
        # Discrete coloring for taxonomy
        marker_config$color <- sapply(node_colors, function(i) color_palette[i])
        marker_config$showscale = FALSE
      } else {
        # Continuous coloring for metrics
        marker_config$color <- node_colors
        marker_config$colorscale = "Viridis"
        marker_config$showscale = TRUE
        marker_config$colorbar = list(title = tools::toTitleCase(input$node_color_by))
      }

      # Create edge list
      edge_list <- get.edgelist(g, names = FALSE)

      # Use signed_weight for color (positive=green, negative=red)
      # Use absolute weight for width (correlation strength)
      edge_signed_weights <- if (!is.null(E(g)$signed_weight)) E(g)$signed_weight else E(g)$weight
      edge_abs_weights <- E(g)$weight  # Already absolute values

      # Prepare edge shapes with transparency
      edge_shapes <- list()
      for (i in 1:nrow(edge_list)) {
        v1 <- as.integer(edge_list[i, 1])
        v2 <- as.integer(edge_list[i, 2])

        # Edge color: Green for positive correlation, Red for negative correlation
        rgba_color <- ifelse(edge_signed_weights[i] > 0,
                            sprintf("rgba(34, 139, 34, %s)", input$edge_opacity),  # forest green (positive)
                            sprintf("rgba(220, 20, 60, %s)", input$edge_opacity))  # crimson red (negative)

        # Edge width: Scale based on correlation magnitude (absolute value)
        # Range from 0.5 (weak) to 6 (strong) pixels
        abs_weight <- edge_abs_weights[i]
        edge_width <- 0.5 + (abs_weight * 5.5)  # 0.5 to 6 pixels

        edge_shapes[[i]] <- list(
          type = "line",
          x0 = layout_coords[v1, 1],
          y0 = layout_coords[v1, 2],
          x1 = layout_coords[v2, 1],
          y1 = layout_coords[v2, 2],
          line = list(
            color = rgba_color,
            width = edge_width
          ),
          layer = "below"  # Draw edges below nodes
        )
      }

      # Create plot
      plot_ly() %>%
        add_trace(
          x = layout_coords[, 1],
          y = layout_coords[, 2],
          mode = "markers",  # Always markers only, labels via hover
          type = "scatter",
          text = hover_labels,  # Use hover_labels for hover text
          marker = marker_config,
          hoverinfo = "text"
        ) %>%
        layout(
          title = paste(group_name,
                       ifelse(!is.null(tax_rank_reactive()) && tax_rank_reactive() != "none",
                              paste0(" (", tax_rank_reactive(), " level)"),
                              "")),
          showlegend = FALSE,
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          shapes = edge_shapes,
          hoverlabel = list(bgcolor = "white", font = list(size = 12))
        )
    }

    # Comparison plot - Group 1
    output$network_plot_group1 <- renderPlotly({
      req(comparison_groups())
      groups <- comparison_groups()
      if (length(groups) >= 1) {
        create_network_plot(groups[1])
      } else {
        plotly_empty()
      }
    })

    # Comparison plot - Group 2
    output$network_plot_group2 <- renderPlotly({
      req(comparison_groups())
      groups <- comparison_groups()
      if (length(groups) >= 2) {
        create_network_plot(groups[2])
      } else {
        plotly_empty()
      }
    })

    # =========================================================================
    # [COMMENTED OUT] Differential Network Analysis (NetCoMi) - To be implemented later
    # =========================================================================

    # # Reactive value to store diffnet results
    # diffnet_results <- reactiveVal(NULL)
    #
    # # Run differential network analysis
    # observeEvent(input$run_diffnet, {
    # req(network_results_reactive(), comparison_groups())
    # 
    # groups <- comparison_groups()
    # if (length(groups) < 2) {
    # showNotification(
    # "Need at least 2 groups for differential analysis",
    # type = "error",
    # duration = 5
    # )
    # return()
    # }
    # 
    # showNotification(
    # "Running differential network analysis... This may take a minute",
    # type = "message",
    # duration = NULL,
    # id = "diffnet_notification"
    # )
    # 
    # tryCatch({
        # Check if NetCoMi is installed
    # if (!requireNamespace("NetCoMi", quietly = TRUE)) {
    # showNotification(
    # "NetCoMi package is required. Install with: install.packages('NetCoMi')",
    # type = "error",
    # duration = 10
    # )
    # return()
    # }
    # 
    # results <- network_results_reactive()
    # net1 <- results[[groups[1]]]$igraph
    # net2 <- results[[groups[2]]]$igraph
    # 
        # Get adjacency matrices from pre-constructed networks
    # adj1_raw <- as_adjacency_matrix(net1, type = "both", attr = "weight", sparse = FALSE)
    # adj2_raw <- as_adjacency_matrix(net2, type = "both", attr = "weight", sparse = FALSE)
    # 
        # Convert to regular matrices
    # adj1_raw <- as.matrix(adj1_raw)
    # adj2_raw <- as.matrix(adj2_raw)
    # 
        # Get node names from both networks
    # nodes1 <- rownames(adj1_raw)
    # nodes2 <- rownames(adj2_raw)
    # 
        # Create union of all nodes (both networks must have same nodes for diffnet)
    # all_nodes <- sort(union(nodes1, nodes2))
    # n_nodes <- length(all_nodes)
    # 
    # message("Aligning adjacency matrices: ", length(nodes1), " nodes in group1, ",
    # length(nodes2), " nodes in group2, ", n_nodes, " nodes in union")
    # 
        # Create aligned adjacency matrices with all nodes
    # adj1 <- matrix(0, nrow = n_nodes, ncol = n_nodes,
    # dimnames = list(all_nodes, all_nodes))
    # adj2 <- matrix(0, nrow = n_nodes, ncol = n_nodes,
    # dimnames = list(all_nodes, all_nodes))
    # 
        # Fill in values from original matrices
    # adj1[nodes1, nodes1] <- adj1_raw
    # adj2[nodes2, nodes2] <- adj2_raw
    # 
    # message("Aligned matrix dimensions: ", nrow(adj1), " x ", ncol(adj1))
    # 
        # Construct NetCoMi network objects using pre-computed association matrices
        # This bypasses network construction and uses our existing networks
    # net_obj <- NetCoMi::netConstruct(
    # data = adj1,
    # data2 = adj2,
    # dataType = "condDependence",  # Treat as association matrices
    # sparsMethod = "none"  # Networks already filtered/thresholded
    # )
    # 
        # Get sample sizes for each group (required when count matrices are missing)
    # n1 <- results[[groups[1]]]$n_samples
    # n2 <- results[[groups[2]]]$n_samples
    # 
        # Validate sample sizes
    # if (is.null(n1) || is.null(n2)) {
    # stop("Sample size information not found in network results. Please re-run network analysis.")
    # }
    # 
    # message("Differential network analysis - n1: ", n1, ", n2: ", n2)
    # 
        # Run differential network analysis
    # diffnet_result <- NetCoMi::diffnet(
    # net_obj,
    # diffMethod = input$diffnet_method,
    # adjust = input$diffnet_adjust,
    # n1 = n1,
    # n2 = n2
    # )
    # 
        # Debug: Print diffnet result structure
    # message("=== DIFFNET RESULT STRUCTURE ===")
    # message("Class: ", class(diffnet_result))
    # message("Names: ", paste(names(diffnet_result), collapse = ", "))
    # if (!is.null(diffnet_result$diffMat)) {
    # message("diffMat dim: ", paste(dim(diffnet_result$diffMat), collapse = " x "))
    # }
    # if (!is.null(diffnet_result$pvalDiff)) {
    # message("pvalDiff available: TRUE")
    # }
    # message("================================")
    # 
        # Store results
    # diffnet_results(list(
    # result = diffnet_result,
    # group1 = groups[1],
    # group2 = groups[2],
    # method = input$diffnet_method,
    # adjust = input$diffnet_adjust
    # ))
    # 
    # removeNotification("diffnet_notification")
    # showNotification(
    # "Differential network analysis completed!",
    # type = "message",
    # duration = 5
    # )
    # 
    # }, error = function(e) {
    # removeNotification("diffnet_notification")
    # showNotification(
    # paste("Error in differential analysis:", e$message),
    # type = "error",
    # duration = 10
    # )
    # })
    # })
    # 
    # Hide warning when differential analysis is run
    # observe({
    # if (!is.null(diffnet_results())) {
    # shinyjs::hide("diffnet_warning")
    # } else {
    # shinyjs::show("diffnet_warning")
    # }
    # })
    # 
    # Display method and adjustment
    # output$diffnet_method_display <- renderText({
    # if (is.null(diffnet_results())) return("-")
    # diffnet_results()$method
    # })
    # 
    # output$diffnet_adjust_display <- renderText({
    # if (is.null(diffnet_results())) return("-")
    # diffnet_results()$adjust
    # })
    # 
    # Differential nodes count - count nodes involved in significant differential edges
    # Uses canonical approach: diffAdjustMat != 0 (NetCoMi standard)
    # output$diffnet_n_diff_nodes <- renderText({
    # if (is.null(diffnet_results())) return("-")
    # res <- diffnet_results()$result
    # 
    # tryCatch({
    # diffAdjustMat <- res$diffAdjustMat
    # pAdjustMat <- res$pAdjustMat
    # 
        # Determine significant edges using canonical approach
    # if (!is.null(diffAdjustMat) && sum(diffAdjustMat != 0) > 0) {
    # sig_edges <- diffAdjustMat != 0
    # } else if (!is.null(pAdjustMat)) {
    # sig_edges <- pAdjustMat < 0.05
    # } else {
    # return("N/A")
    # }
    # diag(sig_edges) <- FALSE
    # 
    # n_nodes <- nrow(sig_edges)
        # Count nodes with at least 1 significant edge
    # nodes_with_sig_edges <- sapply(1:n_nodes, function(i) {
    # any(sig_edges[i, ], na.rm = TRUE) || any(sig_edges[, i], na.rm = TRUE)
    # })
    # as.character(sum(nodes_with_sig_edges))
    # }, error = function(e) "Error")
    # })
    # 
    # Differential edges count - Uses canonical approach: diffAdjustMat != 0 (NetCoMi standard)
    # output$diffnet_n_diff_edges <- renderText({
    # if (is.null(diffnet_results())) return("-")
    # res <- diffnet_results()$result
    # 
    # tryCatch({
    # diffAdjustMat <- res$diffAdjustMat
    # pAdjustMat <- res$pAdjustMat
    # 
        # Determine significant edges using canonical approach
    # if (!is.null(diffAdjustMat) && sum(diffAdjustMat != 0) > 0) {
    # sig_edges <- diffAdjustMat != 0
    # } else if (!is.null(pAdjustMat)) {
    # sig_edges <- pAdjustMat < 0.05
    # } else {
    # return("N/A")
    # }
    # 
        # Count significant edges (upper triangle only to avoid double counting)
    # n_sig <- sum(sig_edges[upper.tri(sig_edges)], na.rm = TRUE)
    # as.character(n_sig)
    # }, error = function(e) "Error")
    # })
    # 
    # Differential nodes table - use diffAdjustMat (NetCoMi's canonical approach)
    # output$diffnet_nodes_table <- renderDT({
    # if (is.null(diffnet_results())) {
    # return(datatable(
    # data.frame(Message = "Run differential analysis to see results"),
    # options = list(dom = 't'),
    # rownames = FALSE
    # ))
    # }
    # res <- diffnet_results()$result
    # 
    # tryCatch({
    # diffAdjustMat <- res$diffAdjustMat
    # pAdjustMat <- res$pAdjustMat
    # assoMat1 <- res$assoMat1
    # assoMat2 <- res$assoMat2
    # 
    # if (is.null(assoMat1) || is.null(assoMat2)) {
    # return(datatable(
    # data.frame(Message = "Association matrices not available"),
    # options = list(dom = 't'),
    # rownames = FALSE
    # ))
    # }
    # 
    # node_names <- rownames(assoMat1)
    # n_nodes <- length(node_names)
    # 
        # Determine significant edges: diffAdjustMat != 0 (canonical) or pAdjustMat < 0.05
    # if (!is.null(diffAdjustMat) && sum(diffAdjustMat != 0) > 0) {
    # sig_edges <- diffAdjustMat != 0
    # method_used <- "diffAdjustMat != 0"
    # } else if (!is.null(pAdjustMat)) {
    # sig_edges <- pAdjustMat < 0.05
    # method_used <- "pAdjust < 0.05"
    # } else {
    # return(datatable(data.frame(Message = "No p-value data available"), options = list(dom = 't'), rownames = FALSE))
    # }
    # diag(sig_edges) <- FALSE
    # 
        # Count significant edges per node
    # sig_edges_count <- sapply(1:n_nodes, function(i) {
    # sum(sig_edges[i, ], na.rm = TRUE)  # count per row (already symmetric)
    # })
    # 
        # Only include nodes with significant edges
    # nodes_with_sig <- which(sig_edges_count > 0)
    # 
    # if (length(nodes_with_sig) == 0) {
    # return(datatable(
    # data.frame(Message = paste0("No nodes with significant differential edges (", method_used, ")")),
    # options = list(dom = 't'),
    # rownames = FALSE
    # ))
    # }
    # 
        # Build node table for nodes with significant edges
    # node_df <- data.frame(
    # Node = node_names[nodes_with_sig],
    # Sig_Edges = sig_edges_count[nodes_with_sig],
    # stringsAsFactors = FALSE
    # )
    # 
        # Sort by significant edges count
    # node_df <- node_df[order(-node_df$Sig_Edges), ]
    # 
    # datatable(
    # node_df,
    # options = list(pageLength = 20, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv')),
    # extensions = 'Buttons',
    # rownames = FALSE,
    # caption = paste0(nrow(node_df), " nodes involved in significant differential edges (", method_used, ")")
    # )
    # 
    # }, error = function(e) {
    # datatable(data.frame(Message = paste("Error:", e$message)), options = list(dom = 't'), rownames = FALSE)
    # })
    # })
    # 
    # Differential edges table - use diffAdjustMat (NetCoMi's canonical approach)
    # output$diffnet_edges_table <- renderDT({
    # if (is.null(diffnet_results())) {
    # return(datatable(
    # data.frame(Message = "Run differential analysis to see results"),
    # options = list(dom = 't'),
    # rownames = FALSE
    # ))
    # }
    # res <- diffnet_results()$result
    # groups_info <- diffnet_results()
    # 
    # tryCatch({
        # diffAdjustMat: non-zero only for significant differential edges
    # diffAdjustMat <- res$diffAdjustMat
    # pAdjustMat <- res$pAdjustMat
    # assoMat1 <- res$assoMat1
    # assoMat2 <- res$assoMat2
    # 
    # if (is.null(assoMat1) || is.null(assoMat2)) {
    # return(datatable(
    # data.frame(Message = "Association matrices not available"),
    # options = list(dom = 't'),
    # rownames = FALSE
    # ))
    # }
    # 
    # node_names <- rownames(assoMat1)
    # g1_name <- if (!is.null(groups_info$group1)) groups_info$group1 else "Group1"
    # g2_name <- if (!is.null(groups_info$group2)) groups_info$group2 else "Group2"
    # 
        # Find significant edges using diffAdjustMat != 0 (NetCoMi standard)
    # if (!is.null(diffAdjustMat) && sum(diffAdjustMat != 0) > 0) {
          # Use diffAdjustMat (canonical method)
    # idx <- which(diffAdjustMat != 0 & upper.tri(diffAdjustMat), arr.ind = TRUE)
    # 
    # edge_df <- data.frame(
    # Taxa1 = node_names[idx[, 1]],
    # Taxa2 = node_names[idx[, 2]],
    # stringsAsFactors = FALSE
    # )
    # edge_df[[paste0("Asso_", g1_name)]] <- assoMat1[idx]
    # edge_df[[paste0("Asso_", g2_name)]] <- assoMat2[idx]
    # edge_df$Difference <- diffAdjustMat[idx]
    # if (!is.null(pAdjustMat)) edge_df$pAdjust <- pAdjustMat[idx]
    # 
    # edge_df <- edge_df[order(edge_df$pAdjust), ]
    # 
    # datatable(
    # edge_df,
    # options = list(pageLength = 20, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv')),
    # extensions = 'Buttons',
    # rownames = FALSE,
    # caption = paste0(nrow(edge_df), " significant differential edges (diffAdjustMat != 0)")
    # ) %>%
    # formatRound(columns = which(sapply(edge_df, is.numeric)), digits = 4)
    # 
    # } else if (!is.null(pAdjustMat)) {
          # Fallback: use pAdjustMat < 0.05
    # idx <- which(pAdjustMat < 0.05 & upper.tri(pAdjustMat), arr.ind = TRUE)
    # 
    # if (nrow(idx) > 0) {
    # edge_df <- data.frame(
    # Taxa1 = node_names[idx[, 1]],
    # Taxa2 = node_names[idx[, 2]],
    # stringsAsFactors = FALSE
    # )
    # edge_df[[paste0("Asso_", g1_name)]] <- assoMat1[idx]
    # edge_df[[paste0("Asso_", g2_name)]] <- assoMat2[idx]
    # edge_df$Difference <- assoMat1[idx] - assoMat2[idx]
    # edge_df$pAdjust <- pAdjustMat[idx]
    # edge_df <- edge_df[order(edge_df$pAdjust), ]
    # 
    # datatable(
    # edge_df,
    # options = list(pageLength = 20, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv')),
    # extensions = 'Buttons',
    # rownames = FALSE,
    # caption = paste0(nrow(edge_df), " significant edges (pAdjust < 0.05)")
    # ) %>%
    # formatRound(columns = which(sapply(edge_df, is.numeric)), digits = 4)
    # } else {
            # No significant - show top 20 by p-value
    # idx <- which(upper.tri(pAdjustMat), arr.ind = TRUE)
    # edge_df <- data.frame(
    # Taxa1 = node_names[idx[, 1]],
    # Taxa2 = node_names[idx[, 2]],
    # stringsAsFactors = FALSE
    # )
    # edge_df[[paste0("Asso_", g1_name)]] <- assoMat1[idx]
    # edge_df[[paste0("Asso_", g2_name)]] <- assoMat2[idx]
    # edge_df$Difference <- assoMat1[idx] - assoMat2[idx]
    # edge_df$pAdjust <- pAdjustMat[idx]
    # edge_df <- edge_df[order(edge_df$pAdjust), ]
    # edge_df <- head(edge_df, 20)
    # 
    # datatable(
    # edge_df,
    # options = list(pageLength = 20, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv')),
    # extensions = 'Buttons',
    # rownames = FALSE,
    # caption = "No significant edges. Showing top 20 by p-value."
    # ) %>%
    # formatRound(columns = which(sapply(edge_df, is.numeric)), digits = 4)
    # }
    # } else {
    # datatable(data.frame(Message = "P-value matrix not available"), options = list(dom = 't'), rownames = FALSE)
    # }
    # 
    # }, error = function(e) {
    # datatable(data.frame(Message = paste("Error:", e$message)), options = list(dom = 't'), rownames = FALSE)
    # })
    # })
    # 
    # Differential network plot - show significant differences
    # Uses canonical approach: diffAdjustMat != 0 (NetCoMi standard)
    # output$diffnet_plot <- renderPlotly({
    # if (is.null(diffnet_results())) {
    # return(plotly_empty() %>%
    # layout(title = list(
    # text = "Run differential analysis to see results",
    # font = list(size = 14, color = "#666")
    # )))
    # }
    # res <- diffnet_results()$result
    # 
    # tryCatch({
    # diffAdjustMat <- res$diffAdjustMat
    # pAdjustMat <- res$pAdjustMat
    # pvalsMat <- res$pvalsMat  # Raw p-values for threshold breakdown
    # 
        # Check required matrices
    # if (is.null(pAdjustMat) && is.null(diffAdjustMat)) {
    # return(plotly_empty() %>%
    # layout(title = "Required matrices not available"))
    # }
    # 
        # Use diffAdjustMat as the canonical source
    # if (!is.null(diffAdjustMat)) {
    # n_nodes <- nrow(diffAdjustMat)
    # 
          # Canonical: diffAdjustMat != 0 indicates significant differential edges
    # sig_edges_canonical <- diffAdjustMat != 0
    # diag(sig_edges_canonical) <- FALSE
    # 
    # n_edges_canonical <- sum(sig_edges_canonical[upper.tri(sig_edges_canonical)], na.rm = TRUE)
    # n_nodes_canonical <- sum(sapply(1:n_nodes, function(i) {
    # any(sig_edges_canonical[i, ], na.rm = TRUE) || any(sig_edges_canonical[, i], na.rm = TRUE)
    # }))
    # 
          # Also show p-value breakdown if available
    # if (!is.null(pAdjustMat)) {
    # count_nodes_with_sig <- function(threshold) {
    # sig_mat <- pAdjustMat < threshold
    # diag(sig_mat) <- FALSE
    # sum(sapply(1:n_nodes, function(i) {
    # any(sig_mat[i, ], na.rm = TRUE) || any(sig_mat[, i], na.rm = TRUE)
    # }))
    # }
    # 
    # n_nodes_001 <- count_nodes_with_sig(0.01)
    # n_nodes_005 <- count_nodes_with_sig(0.05)
    # n_nodes_01 <- count_nodes_with_sig(0.1)
    # 
    # n_edges_001 <- sum(pAdjustMat[upper.tri(pAdjustMat)] < 0.01, na.rm = TRUE)
    # n_edges_005 <- sum(pAdjustMat[upper.tri(pAdjustMat)] < 0.05, na.rm = TRUE)
    # n_edges_01 <- sum(pAdjustMat[upper.tri(pAdjustMat)] < 0.1, na.rm = TRUE)
    # 
            # Create grouped bar chart with canonical + threshold breakdown
    # plot_df <- data.frame(
    # Category = rep(c("Nodes with\nDiff. Edges", "Differential\nEdges"), each = 4),
    # Threshold = rep(c("Canonical (adj.)", "pAdj < 0.01", "pAdj < 0.05", "pAdj < 0.1"), 2),
    # Count = c(n_nodes_canonical, n_nodes_001, n_nodes_005, n_nodes_01,
    # n_edges_canonical, n_edges_001, n_edges_005, n_edges_01)
    # )
    # 
            # Order threshold factor for better display
    # plot_df$Threshold <- factor(plot_df$Threshold,
    # levels = c("Canonical (adj.)", "pAdj < 0.01", "pAdj < 0.05", "pAdj < 0.1"))
    # } else {
            # Only canonical results available
    # plot_df <- data.frame(
    # Category = c("Nodes with\nDiff. Edges", "Differential\nEdges"),
    # Threshold = rep("Canonical (adj.)", 2),
    # Count = c(n_nodes_canonical, n_edges_canonical)
    # )
    # }
    # } else if (!is.null(pAdjustMat)) {
          # Fallback: only pAdjustMat available
    # n_nodes <- nrow(pAdjustMat)
    # 
    # count_nodes_with_sig <- function(threshold) {
    # sig_mat <- pAdjustMat < threshold
    # diag(sig_mat) <- FALSE
    # sum(sapply(1:n_nodes, function(i) {
    # any(sig_mat[i, ], na.rm = TRUE) || any(sig_mat[, i], na.rm = TRUE)
    # }))
    # }
    # 
    # n_nodes_001 <- count_nodes_with_sig(0.01)
    # n_nodes_005 <- count_nodes_with_sig(0.05)
    # n_nodes_01 <- count_nodes_with_sig(0.1)
    # 
    # n_edges_001 <- sum(pAdjustMat[upper.tri(pAdjustMat)] < 0.01, na.rm = TRUE)
    # n_edges_005 <- sum(pAdjustMat[upper.tri(pAdjustMat)] < 0.05, na.rm = TRUE)
    # n_edges_01 <- sum(pAdjustMat[upper.tri(pAdjustMat)] < 0.1, na.rm = TRUE)
    # 
    # plot_df <- data.frame(
    # Category = rep(c("Nodes with\nDiff. Edges", "Differential\nEdges"), each = 3),
    # Threshold = rep(c("pAdj < 0.01", "pAdj < 0.05", "pAdj < 0.1"), 2),
    # Count = c(n_nodes_001, n_nodes_005, n_nodes_01,
    # n_edges_001, n_edges_005, n_edges_01)
    # )
    # }
    # 
    # plot_ly(plot_df, x = ~Category, y = ~Count, color = ~Threshold,
    # type = "bar", text = ~Count, textposition = "outside") %>%
    # layout(
    # title = list(
    # text = "Differential Network Components (NetCoMi)",
    # font = list(size = 14)
    # ),
    # xaxis = list(title = ""),
    # yaxis = list(title = "Count"),
    # barmode = "group",
    # showlegend = TRUE,
    # legend = list(orientation = "h", y = -0.2)
    # )
    # }, error = function(e) {
    # plotly_empty() %>%
    # layout(title = paste("Error:", e$message))
    # })
    # })

    # Node properties table - shows Group column
    output$node_table <- renderDT({
      req(network_results_reactive())

      results <- network_results_reactive()
      groups <- names(results)[names(results) != "comparison"]

      # Build combined node data frame based on view mode
      all_node_data <- data.frame()

      if (input$view_mode == "compare") {
        # Compare mode: show nodes from all groups
        for (grp in groups) {
          g <- results[[grp]]$igraph
          if (!is.null(g) && vcount(g) > 0) {
            node_df <- data.frame(
              Group = grp,
              Node = V(g)$name,
              Degree = degree(g),
              Betweenness = betweenness(g),
              Closeness = closeness(g),
              Eigenvector = eigen_centrality(g)$vector,
              stringsAsFactors = FALSE
            )
            all_node_data <- rbind(all_node_data, node_df)
          }
        }
      } else {
        # Single mode: show nodes from selected group only
        req(input$selected_group)
        grp <- input$selected_group
        g <- results[[grp]]$igraph
        if (!is.null(g) && vcount(g) > 0) {
          all_node_data <- data.frame(
            Group = grp,
            Node = V(g)$name,
            Degree = degree(g),
            Betweenness = betweenness(g),
            Closeness = closeness(g),
            Eigenvector = eigen_centrality(g)$vector,
            stringsAsFactors = FALSE
          )
        }
      }

      if (nrow(all_node_data) == 0) {
        return(datatable(data.frame(Message = "No nodes in network")))
      }

      datatable(
        all_node_data,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          order = list(list(2, 'desc'))  # Sort by degree
        ),
        rownames = FALSE,
        filter = 'top'
      ) %>%
        formatRound(columns = c("Betweenness", "Closeness", "Eigenvector"), digits = 3)
    })

    # Module assignment table - shows Group column
    output$module_table <- renderDT({
      req(network_results_reactive())

      results <- network_results_reactive()
      groups <- names(results)[names(results) != "comparison"]

      tryCatch({
        all_module_data <- data.frame()
        summary_info <- list()

        # Determine which groups to process
        groups_to_process <- if (input$view_mode == "compare") {
          groups
        } else {
          req(input$selected_group)
          input$selected_group
        }

        for (grp in groups_to_process) {
          g <- results[[grp]]$igraph
          if (is.null(g) || vcount(g) == 0) next

          # Calculate community detection using Louvain algorithm
          communities <- cluster_louvain(g, weights = abs(E(g)$weight))
          membership_vec <- as.numeric(membership(communities))

          # Create module assignment data frame
          module_df <- data.frame(
            Group = grp,
            Node = V(g)$name,
            Module = membership_vec,
            Degree = degree(g),
            stringsAsFactors = FALSE
          )

          # Add taxonomy if available
          if (!is.null(V(g)$taxonomy)) {
            module_df$Taxonomy <- V(g)$taxonomy
          }

          # Add module size information
          module_sizes <- table(membership_vec)
          module_df$ModuleSize <- as.numeric(module_sizes[as.character(module_df$Module)])

          # Store summary info
          summary_info[[grp]] <- list(
            n_modules = length(unique(membership_vec)),
            modularity = modularity(communities)
          )

          all_module_data <- rbind(all_module_data, module_df)
        }

        if (nrow(all_module_data) == 0) {
          return(datatable(data.frame(Message = "No nodes in network")))
        }

        # Sort by Group, then Module, then by degree within module
        all_module_data <- all_module_data[order(all_module_data$Group, all_module_data$Module, -all_module_data$Degree), ]

        # Create summary caption
        summary_text <- paste(sapply(names(summary_info), function(grp) {
          paste0(grp, ": ", summary_info[[grp]]$n_modules, " modules (Q=", round(summary_info[[grp]]$modularity, 3), ")")
        }), collapse = " | ")

        datatable(
          all_module_data,
          options = list(
            pageLength = 20,
            scrollX = TRUE,
            order = list(list(0, 'asc'), list(2, 'asc')),  # Sort by Group, then Module
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
          ),
          extensions = 'Buttons',
          rownames = FALSE,
          filter = 'top',
          caption = htmltools::tags$caption(
            style = 'caption-side: top; text-align: left; font-size: 14px; margin-bottom: 10px;',
            htmltools::HTML(paste0("<strong>", summary_text, "</strong>"))
          )
        )
      }, error = function(e) {
        datatable(data.frame(Message = paste("Error calculating modules:", e$message)))
      })
    })

    # Download network plot
    output$download_network <- downloadHandler(
      filename = function() {
        paste0("network_", input$selected_group, "_", Sys.Date(), ".html")
      },
      content = function(file) {
        # Save the plotly plot as HTML
        p <- plotlyOutput("network_plot")
        htmlwidgets::saveWidget(as_widget(p), file)
      }
    )

    # Download GraphML
    output$download_graphml <- downloadHandler(
      filename = function() {
        paste0("network_", input$selected_group, "_", Sys.Date(), ".graphml")
      },
      content = function(file) {
        req(selected_network())
        g <- selected_network()$igraph
        write_graph(g, file, format = "graphml")
      }
    )

  })
}
