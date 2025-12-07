# =============================================================================
# Network Analysis Module
# Displays network analysis results and metrics
# =============================================================================

library(shiny)
library(shinyjs)
library(bslib)
library(plotly)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)

#' Network Analysis UI
network_analysis_UI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),

    card(
      full_screen = TRUE,
      class = "main-card",
      card_header(
        style = "background-color: #2FA4E7; color: white;",
        "Network Analysis Results"
      ),
      card_body(
        # Message before analysis
        div(
          id = ns("analysis_message"),
          style = "color: red; font-weight: bold; margin-bottom: 20px;",
          icon("exclamation-triangle"),
          " Run network analysis using the Master Controller"
        ),

        # Network metrics summary
        h4("Network Topology Metrics by Group"),
        DTOutput(ns("metrics_summary_table")),

        br(),

        # Single comprehensive metrics comparison plot
        h4(
          "Network Metrics Comparison",
          create_plot_download_btn(ns, "metrics", "Download"),
          style = "margin-top: 30px;"
        ),
        p("Comparison of all network topology metrics across groups", style = "margin-bottom: 20px;"),

        div(
          style = "margin-bottom: 30px; padding: 20px; background-color: #ffffff; border: 1px solid #dee2e6; border-radius: 5px;",
          plotlyOutput(ns("all_metrics_plot"), height = "600px")
        ),

        br(),

        # Edge Type Distribution by Taxonomy
        h4(
          "Edge Type Distribution by Taxonomy",
          create_plot_download_btn(ns, "edge_dist", "Download"),
          style = "margin-top: 30px;"
        ),
        p("Comparison of positive vs negative edge distribution across groups", style = "margin-bottom: 20px;"),

        # Taxonomic rank selector and Top N taxa filter
        div(
          style = "margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px;",
          fluidRow(
            column(6,
              selectInput(
                ns("edge_tax_rank"),
                "Taxonomic Rank:",
                choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                selected = "Species"
              )
            ),
            column(6,
              numericInput(
                ns("edge_top_n"),
                "Top N Taxa (by total edges):",
                value = 10,
                min = 5,
                max = 100,
                step = 5
              )
            )
          ),
          # Display total edge count info
          div(
            style = "margin-top: 10px; padding: 10px; background-color: #e9ecef; border-radius: 3px; font-size: 0.9em;",
            uiOutput(ns("edge_count_info"))
          )
        ),

        # Full-width mirrored edge polarity plot (no card wrapper)
        div(
          style = "margin-bottom: 30px; padding: 20px; background-color: #ffffff; border: 1px solid #dee2e6; border-radius: 5px;",
          plotlyOutput(ns("edge_polarity_plot"), height = "700px")
        )
      )
    )
  )
}

#' Network Analysis Server
network_analysis_Server <- function(id, filtered_ps_reactive, network_results_reactive,
                                   analysis_completed_reactive, tax_rank_reactive = reactive("Species")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Taxonomic rank hierarchy (lower index = higher level)
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    # Show/hide analysis message
    observe({
      if (analysis_completed_reactive()) {
        shinyjs::hide("analysis_message")
      } else {
        shinyjs::show("analysis_message")
      }
    })

    # Update taxonomic rank selector based on Master Controller selection
    observe({
      req(tax_rank_reactive())
      master_rank <- tax_rank_reactive()

      # Skip if master_rank is "none" (legacy)
      if (master_rank == "none") return()

      # Get allowed ranks (master rank and above only)
      idx <- which(tax_ranks == master_rank)
      if (length(idx) == 0) idx <- length(tax_ranks)
      allowed_ranks <- tax_ranks[1:idx]

      # Update edge_tax_rank selectInput
      updateSelectInput(session, "edge_tax_rank",
                       choices = allowed_ranks,
                       selected = master_rank)
    })

    # Metrics summary table
    output$metrics_summary_table <- renderDT({
      # Don't use req() - show message instead
      if (!analysis_completed_reactive()) {
        return(datatable(
          data.frame(Message = "Run network analysis using the Master Controller"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      results <- network_results_reactive()
      if (is.null(results) || is.null(results$comparison)) {
        return(datatable(
          data.frame(Message = "No network data available - please run analysis"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      metrics_df <- results$comparison

      # DEBUG: Log what we received
      message("=== METRICS SUMMARY TABLE DEBUG ===")
      message("metrics_df class: ", class(metrics_df))
      message("metrics_df dimensions: ", nrow(metrics_df), " rows x ", ncol(metrics_df), " cols")
      message("metrics_df columns: ", paste(colnames(metrics_df), collapse = ", "))
      if (nrow(metrics_df) > 0) {
        message("First row: ", paste(metrics_df[1, ], collapse = ", "))
      }
      message("===================================")

      # Check if data frame is empty or all edges are 0
      if (nrow(metrics_df) == 0) {
        # Return empty table with message
        return(datatable(
          data.frame(Message = "No network data available"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      # Add warnings for edge cases
      if (all(metrics_df$n_edges == 0, na.rm = TRUE)) {
        # Show notification to user
        showNotification(
          "Warning: No edges detected in the network. Try lowering the correlation threshold (r.threshold) or p-value threshold.",
          type = "warning",
          duration = 10
        )
        metrics_df$Warning <- "⚠ No edges detected"
      }

      # Warn if very few edges
      if (any(metrics_df$n_edges > 0 & metrics_df$n_edges < 10, na.rm = TRUE)) {
        showNotification(
          "Notice: Some networks have very few edges (<10). Consider adjusting filtering or threshold parameters.",
          type = "message",
          duration = 8
        )
      }

      # Find numeric columns
      numeric_cols <- which(sapply(metrics_df, is.numeric))

      dt <- datatable(
        metrics_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          autoWidth = TRUE,
          dom = 't'  # Only table, no search/pagination for small data
        ),
        rownames = FALSE
      )

      # Only format if there are numeric columns
      if (length(numeric_cols) > 0) {
        dt <- dt %>% formatRound(columns = numeric_cols, digits = 3)
      }

      dt
    })

    # Single comprehensive plot with all metrics
    output$all_metrics_plot <- renderPlotly({
      # Don't use req() - show empty plot instead
      if (!analysis_completed_reactive()) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Run network analysis first<br><sub>Use Master Controller to start analysis</sub>",
                  font = list(size = 14)
                )))
      }

      results <- network_results_reactive()
      if (is.null(results) || is.null(results$comparison)) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "No network data available<br><sub>Please run network analysis</sub>",
                  font = list(size = 14)
                )))
      }

      metrics_df <- results$comparison

      # DEBUG: Log what we received
      message("=== NODES_EDGES_PLOT DEBUG ===")
      message("metrics_df class: ", class(metrics_df))
      message("metrics_df dimensions: ", nrow(metrics_df), " rows x ", ncol(metrics_df), " cols")
      message("metrics_df columns: ", paste(colnames(metrics_df), collapse = ", "))
      message("==============================")

      # Check if required columns exist
      if (!all(c("group", "n_nodes", "n_edges") %in% colnames(metrics_df))) {
        message("ERROR: Missing required columns!")
        message("Required: group, n_nodes, n_edges")
        message("Available: ", paste(colnames(metrics_df), collapse = ", "))
        return(plotly_empty() %>%
                layout(title = "Error: Required columns missing from network metrics"))
      }

      # Check if data is empty
      if (nrow(metrics_df) == 0) {
        message("ERROR: metrics_df has 0 rows!")
        return(plotly_empty() %>%
                layout(title = list(
                  text = "No data available<br><sub>Please run network analysis first</sub>",
                  font = list(size = 14)
                )))
      }

      message("Proceeding with faceted subplot creation...")

      # Create faceted subplots with individual Y-axis scales
      tryCatch({
        # Get unique groups
        groups <- unique(metrics_df$group)
        n_groups <- length(groups)

        # Dynamically generate colors for any number of groups
        if (n_groups <= 8) {
          base_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                           "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")
        } else if (n_groups <= 12) {
          base_palette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                           "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                           "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
        } else {
          # For many groups, use HSL color wheel
          base_palette <- sapply(seq(0, 360 - 360/n_groups, length.out = n_groups), function(h) {
            rgb_val <- col2rgb(hsv(h/360, 0.7, 0.8))
            rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255)
          })
        }

        # Create named color vector for groups
        group_colors <- setNames(
          base_palette[((seq_along(groups) - 1) %% length(base_palette)) + 1],
          groups
        )

        # Create individual plots for each metric
        plots <- list()

        # 1. Nodes
        p1 <- plot_ly()
        for (grp in groups) {
          grp_data <- metrics_df %>% filter(group == grp)
          p1 <- p1 %>%
            add_trace(
              x = list(grp),
              y = grp_data$n_nodes,
              type = "bar",
              name = grp,
              marker = list(color = group_colors[[grp]]),
              legendgroup = grp,
              showlegend = TRUE,
              hovertemplate = paste0("<b>", grp, "</b><br>Nodes: %{y}<extra></extra>")
            )
        }
        p1 <- p1 %>% layout(
          xaxis = list(title = "Group"),
          yaxis = list(title = "Nodes"),
          title = list(text = "Nodes", font = list(size = 12))
        )

        # 2. Edges
        p2 <- plot_ly()
        for (grp in groups) {
          grp_data <- metrics_df %>% filter(group == grp)
          p2 <- p2 %>%
            add_trace(
              x = list(grp),
              y = grp_data$n_edges,
              type = "bar",
              name = grp,
              marker = list(color = group_colors[[grp]]),
              legendgroup = grp,
              showlegend = FALSE,
              hovertemplate = paste0("<b>", grp, "</b><br>Edges: %{y}<extra></extra>")
            )
        }
        p2 <- p2 %>% layout(
          xaxis = list(title = "Group"),
          yaxis = list(title = "Edges"),
          title = list(text = "Edges", font = list(size = 12))
        )

        # 3. Density
        p3 <- plot_ly()
        for (grp in groups) {
          grp_data <- metrics_df %>% filter(group == grp)
          p3 <- p3 %>%
            add_trace(
              x = list(grp),
              y = grp_data$density,
              type = "bar",
              name = grp,
              marker = list(color = group_colors[[grp]]),
              legendgroup = grp,
              showlegend = FALSE,
              hovertemplate = paste0("<b>", grp, "</b><br>Density: %{y:.4f}<extra></extra>")
            )
        }
        p3 <- p3 %>% layout(
          xaxis = list(title = "Group"),
          yaxis = list(title = "Density"),
          title = list(text = "Density", font = list(size = 12))
        )

        # 4. Average Degree
        p4 <- plot_ly()
        for (grp in groups) {
          grp_data <- metrics_df %>% filter(group == grp)
          p4 <- p4 %>%
            add_trace(
              x = list(grp),
              y = grp_data$avg_degree,
              type = "bar",
              name = grp,
              marker = list(color = group_colors[[grp]]),
              legendgroup = grp,
              showlegend = FALSE,
              hovertemplate = paste0("<b>", grp, "</b><br>Avg Degree: %{y:.2f}<extra></extra>")
            )
        }
        p4 <- p4 %>% layout(
          xaxis = list(title = "Group"),
          yaxis = list(title = "Avg Degree"),
          title = list(text = "Average Degree", font = list(size = 12))
        )

        # 5. Clustering (Transitivity)
        p5 <- plot_ly()
        for (grp in groups) {
          grp_data <- metrics_df %>% filter(group == grp)
          p5 <- p5 %>%
            add_trace(
              x = list(grp),
              y = grp_data$transitivity,
              type = "bar",
              name = grp,
              marker = list(color = group_colors[[grp]]),
              legendgroup = grp,
              showlegend = FALSE,
              hovertemplate = paste0("<b>", grp, "</b><br>Clustering: %{y:.4f}<extra></extra>")
            )
        }
        p5 <- p5 %>% layout(
          xaxis = list(title = "Group"),
          yaxis = list(title = "Clustering"),
          title = list(text = "Clustering Coefficient", font = list(size = 12))
        )

        # 6. Modularity
        p6 <- plot_ly()
        for (grp in groups) {
          grp_data <- metrics_df %>% filter(group == grp)
          p6 <- p6 %>%
            add_trace(
              x = list(grp),
              y = grp_data$modularity,
              type = "bar",
              name = grp,
              marker = list(color = group_colors[[grp]]),
              legendgroup = grp,
              showlegend = FALSE,
              hovertemplate = paste0("<b>", grp, "</b><br>Modularity: %{y:.4f}<extra></extra>")
            )
        }
        p6 <- p6 %>% layout(
          xaxis = list(title = "Group"),
          yaxis = list(title = "Modularity"),
          title = list(text = "Modularity", font = list(size = 12))
        )

        # Combine into subplot layout (2 rows x 3 columns)
        p <- subplot(
          p1, p2, p3,
          p4, p5, p6,
          nrows = 2,
          shareX = FALSE,
          shareY = FALSE,
          titleX = TRUE,
          titleY = TRUE,
          margin = 0.05
        ) %>%
          layout(
            title = list(
              text = paste0("Network Metrics Comparison by Group",
                          ifelse(!is.null(tax_rank_reactive()) && tax_rank_reactive() != "none",
                                 paste0("<br><span style='font-size:12px; color:#666;'>",
                                       tax_rank_reactive(), " level aggregation</span>"),
                                 "")),
              font = list(size = 16, family = "Arial, sans-serif"),
              y = 0.98
            ),
            showlegend = TRUE,
            legend = list(
              orientation = "h",
              yanchor = "top",
              y = -0.15,
              xanchor = "center",
              x = 0.5
            ),
            margin = list(t = 60, b = 100),
            barmode = "group"
          )

        message("Faceted subplot created successfully with 6 panels")
        p
      }, error = function(e) {
        message("Error in nodes_edges_plot: ", e$message)
        message("Available columns: ", paste(colnames(metrics_df), collapse = ", "))
        message("Data dimensions: ", nrow(metrics_df), " rows")
        plotly_empty() %>%
          layout(title = list(
            text = paste0("Error creating plot<br><sub>", e$message, "</sub>"),
            font = list(size = 12)
          ))
      })
    })

    # Correlation heatmaps
    output$correlation_heatmaps <- renderUI({
      # Don't use req() - show message instead
      if (!analysis_completed_reactive()) {
        return(div(
          style = "color: #666; font-style: italic; text-align: center; padding: 20px;",
          "Run network analysis to see correlation heatmaps"
        ))
      }

      results <- network_results_reactive()
      if (is.null(results)) {
        return(div(
          style = "color: #666; font-style: italic; text-align: center; padding: 20px;",
          "No network data available"
        ))
      }

      groups <- names(results)[names(results) != "comparison"]

      # Create a heatmap card for each group
      heatmap_cards <- lapply(groups, function(g) {
        card(
          full_screen = TRUE,
          card_header(
            style = "background-color: #2FA4E7; color: white;",
            paste("Correlation Heatmap -", g)
          ),
          plotlyOutput(ns(paste0("heatmap_", g)), height = "500px")
        )
      })

      do.call(tagList, heatmap_cards)
    })

    # Render correlation heatmaps dynamically
    observe({
      # Don't use req() here either
      if (!analysis_completed_reactive()) {
        return(NULL)
      }

      results <- network_results_reactive()
      if (is.null(results)) {
        return(NULL)
      }

      groups <- names(results)[names(results) != "comparison"]

      lapply(groups, function(g) {
        local({
          group_name <- g
          output[[paste0("heatmap_", group_name)]] <- renderPlotly({

            cor_mat <- results[[group_name]]$correlation_matrix
            p_mat <- results[[group_name]]$p_matrix

            # Filter by significance
            cor_mat_filtered <- cor_mat
            cor_mat_filtered[p_mat > 0.05] <- 0

            # Convert to long format
            cor_df <- as.data.frame(cor_mat_filtered)
            cor_df$taxa1 <- rownames(cor_df)
            cor_long <- cor_df %>%
              pivot_longer(cols = -taxa1, names_to = "taxa2", values_to = "correlation")

            # Remove zero correlations for cleaner visualization
            cor_long <- cor_long %>%
              filter(correlation != 0)

            # Create heatmap
            plot_ly(
              data = cor_long,
              x = ~taxa2,
              y = ~taxa1,
              z = ~correlation,
              type = "heatmap",
              colors = colorRamp(c("blue", "white", "red")),
              zmid = 0,
              zmin = -1,
              zmax = 1,
              colorbar = list(title = "Correlation")
            ) %>%
              layout(
                title = paste("Significant Correlations -", group_name),
                xaxis = list(title = "", showticklabels = FALSE),
                yaxis = list(title = "", showticklabels = FALSE)
              )
          })
        })
      })
    })

    # Edge count info display
    output$edge_count_info <- renderUI({
      if (!analysis_completed_reactive()) {
        return(HTML("<i>Run network analysis to see edge counts</i>"))
      }

      results <- network_results_reactive()
      if (is.null(results)) {
        return(HTML("<i>No data available</i>"))
      }

      tryCatch({
        ps <- filtered_ps_reactive()
        tax_table <- as.data.frame(as(tax_table(ps), "matrix"))

        groups <- names(results)[names(results) != "comparison"]
        total_edges <- 0

        for (group_name in groups) {
          group_data <- results[[group_name]]
          if (!is.null(group_data$correlation_matrix)) {
            cor_mat <- group_data$correlation_matrix
            edge_list <- which(cor_mat != 0, arr.ind = TRUE)
            edge_list <- edge_list[edge_list[,1] < edge_list[,2], ]
            total_edges <- total_edges + nrow(edge_list)
          }
        }

        HTML(paste0(
          "<b>Total edges across all groups:</b> ", total_edges, "<br>",
          "<i>Showing top ", input$edge_top_n, " taxa ranked by total edge count</i>"
        ))
      }, error = function(e) {
        HTML("<i>Error calculating edge counts</i>")
      })
    })

    # Edge Polarity Plot - Side-by-side groups with mirrored positive/negative
    output$edge_polarity_plot <- renderPlotly({
      if (!analysis_completed_reactive()) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Run network analysis first",
                  font = list(size = 14)
                )))
      }

      results <- network_results_reactive()
      req(!is.null(results))

      tryCatch({
        # Get taxonomy
        ps <- filtered_ps_reactive()
        tax_table <- as.data.frame(as(tax_table(ps), "matrix"))

        # Process both groups
        groups <- names(results)[names(results) != "comparison"]
        all_edge_data <- data.frame()

        for (group_name in groups) {
          group_data <- results[[group_name]]
          if (is.null(group_data$correlation_matrix)) next

          cor_mat <- group_data$correlation_matrix

          # Debug: Check raw correlation matrix
          message("=== EDGE POLARITY DEBUG ===")
          message("Group: ", group_name)
          message("Matrix range: [", min(cor_mat), ", ", max(cor_mat), "]")
          message("Positive entries: ", sum(cor_mat > 0))
          message("Negative entries: ", sum(cor_mat < 0))
          message("===========================")

          # Count positive and negative edges by taxonomy
          edge_list <- which(cor_mat != 0, arr.ind = TRUE)
          edge_list <- edge_list[edge_list[,1] < edge_list[,2], ]  # Remove duplicates

          if (nrow(edge_list) == 0) next

          # Extract correlations properly
          correlations <- numeric(nrow(edge_list))
          for (i in 1:nrow(edge_list)) {
            correlations[i] <- cor_mat[edge_list[i, 1], edge_list[i, 2]]
          }

          edge_df <- data.frame(
            node1 = rownames(cor_mat)[edge_list[,1]],
            node2 = rownames(cor_mat)[edge_list[,2]],
            correlation = correlations,
            stringsAsFactors = FALSE
          )

          # Debug: Check correlation signs
          n_pos <- sum(edge_df$correlation > 0)
          n_neg <- sum(edge_df$correlation < 0)
          message("Edge polarity extraction for ", group_name, ": ", n_pos, " positive, ", n_neg, " negative")

          edge_df$tax1 <- tax_table[[input$edge_tax_rank]][match(edge_df$node1, rownames(tax_table))]
          edge_df$tax2 <- tax_table[[input$edge_tax_rank]][match(edge_df$node2, rownames(tax_table))]
          edge_df$tax1[is.na(edge_df$tax1)] <- "Unknown"
          edge_df$tax2[is.na(edge_df$tax2)] <- "Unknown"

          # Count by taxonomy - count both endpoints of each edge
          # This way, if edge connects Bacteroidetes-Firmicutes, both get counted
          edge_df <- edge_df %>%
            mutate(edge_type = ifelse(correlation > 0, "Positive", "Negative"))

          # Count edges for each taxonomy (both endpoints)
          # An edge between A and B contributes 1 to both A's count and B's count
          edge_summary_tax1 <- edge_df %>%
            group_by(taxonomy = tax1, edge_type) %>%
            summarise(count = n(), .groups = "drop")

          edge_summary_tax2 <- edge_df %>%
            group_by(taxonomy = tax2, edge_type) %>%
            summarise(count = n(), .groups = "drop")

          # Combine and aggregate (edges within same taxonomy counted once, not twice)
          edge_summary <- bind_rows(edge_summary_tax1, edge_summary_tax2) %>%
            group_by(taxonomy, edge_type) %>%
            summarise(count = sum(count), .groups = "drop") %>%
            mutate(group = group_name)

          all_edge_data <- rbind(all_edge_data, edge_summary)
        }

        if (nrow(all_edge_data) == 0) {
          return(plotly_empty() %>%
                  layout(title = "No edges found"))
        }

        # Calculate total edges per taxonomy across all groups
        taxa_totals <- all_edge_data %>%
          group_by(taxonomy) %>%
          summarise(total_edges = sum(count), .groups = "drop") %>%
          arrange(desc(total_edges))

        # Filter to top N taxa
        top_n <- min(input$edge_top_n, nrow(taxa_totals))
        top_taxa <- taxa_totals$taxonomy[1:top_n]

        all_edge_data <- all_edge_data %>%
          filter(taxonomy %in% top_taxa)

        # Reorder taxonomy by total edge count (descending)
        all_edge_data$taxonomy <- factor(all_edge_data$taxonomy,
                                          levels = top_taxa)

        # Create side-by-side grouped bars with stacked positive/negative
        # Use numeric x positions with offset for each group
        n_groups <- length(groups)
        n_taxa <- length(top_taxa)

        # Generate base colors
        if (n_groups <= 8) {
          base_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                           "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")
        } else {
          base_palette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                           "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                           "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
        }

        # Function to darken a color for negative edges
        darken_color <- function(hex_color, factor = 0.6) {
          rgb_val <- col2rgb(hex_color)
          darkened <- rgb_val * factor
          rgb(darkened[1], darkened[2], darkened[3], maxColorValue = 255)
        }

        # Create color mapping for groups
        group_colors <- setNames(
          lapply(seq_along(groups), function(i) {
            base_col <- base_palette[((i - 1) %% length(base_palette)) + 1]
            list(pos = base_col, neg = darken_color(base_col))
          }),
          groups
        )

        # Calculate bar width and offsets
        bar_width <- 0.8 / n_groups
        offsets <- seq(-0.4 + bar_width/2, 0.4 - bar_width/2, length.out = n_groups)
        names(offsets) <- groups

        p <- plot_ly()

        # Add traces for each group
        for (g in groups) {
          pos_color <- group_colors[[g]]$pos
          neg_color <- group_colors[[g]]$neg
          offset <- offsets[[g]]

          group_data <- all_edge_data %>% filter(group == g)

          # Calculate x positions (numeric with offset)
          x_positions <- seq_along(top_taxa) + offset

          # Positive data
          pos_counts <- sapply(top_taxa, function(tax) {
            val <- group_data %>% filter(taxonomy == tax, edge_type == "Positive") %>% pull(count)
            if (length(val) == 0) 0 else val
          })

          # Negative data (as negative values)
          neg_counts <- sapply(top_taxa, function(tax) {
            val <- group_data %>% filter(taxonomy == tax, edge_type == "Negative") %>% pull(count)
            if (length(val) == 0) 0 else -val
          })

          # Add positive bars
          p <- p %>% add_trace(
            x = x_positions,
            y = pos_counts,
            type = "bar",
            width = bar_width,
            name = paste0(g, " Positive"),
            marker = list(color = pos_color),
            text = paste0(pos_counts, " edges"),
            hovertemplate = paste0("<b>", g, "</b><br>", top_taxa, "<br>Positive: ", pos_counts, " edges<extra></extra>"),
            legendgroup = g
          )

          # Add negative bars
          p <- p %>% add_trace(
            x = x_positions,
            y = neg_counts,
            type = "bar",
            width = bar_width,
            name = paste0(g, " Negative"),
            marker = list(color = neg_color),
            text = paste0(abs(neg_counts), " edges"),
            hovertemplate = paste0("<b>", g, "</b><br>", top_taxa, "<br>Negative: ", abs(neg_counts), " edges<extra></extra>"),
            legendgroup = g,
            showlegend = FALSE
          )
        }

        p <- p %>%
          layout(
            title = list(
              text = paste0("Edge Type Distribution by ", input$edge_tax_rank, " (Top ", top_n, " Taxa)"),
              font = list(size = 14),
              y = 0.95
            ),
            xaxis = list(
              title = list(text = input$edge_tax_rank, standoff = 20),
              tickmode = "array",
              tickvals = seq_along(top_taxa),
              ticktext = top_taxa,
              tickangle = -45
            ),
            yaxis = list(
              title = "Edge Count (Positive ↑ / Negative ↓)",
              zeroline = TRUE,
              zerolinewidth = 2,
              zerolinecolor = "#000"
            ),
            barmode = "overlay",
            showlegend = TRUE,
            legend = list(
              orientation = "h",
              yanchor = "top",
              y = -0.35,
              xanchor = "center",
              x = 0.5
            ),
            margin = list(t = 60, b = 150),
            hovermode = "closest"
          )

        p

      }, error = function(e) {
        message("Error in edge_polarity_plot: ", e$message)
        plotly_empty() %>%
          layout(title = list(
            text = paste0("Error<br><sub>", e$message, "</sub>"),
            font = list(size = 12)
          ))
      })
    })

    # =========================================================================
    # Download Handlers with Full Preview
    # =========================================================================

    # Metrics plot base function
    metrics_base_plot <- function() {
      results <- network_results_reactive()
      if (is.null(results) || is.null(results$comparison)) return(NULL)

      metrics_df <- results$comparison

      plot_data <- metrics_df %>%
        tidyr::pivot_longer(
          cols = c(n_nodes, n_edges, density, avg_degree, transitivity, modularity),
          names_to = "metric",
          values_to = "value"
        ) %>%
        filter(!is.na(value))

      if (nrow(plot_data) == 0) return(NULL)

      # Use proper labels matching plotly display
      metric_labels <- c(
        "n_nodes" = "Nodes",
        "n_edges" = "Edges",
        "density" = "Density",
        "avg_degree" = "Avg Degree",
        "transitivity" = "Transitivity",
        "modularity" = "Modularity"
      )

      plot_data$metric <- factor(plot_data$metric,
                                  levels = names(metric_labels),
                                  labels = metric_labels)

      p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~metric, scales = "free_y", ncol = 3) +
        labs(title = "Network Topology Metrics Comparison",
             x = "Group", y = "Value", fill = "Group")

      return(p)
    }

    # Setup metrics download with full preview
    setup_full_download_handler(
      input, output, session,
      plot_id = "metrics",
      plot_fn = metrics_base_plot,
      filename_prefix = "network_metrics",
      modal_title = "Download Network Metrics Plot"
    )

    # Edge distribution plot base function
    edge_dist_base_plot <- function() {
      results <- network_results_reactive()
      if (is.null(results)) return(NULL)

      groups <- names(results)[!names(results) %in% c("comparison", "braingraph_robustness")]

      all_edge_data <- data.frame()
      for (g in groups) {
        if (is.null(results[[g]]$igraph)) next
        graph <- results[[g]]$igraph
        if (igraph::ecount(graph) == 0) next

        edge_signs <- igraph::E(graph)$sign
        if (is.null(edge_signs)) next

        pos_count <- sum(edge_signs == "positive", na.rm = TRUE)
        neg_count <- sum(edge_signs == "negative", na.rm = TRUE)

        all_edge_data <- rbind(all_edge_data, data.frame(
          group = g,
          positive = pos_count,
          negative = neg_count
        ))
      }

      if (nrow(all_edge_data) == 0) return(NULL)

      plot_data <- all_edge_data %>%
        tidyr::pivot_longer(cols = c(positive, negative),
                            names_to = "edge_type", values_to = "count")

      # Use proper labels matching plotly display
      plot_data$edge_type <- factor(plot_data$edge_type,
                                     levels = c("positive", "negative"),
                                     labels = c("Positive", "Negative"))

      p <- ggplot(plot_data, aes(x = group, y = count, fill = edge_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Positive" = "#2ca02c", "Negative" = "#d62728"),
                          name = "Edge Type") +
        labs(title = "Edge Type Distribution by Group",
             x = "Group", y = "Edge Count")

      return(p)
    }

    # Setup edge distribution download with full preview
    setup_full_download_handler(
      input, output, session,
      plot_id = "edge_dist",
      plot_fn = edge_dist_base_plot,
      filename_prefix = "edge_distribution",
      modal_title = "Download Edge Distribution Plot"
    )

  })
}
