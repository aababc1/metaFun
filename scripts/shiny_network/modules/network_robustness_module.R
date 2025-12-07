# =============================================================================
# Network Robustness Module
# Analyzes network resilience to node removal attacks
# =============================================================================

library(shiny)
library(shinyjs)
library(bslib)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)

#' Network Robustness UI
network_robustness_UI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),

    card(
      full_screen = TRUE,
      class = "main-card",
      card_header(
        style = "background-color: #2FA4E7; color: white;",
        "Network Robustness Analysis"
      ),
      card_body(
        # Message before analysis
        div(
          id = ns("analysis_message"),
          style = "color: red; font-weight: bold; margin-bottom: 20px;",
          icon("exclamation-triangle"),
          " Run network analysis using the Master Controller"
        ),

        # brainGraph Robustness Analysis (New Section)
        h4("brainGraph Network Stability Test"),
        p("Comprehensive attack simulation using brainGraph package (iterations configurable below)"),

        # Settings bar
        div(
          style = "margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px;",
          layout_column_wrap(
            width = 1/3,
            div(
              actionButton(
                ns("run_stability_test"),
                "Run Stability Test (Node Removal)",
                icon = icon("play"),
                class = "btn-primary btn-lg",
                style = "width: 100%; margin-bottom: 10px;"
              ),
              uiOutput(ns("stability_status")),
              div(
                style = "margin-top: 10px; font-size: 12px; color: #6c757d;",
                icon("info-circle"),
                " Attack type: Node (Vertex) removal only"
              )
            ),
            sliderInput(
              ns("n_iterations"),
              "Random Attack Iterations:",
              min = 10,
              max = 500,
              value = 100,
              step = 10
            ),
            checkboxGroupInput(
              ns("braingraph_strategies"),
              "Attack Strategies:",
              choices = c("Degree" = "Degree",
                         "Betweenness" = "Betweenness",
                         "IVI" = "IVI",
                         "Random" = "Random"),
              selected = c("Degree", "Betweenness", "Random"),
              inline = TRUE
            )
          )
        ),

        br(),

        # brainGraph Robustness Plot - full width
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header(
            "Robustness Curves (brainGraph)",
            create_plot_download_btn(ns, "robustness", "Download Plot")
          ),
          plotlyOutput(ns("braingraph_plot"), height = "650px")
        ),

        br(),

        # brainGraph Summary Table
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header("Robustness Metrics Summary (AUC & R50)"),
          DTOutput(ns("braingraph_summary_table"))
        ),

      )
    )
  )
}

#' Network Robustness Server
network_robustness_Server <- function(id, filtered_ps_reactive, network_results_reactive,
                                     analysis_completed_reactive,
                                     reset_trigger_reactive = reactive(NULL)) {
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


    # Reactive value to store brainGraph results
    braingraph_results <- reactiveVal(NULL)

    # Reset stability test results when threshold changes
    observeEvent(reset_trigger_reactive(), {
      if (!is.null(reset_trigger_reactive())) {
        braingraph_results(NULL)
        showNotification(
          "Stability test results cleared. Run the test again with new network.",
          type = "warning",
          duration = 4
        )
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # Status message for stability test
    output$stability_status <- renderUI({
      if (is.null(braingraph_results())) {
        div(
          style = "color: #666; font-size: 12px; margin-top: 10px;",
          icon("info-circle"),
          " Click 'Run Stability Test' to analyze network robustness"
        )
      } else {
        div(
          style = "color: green; font-size: 12px; margin-top: 10px;",
          icon("check-circle"),
          " Stability test completed successfully"
        )
      }
    })

    # Run stability test button
    observeEvent(input$run_stability_test, {
      req(analysis_completed_reactive(), network_results_reactive())

      results <- network_results_reactive()

      showNotification(
        "Running brainGraph stability test... This may take 1-2 minutes",
        type = "message",
        duration = 3
      )

      tryCatch({
        # Extract igraph objects for all groups
        igraph_list <- list()
        ivi_scores_list <- list()
        groups <- names(results)[!names(results) %in% c("comparison", "braingraph_robustness")]

        for (g in groups) {
          if (!is.null(results[[g]]) && !is.null(results[[g]]$igraph)) {
            igraph_list[[g]] <- results[[g]]$igraph

            # Extract IVI scores if available
            if (!is.null(results[[g]]$ivi)) {
              ivi_df <- results[[g]]$ivi
              if (nrow(ivi_df) > 0 && "IVI" %in% colnames(ivi_df)) {
                ivi_scores_list[[g]] <- ivi_df$IVI
              }
            }
          }
        }

        # Check if we have any networks
        if (length(igraph_list) == 0) {
          showNotification(
            "No network objects found in results. Please ensure network analysis completed successfully.",
            type = "error",
            duration = 10
          )
          return()
        }

        # Run brainGraph robustness analysis (Node/Vertex removal only)
        bg_results <- run_network_robustness_analysis(
          igraph_list = igraph_list,
          n_iterations = input$n_iterations,
          centrality_methods = c("random", "degree", "betweenness"),
          optional_ivi_scores = if (length(ivi_scores_list) > 0) ivi_scores_list else NULL,
          edge_removal = FALSE,  # Only node (vertex) removal
          group_names = groups
        )

        braingraph_results(bg_results)

        showNotification(
          "brainGraph stability test completed successfully!",
          type = "message",
          duration = 5
        )

      }, error = function(e) {
        showNotification(
          paste("Error running stability test:", e$message),
          type = "error",
          duration = 10
        )
      })
    })

    # brainGraph Robustness Plot
    output$braingraph_plot <- renderPlotly({
      if (!analysis_completed_reactive()) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Run network analysis first<br><sub>Use Master Controller</sub>",
                  font = list(size = 14)
                )))
      }

      # Check if brainGraph results exist
      if (is.null(braingraph_results())) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Click 'Run Stability Test' to analyze network robustness<br><sub>100 iterations will be performed</sub>",
                  font = list(size = 14, color = "#666")
                )))
      }

      braingraph_data <- braingraph_results()

      if (is.null(braingraph_data$group_results)) {
        return(plotly_empty() %>%
                layout(title = "No group results found"))
      }

      tryCatch({
        # 1. Combine data from ALL groups with preprocessing
        all_plot_data <- data.frame()

        for (grp in names(braingraph_data$group_results)) {
          group_data <- braingraph_data$group_results[[grp]]

          if (!is.null(group_data) && !is.null(group_data$data)) {
            temp_data <- group_data$data

            # (A) Filter: Only Vertex (Node) removal
            # Exception: Include Random attacks even if attack_type is NA or different
            # (we only request vertex removal anyway - see edge_removal = FALSE above)
            if ("attack_type" %in% colnames(temp_data) && "strategy" %in% colnames(temp_data)) {
              temp_data <- temp_data[
                temp_data$attack_type == "Vertex" |
                tolower(temp_data$strategy) == "random",
              ]
            } else if ("attack_type" %in% colnames(temp_data)) {
              temp_data <- temp_data[temp_data$attack_type == "Vertex", ]
            }

            # (B) [KEY FIX] Convert strategy names to Title Case (random -> Random)
            # This ensures matching with UI checkbox values ("Random", "Degree", etc.)
            if ("strategy" %in% colnames(temp_data)) {
              temp_data$strategy <- tools::toTitleCase(as.character(temp_data$strategy))
            }

            temp_data$group <- grp
            all_plot_data <- rbind(all_plot_data, temp_data)
          }
        }

        # Debug: Print available strategies (for verification in console)
        print(paste("Available strategies in data:", paste(unique(all_plot_data$strategy), collapse=", ")))
        print(paste("Selected strategies from UI:", paste(input$braingraph_strategies, collapse=", ")))

        # 2. Filter by UI-selected strategies
        if (nrow(all_plot_data) > 0) {
          all_plot_data <- all_plot_data[all_plot_data$strategy %in% input$braingraph_strategies, ]
        }

        if (nrow(all_plot_data) == 0) {
          return(plotly_empty() %>%
                  layout(title = "No data matches selected strategies"))
        }

        # 3. [KEY FIX] Aggregate Random data (mean across 100 iterations)
        # Random has many iterations, so we average by fraction_removed
        aggregated_data <- all_plot_data %>%
          group_by(group, strategy, fraction_removed) %>%
          summarise(
            largest_component = mean(largest_component, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          arrange(group, strategy, fraction_removed)

        print(paste("Aggregated data rows:", nrow(aggregated_data)))
        print(paste("Final strategies:", paste(unique(aggregated_data$strategy), collapse=", ")))

        # 4. Create plotly plot
        p <- plot_ly()

        # Get unique groups and strategies
        unique_groups <- unique(aggregated_data$group)
        n_groups <- length(unique_groups)

        # Define strategy order: targeted attacks first, then random baseline
        strategy_order <- c("Degree", "Betweenness", "IVI", "Random")
        unique_strategies <- strategy_order[strategy_order %in% unique(aggregated_data$strategy)]

        # COLOR BY STRATEGY (same attack type = same color across groups)
        strategy_colors <- c(
          "Random" = "#7f7f7f",       # Gray - baseline
          "Degree" = "#d62728",       # Red - degree attack
          "Betweenness" = "#1f77b4",  # Blue - betweenness attack
          "IVI" = "#2ca02c"           # Green - IVI attack
        )

        # LINE STYLE BY GROUP (different groups = different dash patterns)
        # solid, dash, dot, dashdot, longdash, longdashdot
        line_styles <- c("solid", "dash", "dot", "dashdot", "longdash", "longdashdot")
        group_dashes <- setNames(
          line_styles[((seq_along(unique_groups) - 1) %% length(line_styles)) + 1],
          unique_groups
        )

        # MARKER SHAPE BY GROUP (for better distinction)
        marker_shapes <- c("circle", "square", "diamond", "triangle-up", "cross", "x")
        group_markers <- setNames(
          marker_shapes[((seq_along(unique_groups) - 1) %% length(marker_shapes)) + 1],
          unique_groups
        )

        # Plot each strategy + group combination (strategy first for legend grouping)
        for (strat in unique_strategies) {
          strat_col <- if(strat %in% names(strategy_colors)) strategy_colors[[strat]] else "#333333"

          for (grp in unique_groups) {
            plot_subset <- aggregated_data[aggregated_data$group == grp &
                                          aggregated_data$strategy == strat, ]

            if (nrow(plot_subset) > 0) {
              line_dash <- group_dashes[[grp]]
              marker_sym <- group_markers[[grp]]

              p <- p %>% add_trace(
                data = plot_subset,
                x = ~fraction_removed,
                y = ~largest_component,
                name = paste0(strat, " (", grp, ")"),
                legendgroup = strat,  # Group by strategy in legend
                type = "scatter",
                mode = "lines+markers",
                line = list(color = strat_col, width = 2.5, dash = line_dash),
                marker = list(size = 7, symbol = marker_sym, color = strat_col),
                hovertemplate = paste0(
                  "<b>", strat, " - ", grp, "</b><br>",
                  "Removed: %{x:.1%}<br>",
                  "LCC Size: %{y:.1%}<br>",
                  "<extra></extra>"
                )
              )
            }
          }
        }

        p %>% layout(
          title = list(
            text = "<b>Network Robustness Comparison - Node (Vertex) Removal</b>",
            font = list(size = 16)
          ),
          xaxis = list(title = "Fraction of Nodes Removed"),
          yaxis = list(title = "Largest Component Size (fraction)"),
          hovermode = "closest",
          legend = list(x = 1.02, y = 1, xanchor = "left"),
          plot_bgcolor = "#f8f9fa"
        )

      }, error = function(e) {
        return(plotly_empty() %>%
                layout(title = paste("Error:", e$message)))
      })
    })

    # brainGraph Summary Table
    output$braingraph_summary_table <- renderDT({
      if (!analysis_completed_reactive()) {
        return(datatable(
          data.frame(Message = "Run network analysis first using the Master Controller"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      # Use the braingraph_results reactive
      braingraph_data <- braingraph_results()

      if (is.null(braingraph_data) || is.null(braingraph_data$summary_table)) {
        return(datatable(
          data.frame(Message = "Click 'Run Stability Test' button above to generate robustness metrics"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      tryCatch({
        summary_df <- braingraph_data$summary_table

        if (nrow(summary_df) == 0) {
          return(datatable(
            data.frame(Message = "No robustness data available - try running the stability test again"),
            options = list(dom = 't'),
            rownames = FALSE
          ))
        }

        datatable(
          summary_df,
          extensions = "Buttons",
          options = list(
            pageLength = 20,
            dom = "Bfrtip",
            buttons = c("copy", "csv", "excel"),
            scrollX = TRUE
          ),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("AUC", "R50"), digits = 3)

      }, error = function(e) {
        return(datatable(
          data.frame(Error = paste("Error:", e$message)),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      })
    })

    # =========================================================================
    # Download Handlers with Full Preview
    # =========================================================================

    # Robustness plot data reactive
    robustness_plot_data <- reactive({
      braingraph_data <- braingraph_results()
      if (is.null(braingraph_data) || is.null(braingraph_data$group_results)) {
        return(NULL)
      }

      # Combine all data
      all_plot_data <- data.frame()
      for (grp in names(braingraph_data$group_results)) {
        group_data <- braingraph_data$group_results[[grp]]
        if (!is.null(group_data) && !is.null(group_data$data)) {
          temp_data <- group_data$data
          if ("strategy" %in% colnames(temp_data)) {
            temp_data$strategy <- tools::toTitleCase(as.character(temp_data$strategy))
          }
          temp_data$group <- grp
          all_plot_data <- rbind(all_plot_data, temp_data)
        }
      }

      if (nrow(all_plot_data) > 0) {
        all_plot_data <- all_plot_data[all_plot_data$strategy %in% input$braingraph_strategies, ]
      }

      all_plot_data
    })

    # Robustness plot base function
    robustness_base_plot <- function() {
      data <- robustness_plot_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)

      agg_data <- data %>%
        dplyr::group_by(group, strategy, fraction_removed) %>%
        dplyr::summarise(largest_component = mean(largest_component, na.rm = TRUE), .groups = 'drop')

      if (nrow(agg_data) == 0) return(NULL)

      # Use colors matching plotly display
      strategy_colors <- c(
        "Random" = "#7f7f7f",
        "Degree" = "#d62728",
        "Betweenness" = "#1f77b4",
        "Ivi" = "#2ca02c"
      )

      p <- ggplot(agg_data, aes(x = fraction_removed, y = largest_component,
                                color = strategy, linetype = group)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_color_manual(values = strategy_colors, name = "Attack Strategy") +
        labs(
          title = "Network Robustness Comparison - Node Removal",
          x = "Fraction of Nodes Removed",
          y = "Largest Component Size (fraction)",
          linetype = "Group"
        )

      return(p)
    }

    # Setup robustness download with full preview
    setup_full_download_handler(
      input, output, session,
      plot_id = "robustness",
      plot_fn = robustness_base_plot,
      filename_prefix = "robustness_curves",
      modal_title = "Download Robustness Plot"
    )

  })
}
