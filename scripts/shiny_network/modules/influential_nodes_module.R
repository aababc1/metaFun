# =============================================================================
# Influential Nodes Module
# Analyzes key network players using Zi-Pi plots and hub metrics
# =============================================================================

library(shiny)
library(shinyjs)
library(bslib)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)

#' Influential Nodes UI
influential_nodes_UI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),

    card(
      full_screen = TRUE,
      class = "main-card",
      card_header(
        style = "background-color: #2FA4E7; color: white;",
        "Influential Nodes Analysis"
      ),
      card_body(
        # Message before analysis
        div(
          id = ns("analysis_message"),
          style = "color: red; font-weight: bold; margin-bottom: 20px;",
          icon("exclamation-triangle"),
          " Run network analysis using the Master Controller"
        ),

        # Zi-Pi Plot Section
        h4("Module Roles: Zi-Pi Analysis"),
        p("Identifies network roles based on within-module (Zi) and among-module (Pi) connectivity"),

        # Horizontal control bar (no card wrapper)
        div(
          style = "margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px;",
          layout_column_wrap(
            width = 1/3,
            selectInput(
              ns("tax_rank"),
              "Color by Taxonomic Rank:",
              choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
              selected = "Species"
            ),
            selectInput(
              ns("group_select"),
              "Select Group:",
              choices = NULL
            ),
            checkboxInput(
              ns("show_labels"),
              "Show node labels (top 10 by IVI)",
              value = FALSE
            )
          )
        ),

        # Full-width Zi-Pi plot
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header(
            "Zi-Pi Plot (Point size = IVI)",
            create_plot_download_btn(ns, "zipi", "Download")
          ),
          plotlyOutput(ns("zipi_plot"), height = "650px")
        ),

        br(),

        # Differential Zi-Pi Plot Section
        h4("Differential Zi-Pi Analysis"),
        p("Compare how taxa network roles change between two groups"),

        # Group comparison controls
        div(
          style = "margin-bottom: 20px; padding: 15px; background-color: #fff3cd; border-radius: 5px; border-left: 4px solid #ffc107;",
          layout_column_wrap(
            width = 1/3,
            selectInput(
              ns("diff_group_a"),
              "Group A:",
              choices = NULL
            ),
            selectInput(
              ns("diff_group_b"),
              "Group B:",
              choices = NULL
            ),
            selectInput(
              ns("diff_tax_rank"),
              "Color by Taxonomy:",
              choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
              selected = "Species"
            )
          ),
          div(
            style = "margin-top: 10px; font-size: 12px; color: #856404;",
            icon("info-circle"),
            HTML(" <strong>Interpretation:</strong> ΔZi = Zi(A) - Zi(B), ΔPi = Pi(A) - Pi(B).
                 Positive values = higher in Group A, Negative = higher in Group B")
          )
        ),

        # Differential Zi-Pi plot
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header(
            "Differential Zi-Pi Plot: Group A vs Group B",
            create_plot_download_btn(ns, "diff_zipi", "Download")
          ),
          plotlyOutput(ns("diff_zipi_plot"), height = "650px")
        ),

        br(),

        # Top Influential Nodes Table
        h4("Top Influential Nodes"),
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header("Hub Metrics (sortable)"),
          DTOutput(ns("hub_table"))
        ),

        br(),

        # Group-Specific Taxa Section
        h4("Group-Specific Taxa Analysis"),
        p("Identify taxa that are unique to each group's network (present in one network but not the other)"),

        # Group comparison controls for unique taxa
        div(
          style = "margin-bottom: 20px; padding: 15px; background-color: #e8f5e9; border-radius: 5px; border-left: 4px solid #4caf50;",
          layout_column_wrap(
            width = 1/2,
            selectInput(
              ns("unique_group_a"),
              "Group A:",
              choices = NULL
            ),
            selectInput(
              ns("unique_group_b"),
              "Group B:",
              choices = NULL
            )
          ),
          div(
            style = "margin-top: 10px; font-size: 12px; color: #2e7d32;",
            icon("info-circle"),
            HTML(" <strong>Note:</strong> Shows taxa present in one group's network but absent in the other.
                 This reveals group-specific microbial associations.")
          )
        ),

        # Network Taxa Overlap
        card(
          fill = FALSE,
          card_header("Network Taxa Overlap"),
          card_body(
            plotlyOutput(ns("venn_plot"), height = "350px")
          )
        ),

        br(),

        # Summary Statistics
        card(
          fill = FALSE,
          card_header("Summary Statistics"),
          card_body(
            uiOutput(ns("venn_summary"))
          )
        ),

        br(),

        # Group-specific taxa table
        card(
          fill = FALSE,
          card_header("Group-Specific Taxa Details"),
          card_body(
            DTOutput(ns("unique_taxa_table"))
          )
        )

      )
    )
  )
}

#' Influential Nodes Server
influential_nodes_Server <- function(id, filtered_ps_reactive, network_results_reactive,
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

    # Update taxonomic rank selectors based on Master Controller selection
    observe({
      req(tax_rank_reactive())
      master_rank <- tax_rank_reactive()

      # Get allowed ranks (master rank and above only)
      idx <- which(tax_ranks == master_rank)
      if (length(idx) == 0) idx <- length(tax_ranks)
      allowed_ranks <- tax_ranks[1:idx]

      # Update tax_rank selectInput
      updateSelectInput(session, "tax_rank",
                       choices = allowed_ranks,
                       selected = master_rank)

      # Update diff_tax_rank selectInput
      updateSelectInput(session, "diff_tax_rank",
                       choices = allowed_ranks,
                       selected = master_rank)
    })

    # Update group selectors
    observe({
      req(analysis_completed_reactive(), network_results_reactive())

      results <- network_results_reactive()
      groups <- names(results)[names(results) != "comparison"]

      updateSelectInput(session, "group_select", choices = groups, selected = groups[1])

      # Update differential group selectors
      if (length(groups) >= 2) {
        updateSelectInput(session, "diff_group_a", choices = groups, selected = groups[1])
        updateSelectInput(session, "diff_group_b", choices = groups, selected = groups[2])
        # Update unique taxa group selectors
        updateSelectInput(session, "unique_group_a", choices = groups, selected = groups[1])
        updateSelectInput(session, "unique_group_b", choices = groups, selected = groups[2])
      }
    })

    # Zi-Pi Plot
    output$zipi_plot <- renderPlotly({
      if (!analysis_completed_reactive()) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Run network analysis first<br><sub>Use Master Controller</sub>",
                  font = list(size = 14)
                )))
      }

      results <- network_results_reactive()
      req(!is.null(results), input$group_select)

      group_data <- results[[input$group_select]]
      req(!is.null(group_data$zi_pi), !is.null(group_data$ivi))

      # Check if data frames have rows (network may have 0 edges)
      if (nrow(group_data$zi_pi) == 0 || nrow(group_data$ivi) == 0) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "No network edges found<br><sub>Try lowering correlation threshold (r) or increasing p-value threshold</sub>",
                  font = list(size = 14)
                )))
      }

      tryCatch({
        # Get Zi-Pi data
        zipi_df <- group_data$zi_pi

        # Get taxonomy
        ps <- filtered_ps_reactive()
        tax_table <- as.data.frame(as(tax_table(ps), "matrix"))

        # Match taxonomy to nodes
        zipi_df$taxonomy <- tax_table[[input$tax_rank]][match(zipi_df$node, rownames(tax_table))]
        zipi_df$taxonomy[is.na(zipi_df$taxonomy)] <- "Unknown"

        # Add IVI scores
        ivi_df <- group_data$ivi
        zipi_df$IVI <- ivi_df$IVI[match(zipi_df$node, ivi_df$node)]

        # Create plot
        p <- plot_ly(
          data = zipi_df,
          x = ~Pi,
          y = ~Zi,
          type = "scatter",
          mode = "markers",
          color = ~taxonomy,
          size = ~IVI,
          sizes = c(10, 100),
          text = ~paste0(
            "<b>Node:</b> ", node, "<br>",
            "<b>Taxonomy:</b> ", taxonomy, "<br>",
            "<b>Zi:</b> ", round(Zi, 3), "<br>",
            "<b>Pi:</b> ", round(Pi, 3), "<br>",
            "<b>IVI:</b> ", round(IVI, 2), "<br>",
            "<b>Role:</b> ", role
          ),
          hoverinfo = "text",
          marker = list(
            opacity = 0.7,
            line = list(width = 1, color = "white")
          )
        ) %>%
          layout(
            title = list(
              text = paste0("Module Roles - ", input$group_select),
              font = list(size = 14)
            ),
            xaxis = list(
              title = "Participation Coefficient (Pi)",
              zeroline = TRUE,
              zerolinewidth = 2
            ),
            yaxis = list(
              title = "Within-Module Degree (Zi)",
              zeroline = TRUE,
              zerolinewidth = 2
            ),
            shapes = list(
              # Vertical line at Pi = 0.62
              list(
                type = "line",
                x0 = 0.62, x1 = 0.62,
                y0 = -3, y1 = 3,
                line = list(dash = "dash", color = "gray", width = 1)
              ),
              # Horizontal line at Zi = 2.5
              list(
                type = "line",
                x0 = 0, x1 = 1,
                y0 = 2.5, y1 = 2.5,
                line = list(dash = "dash", color = "gray", width = 1)
              )
            ),
            annotations = list(
              list(x = 0.3, y = 2.8, text = "Module Hubs", showarrow = FALSE, font = list(size = 10, color = "gray")),
              list(x = 0.8, y = 2.8, text = "Network Hubs", showarrow = FALSE, font = list(size = 10, color = "gray")),
              list(x = 0.3, y = -2.5, text = "Peripherals", showarrow = FALSE, font = list(size = 10, color = "gray")),
              list(x = 0.8, y = -2.5, text = "Connectors", showarrow = FALSE, font = list(size = 10, color = "gray"))
            )
          )

        # Add labels for top nodes if requested
        if (input$show_labels) {
          top_nodes <- zipi_df %>% arrange(desc(IVI)) %>% head(10)
          p <- p %>% add_text(
            data = top_nodes,
            x = ~Pi,
            y = ~Zi,
            text = ~node,
            textposition = "top center",
            textfont = list(size = 8),
            showlegend = FALSE
          )
        }

        p

      }, error = function(e) {
        message("Error in zipi_plot: ", e$message)
        plotly_empty() %>%
          layout(title = list(
            text = paste0("Error creating Zi-Pi plot<br><sub>", e$message, "</sub>"),
            font = list(size = 12)
          ))
      })
    })

    # Differential Zi-Pi Plot
    output$diff_zipi_plot <- renderPlotly({
      if (!analysis_completed_reactive()) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Run network analysis first<br><sub>Use Master Controller</sub>",
                  font = list(size = 14)
                )))
      }

      results <- network_results_reactive()
      req(!is.null(results), input$diff_group_a, input$diff_group_b)

      # Validate that two different groups are selected
      if (input$diff_group_a == input$diff_group_b) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Please select two different groups to compare",
                  font = list(size = 14)
                )))
      }

      tryCatch({
        # Get Zi-Pi data for both groups
        group_a_data <- results[[input$diff_group_a]]
        group_b_data <- results[[input$diff_group_b]]

        req(!is.null(group_a_data$zi_pi), !is.null(group_b_data$zi_pi))

        zipi_a <- group_a_data$zi_pi
        zipi_b <- group_b_data$zi_pi

        # Check if data frames have rows
        if (nrow(zipi_a) == 0 || nrow(zipi_b) == 0) {
          return(plotly_empty() %>%
                  layout(title = list(
                    text = "One or both groups have no network edges<br><sub>Adjust thresholds to get more edges</sub>",
                    font = list(size = 14)
                  )))
        }

        # Merge by node (taxa present in both groups)
        diff_df <- merge(
          zipi_a[, c("node", "Zi", "Pi")],
          zipi_b[, c("node", "Zi", "Pi")],
          by = "node",
          suffixes = c("_A", "_B")
        )

        # Calculate differences: A - B
        diff_df$delta_Zi <- diff_df$Zi_A - diff_df$Zi_B
        diff_df$delta_Pi <- diff_df$Pi_A - diff_df$Pi_B

        # Get taxonomy
        ps <- filtered_ps_reactive()
        tax_table <- as.data.frame(as(tax_table(ps), "matrix"))

        # Match taxonomy to nodes
        diff_df$taxonomy <- tax_table[[input$diff_tax_rank]][match(diff_df$node, rownames(tax_table))]
        diff_df$taxonomy[is.na(diff_df$taxonomy)] <- "Unknown"

        # Calculate magnitude of change
        diff_df$magnitude <- sqrt(diff_df$delta_Zi^2 + diff_df$delta_Pi^2)

        # Create plot
        p <- plot_ly(
          data = diff_df,
          x = ~delta_Pi,
          y = ~delta_Zi,
          type = "scatter",
          mode = "markers",
          color = ~taxonomy,
          size = ~magnitude,
          sizes = c(10, 100),
          text = ~paste0(
            "<b>Taxon:</b> ", node, "<br>",
            "<b>Taxonomy:</b> ", taxonomy, "<br>",
            "<b>ΔZi:</b> ", round(delta_Zi, 3), " (", input$diff_group_a, " - ", input$diff_group_b, ")<br>",
            "<b>ΔPi:</b> ", round(delta_Pi, 3), " (", input$diff_group_a, " - ", input$diff_group_b, ")<br>",
            "<b>Zi (", input$diff_group_a, "):</b> ", round(Zi_A, 3), "<br>",
            "<b>Zi (", input$diff_group_b, "):</b> ", round(Zi_B, 3), "<br>",
            "<b>Pi (", input$diff_group_a, "):</b> ", round(Pi_A, 3), "<br>",
            "<b>Pi (", input$diff_group_b, "):</b> ", round(Pi_B, 3), "<br>",
            "<b>Magnitude:</b> ", round(magnitude, 3)
          ),
          hoverinfo = "text",
          marker = list(
            opacity = 0.7,
            line = list(width = 1, color = "white")
          )
        ) %>%
          layout(
            title = list(
              text = paste0("Differential Zi-Pi: ", input$diff_group_a, " vs ", input$diff_group_b,
                           "<br><sub>", nrow(diff_df), " taxa present in both groups</sub>"),
              font = list(size = 14)
            ),
            xaxis = list(
              title = paste0("ΔPi = Pi(", input$diff_group_a, ") - Pi(", input$diff_group_b, ")"),
              zeroline = TRUE,
              zerolinewidth = 2,
              zerolinecolor = "black"
            ),
            yaxis = list(
              title = paste0("ΔZi = Zi(", input$diff_group_a, ") - Zi(", input$diff_group_b, ")"),
              zeroline = TRUE,
              zerolinewidth = 2,
              zerolinecolor = "black"
            ),
            shapes = list(
              # Vertical line at ΔPi = 0
              list(
                type = "line",
                x0 = 0, x1 = 0,
                y0 = -5, y1 = 5,
                line = list(dash = "solid", color = "black", width = 2)
              ),
              # Horizontal line at ΔZi = 0
              list(
                type = "line",
                x0 = -5, x1 = 5,
                y0 = 0, y1 = 0,
                line = list(dash = "solid", color = "black", width = 2)
              )
            ),
            annotations = list(
              # Quadrant labels
              list(
                x = 0.8, y = 0.8, xref = "paper", yref = "paper",
                text = paste0("Higher Zi & Pi<br>in ", input$diff_group_a),
                showarrow = FALSE,
                font = list(size = 11, color = "darkgreen"),
                bgcolor = "rgba(144, 238, 144, 0.3)",
                borderpad = 5
              ),
              list(
                x = 0.2, y = 0.8, xref = "paper", yref = "paper",
                text = paste0("Higher Zi in ", input$diff_group_a, "<br>Higher Pi in ", input$diff_group_b),
                showarrow = FALSE,
                font = list(size = 11, color = "darkorange"),
                bgcolor = "rgba(255, 165, 0, 0.2)",
                borderpad = 5
              ),
              list(
                x = 0.2, y = 0.2, xref = "paper", yref = "paper",
                text = paste0("Higher Zi & Pi<br>in ", input$diff_group_b),
                showarrow = FALSE,
                font = list(size = 11, color = "darkred"),
                bgcolor = "rgba(255, 99, 71, 0.2)",
                borderpad = 5
              ),
              list(
                x = 0.8, y = 0.2, xref = "paper", yref = "paper",
                text = paste0("Higher Pi in ", input$diff_group_a, "<br>Higher Zi in ", input$diff_group_b),
                showarrow = FALSE,
                font = list(size = 11, color = "darkorange"),
                bgcolor = "rgba(255, 165, 0, 0.2)",
                borderpad = 5
              )
            )
          )

        p

      }, error = function(e) {
        message("Error in diff_zipi_plot: ", e$message)
        plotly_empty() %>%
          layout(title = list(
            text = paste0("Error creating differential Zi-Pi plot<br><sub>", e$message, "</sub>"),
            font = list(size = 12)
          ))
      })
    })

    # Top Hub Nodes Table
    output$hub_table <- renderDT({
      if (!analysis_completed_reactive()) {
        return(datatable(
          data.frame(Message = "Run network analysis using the Master Controller"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      results <- network_results_reactive()
      req(!is.null(results), input$group_select)

      group_data <- results[[input$group_select]]
      req(!is.null(group_data$node_metrics), !is.null(group_data$ivi))

      # Check if data frames have rows
      if (nrow(group_data$node_metrics) == 0 || nrow(group_data$ivi) == 0) {
        return(datatable(
          data.frame(Message = "No network edges found. Try lowering correlation threshold (r) or increasing p-value threshold."),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      tryCatch({
        # Combine node metrics
        node_df <- group_data$node_metrics
        ivi_df <- group_data$ivi
        zipi_df <- group_data$zi_pi

        # Get taxonomy
        ps <- filtered_ps_reactive()
        tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
        selected_rank <- tax_rank_reactive()

        # Helper function to safely get taxonomy column
        get_tax_col <- function(rank) {
          if (rank %in% colnames(tax_df)) {
            vals <- tax_df[[rank]][match(node_df$node, rownames(tax_df))]
            vals[is.na(vals)] <- "Unknown"
            return(vals)
          }
          return(rep("Unknown", nrow(node_df)))
        }

        hub_df <- data.frame(
          Node = node_df$node,
          Group = input$group_select,
          Phylum = get_tax_col("Phylum"),
          stringsAsFactors = FALSE
        )

        # Add selected rank if different from Phylum
        if (selected_rank != "Phylum") {
          hub_df[[selected_rank]] <- get_tax_col(selected_rank)
        }

        # Add metrics
        hub_df$Degree <- node_df$degree
        hub_df$Betweenness <- node_df$betweenness
        hub_df$Closeness <- node_df$closeness
        hub_df$IVI <- ivi_df$IVI[match(node_df$node, ivi_df$node)]
        hub_df$Zi <- zipi_df$Zi[match(node_df$node, zipi_df$node)]
        hub_df$Pi <- zipi_df$Pi[match(node_df$node, zipi_df$node)]
        hub_df$Role <- zipi_df$role[match(node_df$node, zipi_df$node)]

        # Sort by IVI descending
        hub_df <- hub_df %>% arrange(desc(IVI))

        # Calculate IVI column index dynamically
        ivi_col_idx <- which(colnames(hub_df) == "IVI") - 1  # 0-based

        datatable(
          hub_df,
          extensions = 'Buttons',
          options = list(
            pageLength = 20,
            scrollX = TRUE,
            order = list(list(ivi_col_idx, 'desc')),
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'csv', filename = paste0('hub_metrics_', input$group_select)),
              list(extend = 'excel', filename = paste0('hub_metrics_', input$group_select))
            )
          ),
          rownames = FALSE,
          filter = 'top'
        ) %>%
          formatRound(columns = c("Betweenness", "Closeness", "IVI", "Zi", "Pi"), digits = 3)

      }, error = function(e) {
        message("Error in hub_table: ", e$message)
        datatable(
          data.frame(Error = paste("Error:", e$message)),
          options = list(dom = 't'),
          rownames = FALSE
        )
      })
    })

    # =========================================================================
    # Group-Specific Taxa Analysis
    # =========================================================================

    # Reactive: Calculate unique and shared taxa between groups
    taxa_comparison <- reactive({
      req(analysis_completed_reactive(), network_results_reactive())
      req(input$unique_group_a, input$unique_group_b)

      if (input$unique_group_a == input$unique_group_b) {
        return(NULL)
      }

      results <- network_results_reactive()
      group_a_data <- results[[input$unique_group_a]]
      group_b_data <- results[[input$unique_group_b]]

      req(!is.null(group_a_data$igraph), !is.null(group_b_data$igraph))

      # Get node names from each network
      nodes_a <- V(group_a_data$igraph)$name
      nodes_b <- V(group_b_data$igraph)$name

      # Calculate set operations
      unique_a <- setdiff(nodes_a, nodes_b)
      unique_b <- setdiff(nodes_b, nodes_a)
      shared <- intersect(nodes_a, nodes_b)

      list(
        unique_a = unique_a,
        unique_b = unique_b,
        shared = shared,
        all_a = nodes_a,
        all_b = nodes_b,
        group_a_name = input$unique_group_a,
        group_b_name = input$unique_group_b
      )
    })

    # Venn Diagram Plot
    output$venn_plot <- renderPlotly({
      if (!analysis_completed_reactive()) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Run network analysis first",
                  font = list(size = 14)
                )))
      }

      comp <- taxa_comparison()
      if (is.null(comp)) {
        return(plotly_empty() %>%
                layout(title = list(
                  text = "Select two different groups to compare",
                  font = list(size = 14)
                )))
      }

      # Create a simple Venn-like visualization using plotly shapes
      n_a_only <- length(comp$unique_a)
      n_b_only <- length(comp$unique_b)
      n_shared <- length(comp$shared)

      # Create bar chart as simple alternative to Venn
      plot_ly() %>%
        add_trace(
          type = "bar",
          x = c(paste0(comp$group_a_name, "\nOnly"),
                "Shared",
                paste0(comp$group_b_name, "\nOnly")),
          y = c(n_a_only, n_shared, n_b_only),
          marker = list(
            color = c("#e74c3c", "#9b59b6", "#3498db")
          ),
          text = c(n_a_only, n_shared, n_b_only),
          textposition = "outside",
          hovertemplate = paste0(
            "<b>%{x}</b><br>",
            "Taxa: %{y}<br>",
            "<extra></extra>"
          )
        ) %>%
        layout(
          title = list(
            text = "Network Taxa Distribution",
            font = list(size = 14)
          ),
          xaxis = list(title = ""),
          yaxis = list(title = "Number of Taxa"),
          showlegend = FALSE,
          bargap = 0.3
        )
    })

    # Venn Summary Statistics
    output$venn_summary <- renderUI({
      if (!analysis_completed_reactive()) {
        return(div(style = "color: gray;", "Run network analysis first"))
      }

      comp <- taxa_comparison()
      if (is.null(comp)) {
        return(div(style = "color: gray;", "Select two different groups"))
      }

      n_a_only <- length(comp$unique_a)
      n_b_only <- length(comp$unique_b)
      n_shared <- length(comp$shared)
      n_total <- n_a_only + n_b_only + n_shared

      # Calculate percentages
      pct_a <- round(n_a_only / n_total * 100, 1)
      pct_b <- round(n_b_only / n_total * 100, 1)
      pct_shared <- round(n_shared / n_total * 100, 1)

      # Jaccard similarity
      jaccard <- round(n_shared / n_total, 3)

      div(
        style = "padding: 15px;",
        h5("Summary"),
        tags$table(
          style = "width: 100%; border-collapse: collapse;",
          tags$tr(
            tags$td(style = "padding: 8px; border-bottom: 1px solid #ddd;",
                   tags$span(style = "color: #e74c3c;", icon("circle")),
                   paste0(" ", comp$group_a_name, " only:")),
            tags$td(style = "padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;",
                   paste0(n_a_only, " (", pct_a, "%)"))
          ),
          tags$tr(
            tags$td(style = "padding: 8px; border-bottom: 1px solid #ddd;",
                   tags$span(style = "color: #9b59b6;", icon("circle")),
                   " Shared:"),
            tags$td(style = "padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;",
                   paste0(n_shared, " (", pct_shared, "%)"))
          ),
          tags$tr(
            tags$td(style = "padding: 8px; border-bottom: 1px solid #ddd;",
                   tags$span(style = "color: #3498db;", icon("circle")),
                   paste0(" ", comp$group_b_name, " only:")),
            tags$td(style = "padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;",
                   paste0(n_b_only, " (", pct_b, "%)"))
          ),
          tags$tr(
            tags$td(style = "padding: 8px; border-bottom: 2px solid #333;", "Total unique taxa:"),
            tags$td(style = "padding: 8px; border-bottom: 2px solid #333; text-align: right; font-weight: bold;",
                   n_total)
          )
        ),
        br(),
        div(
          style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px; margin-top: 10px;",
          tags$b("Jaccard Similarity Index: "),
          tags$span(style = "font-size: 18px; color: #2c3e50;", jaccard),
          br(),
          tags$small(style = "color: #7f8c8d;",
                    "0 = no overlap, 1 = identical networks")
        )
      )
    })

    # Group-Specific Taxa Table
    output$unique_taxa_table <- renderDT({
      if (!analysis_completed_reactive()) {
        return(datatable(
          data.frame(Message = "Run network analysis first"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      comp <- taxa_comparison()
      if (is.null(comp)) {
        return(datatable(
          data.frame(Message = "Select two different groups to compare"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }

      results <- network_results_reactive()
      ps <- filtered_ps_reactive()
      tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
      selected_rank <- tax_rank_reactive()

      # Build table with all taxa and their group membership
      all_taxa <- unique(c(comp$unique_a, comp$unique_b, comp$shared))

      # Helper function to safely get taxonomy column
      get_tax_col <- function(rank, taxa) {
        if (rank %in% colnames(tax_df)) {
          vals <- tax_df[[rank]][match(taxa, rownames(tax_df))]
          vals[is.na(vals)] <- "Unknown"
          return(vals)
        }
        return(rep("Unknown", length(taxa)))
      }

      # Get metrics for each taxon from a specific group
      get_metrics <- function(taxon, group_name) {
        group_data <- results[[group_name]]
        if (is.null(group_data)) return(list(degree = NA, ivi = NA, role = NA))

        node_idx <- which(group_data$node_metrics$node == taxon)
        ivi_idx <- which(group_data$ivi$node == taxon)
        zipi_idx <- which(group_data$zi_pi$node == taxon)

        list(
          degree = if (length(node_idx) > 0) group_data$node_metrics$degree[node_idx] else NA,
          ivi = if (length(ivi_idx) > 0) group_data$ivi$IVI[ivi_idx] else NA,
          role = if (length(zipi_idx) > 0) group_data$zi_pi$role[zipi_idx] else NA
        )
      }

      # Create column names with group names
      grp_a <- comp$group_a_name
      grp_b <- comp$group_b_name

      # Start building table
      table_data <- data.frame(
        Taxon = all_taxa,
        Phylum = get_tax_col("Phylum", all_taxa),
        stringsAsFactors = FALSE
      )

      # Add selected rank if different from Phylum
      if (selected_rank != "Phylum") {
        table_data[[selected_rank]] <- get_tax_col(selected_rank, all_taxa)
      }

      # Determine group membership
      table_data$Category <- ifelse(
        all_taxa %in% comp$shared, "Shared",
        ifelse(all_taxa %in% comp$unique_a,
               paste0(grp_a, " only"),
               paste0(grp_b, " only"))
      )

      # Get metrics from BOTH groups for all taxa
      table_data[[paste0("Degree_", grp_a)]] <- NA
      table_data[[paste0("IVI_", grp_a)]] <- NA
      table_data[[paste0("Role_", grp_a)]] <- NA
      table_data[[paste0("Degree_", grp_b)]] <- NA
      table_data[[paste0("IVI_", grp_b)]] <- NA
      table_data[[paste0("Role_", grp_b)]] <- NA

      for (i in 1:nrow(table_data)) {
        taxon <- table_data$Taxon[i]

        # Get metrics from Group A if taxon exists there
        if (taxon %in% comp$all_a) {
          metrics_a <- get_metrics(taxon, grp_a)
          table_data[[paste0("Degree_", grp_a)]][i] <- metrics_a$degree
          table_data[[paste0("IVI_", grp_a)]][i] <- metrics_a$ivi
          table_data[[paste0("Role_", grp_a)]][i] <- metrics_a$role
        }

        # Get metrics from Group B if taxon exists there
        if (taxon %in% comp$all_b) {
          metrics_b <- get_metrics(taxon, grp_b)
          table_data[[paste0("Degree_", grp_b)]][i] <- metrics_b$degree
          table_data[[paste0("IVI_", grp_b)]][i] <- metrics_b$ivi
          table_data[[paste0("Role_", grp_b)]][i] <- metrics_b$role
        }
      }

      # Calculate max IVI for sorting (use max of both groups)
      ivi_a_col <- paste0("IVI_", grp_a)
      ivi_b_col <- paste0("IVI_", grp_b)
      table_data$max_IVI <- pmax(
        ifelse(is.na(table_data[[ivi_a_col]]), 0, table_data[[ivi_a_col]]),
        ifelse(is.na(table_data[[ivi_b_col]]), 0, table_data[[ivi_b_col]])
      )

      # Sort: Group-specific first (by max IVI), then shared
      table_data <- table_data %>%
        arrange(
          factor(Category, levels = c(
            paste0(grp_a, " only"),
            paste0(grp_b, " only"),
            "Shared"
          )),
          desc(max_IVI)
        ) %>%
        select(-max_IVI)  # Remove sorting column from display

      # Column indices for formatting (0-based)
      ivi_cols <- c(paste0("IVI_", grp_a), paste0("IVI_", grp_b))

      datatable(
        table_data,
        extensions = 'Buttons',
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = list(
            list(extend = 'csv', filename = paste0('group_specific_taxa_', grp_a, '_vs_', grp_b)),
            list(extend = 'excel', filename = paste0('group_specific_taxa_', grp_a, '_vs_', grp_b))
          )
        ),
        rownames = FALSE,
        filter = "top"
      ) %>%
        formatRound(columns = ivi_cols, digits = 3) %>%
        formatStyle(
          "Category",
          backgroundColor = styleEqual(
            c(paste0(grp_a, " only"),
              paste0(grp_b, " only"),
              "Shared"),
            c("#ffebee", "#e3f2fd", "#f3e5f5")  # Light red, light blue, light purple
          )
        )
    })

    # =========================================================================
    # Download Handlers with Full Preview
    # =========================================================================

    # Zi-Pi plot base function
    zipi_base_plot <- function() {
      results <- network_results_reactive()
      grp <- input$group_select
      tax_rank <- input$tax_rank
      if (is.null(results) || is.null(grp) || is.null(results[[grp]])) return(NULL)

      zipi_df <- results[[grp]]$zi_pi
      ivi_df <- results[[grp]]$ivi
      g <- results[[grp]]$igraph

      if (is.null(zipi_df)) return(NULL)

      # Add IVI if available
      if (!is.null(ivi_df)) {
        zipi_df$IVI <- ivi_df$IVI[match(zipi_df$node, ivi_df$node)]
      } else {
        zipi_df$IVI <- 1
      }

      # Add taxonomy coloring (matching plotly)
      if (!is.null(g) && !is.null(igraph::V(g)$taxonomy)) {
        taxonomy_vals <- igraph::V(g)$taxonomy
        node_names <- igraph::V(g)$name
        zipi_df$Taxonomy <- taxonomy_vals[match(zipi_df$node, node_names)]

        # Extract the selected taxonomic rank from full taxonomy string
        if (!is.null(tax_rank) && tax_rank != "Species") {
          rank_prefixes <- c("Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
                            "Order" = "o__", "Family" = "f__", "Genus" = "g__", "Species" = "s__")
          prefix <- rank_prefixes[[tax_rank]]
          if (!is.null(prefix)) {
            zipi_df$Taxonomy <- sapply(zipi_df$Taxonomy, function(tax) {
              if (is.na(tax)) return("Unknown")
              parts <- strsplit(tax, ";")[[1]]
              matched <- grep(paste0("^", prefix), parts, value = TRUE)
              if (length(matched) > 0) {
                gsub(paste0("^", prefix), "", matched[1])
              } else {
                "Unknown"
              }
            })
          }
        }

        # Limit to top taxa for cleaner legend
        taxa_counts <- table(zipi_df$Taxonomy)
        top_taxa <- names(sort(taxa_counts, decreasing = TRUE))[1:min(15, length(taxa_counts))]
        zipi_df$Taxonomy <- ifelse(zipi_df$Taxonomy %in% top_taxa, zipi_df$Taxonomy, "Other")

      } else {
        zipi_df$Taxonomy <- "Unknown"
      }

      # Get Pi and Zi ranges for annotation positioning
      pi_range <- range(zipi_df$Pi, na.rm = TRUE)
      zi_range <- range(zipi_df$Zi, na.rm = TRUE)

      p <- ggplot(zipi_df, aes(x = Pi, y = Zi, size = IVI, color = Taxonomy)) +
        geom_point(alpha = 0.7) +
        geom_hline(yintercept = 2.5, linetype = "dashed", color = "#d62728", linewidth = 0.8) +
        geom_vline(xintercept = 0.62, linetype = "dashed", color = "#d62728", linewidth = 0.8) +
        scale_size_continuous(range = c(2, 10), name = "IVI Score") +
        labs(
          title = paste0("Zi-Pi Plot: ", grp),
          x = "Participation Coefficient (Pi)",
          y = "Within-module Degree (Zi)",
          color = tax_rank
        ) +
        # Add quadrant labels matching plotly style
        annotate("text", x = 0.31, y = max(3.5, zi_range[2] * 0.9),
                 label = "Module Hubs", size = 3.5, color = "#666666", fontface = "bold") +
        annotate("text", x = 0.81, y = max(3.5, zi_range[2] * 0.9),
                 label = "Network Hubs", size = 3.5, color = "#666666", fontface = "bold") +
        annotate("text", x = 0.31, y = min(0.5, zi_range[1] + 0.5),
                 label = "Peripherals", size = 3.5, color = "#666666", fontface = "bold") +
        annotate("text", x = 0.81, y = min(0.5, zi_range[1] + 0.5),
                 label = "Connectors", size = 3.5, color = "#666666", fontface = "bold") +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

      return(p)
    }

    # Setup Zi-Pi download with full preview
    setup_full_download_handler(
      input, output, session,
      plot_id = "zipi",
      plot_fn = zipi_base_plot,
      filename_prefix = paste0("zipi_plot_", isolate(input$group_select)),
      modal_title = "Download Zi-Pi Plot"
    )

    # Differential Zi-Pi plot base function
    diff_zipi_base_plot <- function() {
      results <- network_results_reactive()
      grp_a <- input$diff_group_a
      grp_b <- input$diff_group_b
      tax_rank <- input$diff_tax_rank

      if (is.null(results) || is.null(grp_a) || is.null(grp_b)) return(NULL)

      zipi_a <- results[[grp_a]]$zi_pi
      zipi_b <- results[[grp_b]]$zi_pi
      g_a <- results[[grp_a]]$igraph

      if (is.null(zipi_a) || is.null(zipi_b)) return(NULL)

      # Calculate differences
      common_taxa <- intersect(zipi_a$node, zipi_b$node)
      if (length(common_taxa) == 0) return(NULL)

      diff_data <- data.frame(
        node = common_taxa,
        delta_Zi = zipi_a$Zi[match(common_taxa, zipi_a$node)] -
                   zipi_b$Zi[match(common_taxa, zipi_b$node)],
        delta_Pi = zipi_a$Pi[match(common_taxa, zipi_a$node)] -
                   zipi_b$Pi[match(common_taxa, zipi_b$node)]
      )

      # Add taxonomy coloring
      if (!is.null(g_a) && !is.null(igraph::V(g_a)$taxonomy)) {
        taxonomy_vals <- igraph::V(g_a)$taxonomy
        node_names <- igraph::V(g_a)$name
        diff_data$Taxonomy <- taxonomy_vals[match(diff_data$node, node_names)]

        # Extract the selected taxonomic rank
        if (!is.null(tax_rank) && tax_rank != "Species") {
          rank_prefixes <- c("Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
                            "Order" = "o__", "Family" = "f__", "Genus" = "g__", "Species" = "s__")
          prefix <- rank_prefixes[[tax_rank]]
          if (!is.null(prefix)) {
            diff_data$Taxonomy <- sapply(diff_data$Taxonomy, function(tax) {
              if (is.na(tax)) return("Unknown")
              parts <- strsplit(tax, ";")[[1]]
              matched <- grep(paste0("^", prefix), parts, value = TRUE)
              if (length(matched) > 0) {
                gsub(paste0("^", prefix), "", matched[1])
              } else {
                "Unknown"
              }
            })
          }
        }

        # Limit to top taxa
        taxa_counts <- table(diff_data$Taxonomy)
        top_taxa <- names(sort(taxa_counts, decreasing = TRUE))[1:min(15, length(taxa_counts))]
        diff_data$Taxonomy <- ifelse(diff_data$Taxonomy %in% top_taxa, diff_data$Taxonomy, "Other")
      } else {
        diff_data$Taxonomy <- "Unknown"
      }

      p <- ggplot(diff_data, aes(x = delta_Pi, y = delta_Zi, color = Taxonomy)) +
        geom_point(alpha = 0.7, size = 3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#7f7f7f", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#7f7f7f", linewidth = 0.8) +
        labs(
          title = paste0("Differential Zi-Pi: ", grp_a, " vs ", grp_b),
          subtitle = paste0("Positive = higher in ", grp_a, ", Negative = higher in ", grp_b),
          x = expression(Delta ~ "Pi (Participation Coefficient)"),
          y = expression(Delta ~ "Zi (Within-module Degree)"),
          color = tax_rank
        ) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

      return(p)
    }

    # Setup Differential Zi-Pi download with full preview
    setup_full_download_handler(
      input, output, session,
      plot_id = "diff_zipi",
      plot_fn = diff_zipi_base_plot,
      filename_prefix = "diff_zipi_plot",
      modal_title = "Download Differential Zi-Pi Plot"
    )

  })
}
