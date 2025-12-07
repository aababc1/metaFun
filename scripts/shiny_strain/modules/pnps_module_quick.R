# =============================================================================
# pN/pS Analysis Module - Quick Functional Version
# =============================================================================

library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)

# Load statistical helpers
if(!exists("apply_bh_correction")) {
  source("modules/nucleotide_diversity_stats_helpers.R", local = TRUE)
}

# ==============================================================================
# UI Function
# ==============================================================================

pnps_UI <- function(id) {
  ns <- NS(id)

  tagList(
    card(
      full_screen = TRUE,
      card_header(
        style = "background-color: #97C6F1; color: white;",
        "pN/pS Analysis"
      ),
      card_body(
        layout_sidebar(
          sidebar = sidebar(
            title = "Analysis Settings",

            selectInput(ns("tax_level"),
                       "Taxonomic Level:",
                       choices = c("Phylum" = "phylum",
                                 "Class" = "class",
                                 "Order" = "order",
                                 "Family" = "family",
                                 "Genus" = "genus",
                                 "Species" = "species"),
                       selected = "genus"),

            selectInput(ns("metadata_var"),
                       "Metadata Variable:",
                       choices = NULL),

            # Instructional text
            div(
              style = "background-color: #e7f3ff; padding: 8px; margin: 10px 0; border-radius: 4px; border-left: 3px solid #2196F3;",
              tags$small(
                icon("info-circle"),
                " After selecting conditions, click 'Run Analysis' to update the plot and statistics."
              )
            ),

            radioButtons(ns("transformation"),
                        "Transformation:",
                        choices = c("Log10" = "log", "Raw" = "raw"),
                        selected = "raw"),

            sliderInput(ns("min_samples"),
                       "Min Samples per Taxon:",
                       min = 1, max = 20, value = 3, step = 1),

            checkboxInput(ns("show_significant_only"),
                         "Show only significant taxa",
                         value = TRUE),

            hr(),

            # Quality Filters (unified with Diversity vs pN/pS)
            tags$b("Quality Filters"),
            sliderInput(ns("min_coverage"),
                       "Minimum Coverage:",
                       min = 1, max = 50, value = 5, step = 1),

            sliderInput(ns("min_breadth"),
                       "Minimum Breadth:",
                       min = 0.1, max = 1.0, value = 0.5, step = 0.1),

            # Filter status display
            div(
              id = ns("filter_status"),
              style = "background-color: #f8f9fa; padding: 8px; margin: 10px 0; border-radius: 4px; border: 1px solid #dee2e6; font-size: 12px;",
              tags$b(icon("filter"), " Applied Filters:"),
              br(),
              uiOutput(ns("filter_info"))
            ),

            hr(),

            actionButton(ns("run_analysis"),
                        "Run Analysis",
                        class = "btn-primary btn-block"),

            br(),

            actionButton(ns("regenerate_pnps"),
                        "Regenerate pN/pS Data",
                        class = "btn-warning btn-block",
                        icon = icon("refresh")) %>%
              tooltip("Delete cached pN/pS data and recalculate from raw files")
          ),

          # Main content
          div(
            div(
              style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
              h4("Group Comparison", style = "margin: 0;"),
              create_plot_download_btn(ns, "comparison")
            ),
            plotlyOutput(ns("comparison_plot"), height = "600px"),
            br(),

            h4("Statistical Results"),
            DT::dataTableOutput(ns("stats_table")),
            br(),

            h4("Data Summary"),
            verbatimTextOutput(ns("summary_text"))
          )
        )
      )
    )
  )
}

# ==============================================================================
# Server Function
# ==============================================================================

pnps_Server <- function(id, pnps_data, diversity_data = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values
    rv <- reactiveValues(
      genome_wide = NULL,
      genome_wide_raw = NULL,  # Before quality filters
      filtered_data = NULL,
      stats_results = NULL,
      metadata_vars = NULL,
      filter_stats = list(before = 0, after = 0)
    )

    # Load data
    observe({
      req(pnps_data())

      data <- pnps_data()

      if(!is.null(data$genome_wide_pnps)) {
        rv$genome_wide_raw <- data$genome_wide_pnps
        cat("Loaded genome-wide pN/pS:", nrow(data$genome_wide_pnps), "records\n")

        # Detect metadata variables
        exclude_cols <- c("sample_id", "clean_genome", "genome", "total_SNV_N", "total_SNV_S",
                         "total_N_sites", "total_S_sites", "n_genes", "pN", "pS", "pNpS",
                         "domain", "phylum", "class", "order", "family", "genus", "species",
                         "coverage", "breadth", "nucl_diversity")

        meta_cols <- setdiff(names(data$genome_wide_pnps), exclude_cols)
        rv$metadata_vars <- meta_cols

        # Update metadata variable choices
        updateSelectInput(session, "metadata_var",
                         choices = c("None" = "none", meta_cols),
                         selected = "none")
      }
    })

    # Apply quality filters (unified with Diversity vs pN/pS module)
    observe({
      req(rv$genome_wide_raw)

      pnps_raw <- rv$genome_wide_raw
      rv$filter_stats$before <- nrow(pnps_raw)

      # Check if diversity_data is available for quality filtering
      div_data <- tryCatch(diversity_data(), error = function(e) NULL)

      if(!is.null(div_data) && nrow(div_data) > 0) {
        # Join with diversity_data to get coverage/breadth
        div_select <- div_data %>%
          select(sample_id, clean_genome, coverage, breadth) %>%
          distinct()

        # Inner join to apply quality filters
        filtered_pnps <- pnps_raw %>%
          inner_join(div_select, by = c("sample_id", "clean_genome")) %>%
          filter(coverage >= input$min_coverage,
                 breadth >= input$min_breadth)

        rv$genome_wide <- filtered_pnps
        rv$filter_stats$after <- nrow(filtered_pnps)

        cat("Quality filter applied: coverage >=", input$min_coverage,
            ", breadth >=", input$min_breadth, "\n")
        cat("Records:", rv$filter_stats$before, "->", rv$filter_stats$after, "\n")
      } else {
        # No diversity_data available, use raw pnps data
        rv$genome_wide <- pnps_raw
        rv$filter_stats$after <- nrow(pnps_raw)
        cat("No diversity_data available for quality filtering, using all pN/pS records\n")
      }
    })

    # Filter info display
    output$filter_info <- renderUI({
      div(
        tags$span(paste0("Coverage >= ", input$min_coverage)),
        br(),
        tags$span(paste0("Breadth >= ", input$min_breadth)),
        br(),
        tags$span(
          style = "color: #6c757d;",
          paste0("Records: ", rv$filter_stats$before, " -> ", rv$filter_stats$after,
                 " (", round(rv$filter_stats$after / max(rv$filter_stats$before, 1) * 100, 1), "%)")
        )
      )
    })

    # Regenerate pN/pS data
    observeEvent(input$regenerate_pnps, {
      showModal(modalDialog(
        title = "Regenerate pN/pS Data",
        "This will delete the cached pN/pS data and recalculate from raw InStrain files. This may take several minutes. Continue?",
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_regenerate"), "Yes, Regenerate", class = "btn-danger")
        )
      ))
    })

    # Confirm regenerate
    observeEvent(input$confirm_regenerate, {
      removeModal()

      pnps_file <- "data/pN_pS_integrated_data.rds"

      if(file.exists(pnps_file)) {
        file.remove(pnps_file)
        showNotification("Cached pN/pS data deleted. Please click 'Load All Data' to regenerate.",
                        type = "warning", duration = 10)
      } else {
        showNotification("No cached pN/pS data found.", type = "message")
      }

      # Clear current data
      rv$genome_wide <- NULL
      rv$filtered_data <- NULL
      rv$stats_results <- NULL
    })

    # Run analysis
    observeEvent(input$run_analysis, {
      req(rv$genome_wide, input$tax_level, input$metadata_var)

      taxonomy_level <- input$tax_level
      metadata_var <- input$metadata_var
      transformation <- input$transformation
      min_samples <- input$min_samples

      cat("\n=== Running pN/pS Analysis ===\n")
      cat("Taxonomy level:", taxonomy_level, "\n")
      cat("Metadata var:", metadata_var, "\n")

      # Filter data
      filtered <- rv$genome_wide %>%
        filter(!is.na(.data[[taxonomy_level]]),
               !is.na(pNpS),
               pNpS > 0)

      # Apply transformation
      filtered$pNpS_transformed <- if(transformation == "log") {
        log10(filtered$pNpS + 1e-10)
      } else {
        filtered$pNpS
      }

      # Filter by min samples per taxon
      taxon_counts <- filtered %>%
        group_by(.data[[taxonomy_level]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n >= min_samples)

      filtered <- filtered %>%
        filter(.data[[taxonomy_level]] %in% taxon_counts[[taxonomy_level]])

      cat("Filtered data:", nrow(filtered), "records,", nrow(taxon_counts), "taxa\n")

      rv$filtered_data <- filtered

      # Run statistical tests
      if(metadata_var != "none") {
        cat("Running statistical tests...\n")

        # Filter out NA metadata values
        test_data <- filtered %>%
          filter(!is.na(.data[[metadata_var]]))

        if(nrow(test_data) > 0) {
          rv$stats_results <- run_pnps_stats(
            test_data,
            taxonomy_level,
            metadata_var
          )

          if(!is.null(rv$stats_results)) {
            cat("Stats completed:", nrow(rv$stats_results), "taxa tested\n")
          }
        }
      } else {
        rv$stats_results <- NULL
      }
    })

    # Group comparison plot
    output$comparison_plot <- renderPlotly({
      req(rv$filtered_data, input$tax_level)

      data <- rv$filtered_data
      taxonomy_level <- input$tax_level
      metadata_var <- input$metadata_var
      transformation <- input$transformation
      stats_results <- rv$stats_results
      show_significant_only <- input$show_significant_only

      y_label <- if(transformation == "log") "pN/pS (log10)" else "pN/pS"

      # Filter to significant taxa only if requested
      if(show_significant_only && !is.null(stats_results) && nrow(stats_results) > 0) {
        significant_taxa <- stats_results %>%
          filter(significance != "ns") %>%
          pull(taxon)

        if(length(significant_taxa) > 0) {
          data <- data %>% filter(.data[[taxonomy_level]] %in% significant_taxa)
        } else {
          # No significant taxa, show message
          return(plotly::plot_ly() %>%
                  plotly::add_text(x = 0.5, y = 0.5, text = "No significant taxa found",
                                  textfont = list(size = 20)) %>%
                  plotly::layout(xaxis = list(visible = FALSE),
                                yaxis = list(visible = FALSE)))
        }
      }

      # Order taxa by significance if stats available
      if(!is.null(stats_results) && nrow(stats_results) > 0) {
        taxa_order <- stats_results %>%
          arrange(p.value) %>%
          pull(taxon)

        data[[taxonomy_level]] <- factor(data[[taxonomy_level]], levels = taxa_order)
      }

      # Determine if metadata variable is numerical
      is_numerical <- suppressWarnings({
        var_data <- data[[metadata_var]]
        var_data <- var_data[!is.na(var_data)]
        if(length(var_data) == 0) return(FALSE)
        numeric_data <- as.numeric(as.character(var_data))
        !all(is.na(numeric_data))
      })

      # Create plot based on variable type
      if(metadata_var != "none" && metadata_var %in% names(data)) {

        if(is_numerical) {
          # Scatter plot with trend lines for numerical variables
          data$meta_numeric <- as.numeric(as.character(data[[metadata_var]]))

          # Determine number of unique taxa for color palette selection
          n_taxa <- length(unique(data[[taxonomy_level]]))

          p <- ggplot(data, aes(x = meta_numeric, y = pNpS_transformed,
                               color = .data[[taxonomy_level]], fill = .data[[taxonomy_level]])) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
            labs(title = paste("pN/pS by", taxonomy_level, "and", metadata_var),
                 x = metadata_var,
                 y = y_label,
                 color = taxonomy_level,
                 fill = taxonomy_level) +
            theme_minimal() +
            theme(legend.position = "right",
                  legend.title = element_text(size = 10, face = "bold"),
                  legend.text = element_text(size = 9))

          # Choose appropriate color palette based on number of taxa
          if(n_taxa <= 9) {
            p <- p +
              scale_color_brewer(palette = "Set1") +
              scale_fill_brewer(palette = "Set1")
          } else if(n_taxa <= 12) {
            p <- p +
              scale_color_brewer(palette = "Set3") +
              scale_fill_brewer(palette = "Set3")
          } else {
            # Use viridis for many taxa (unlimited colors)
            p <- p +
              scale_color_viridis_d(option = "turbo") +
              scale_fill_viridis_d(option = "turbo")
          }

        } else {
          # Grouped boxplot for categorical variables - match Nucleotide Diversity style
          # 1) Convert metadata variable to factor explicitly and drop unused levels
          data[[metadata_var]] <- as.factor(data[[metadata_var]])
          data[[metadata_var]] <- forcats::fct_drop(data[[metadata_var]])

          # Determine number of groups for color palette selection
          n_groups <- length(unique(data[[metadata_var]]))

          # 2) Define dodge parameter for reuse (optimized width for side-by-side display)
          pd <- position_dodge2(width = 0.65, preserve = "single")

          # 3) Create boxplot with explicit group aesthetic (critical for ggplotly)
          # Use fill only (not color) to avoid duplicate legends
          p <- ggplot(data, aes_string(
            x = taxonomy_level,
            y = "pNpS_transformed",
            fill = metadata_var,
            group = metadata_var  # Explicit grouping for proper dodge behavior
          )) +
            geom_boxplot(
              outlier.shape = NA,
              position = pd,
              width = 0.6
            ) +
            geom_point(
              position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65),
              alpha = 0.35,
              size = 0.9,
              show.legend = FALSE  # Hide points from legend to avoid duplication
            ) +
            labs(
              title = paste("pN/pS by", taxonomy_level, "and", metadata_var),
              x = taxonomy_level,
              y = y_label,
              fill = metadata_var
            ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right",
              legend.text = element_text(size = 9)
            )

          # Choose appropriate color palette based on number of groups
          if(n_groups <= 9) {
            p <- p + scale_fill_brewer(palette = "Set1")
          } else if(n_groups <= 12) {
            p <- p + scale_fill_brewer(palette = "Set3")
          } else {
            # Use viridis for many groups (unlimited colors)
            p <- p + scale_fill_viridis_d(option = "turbo")
          }
        }

        # Add significance annotations (only for categorical/boxplot)
        if(!is_numerical && !is.null(stats_results) && nrow(stats_results) > 0) {
          max_y <- max(data$pNpS_transformed, na.rm = TRUE)

          annotations <- stats_results %>%
            select(taxon, significance) %>%
            filter(!is.na(significance), significance != "ns")

          if(nrow(annotations) > 0) {
            annotations[[taxonomy_level]] <- annotations$taxon

            p <- p + geom_text(data = annotations,
                              aes_string(x = taxonomy_level, label = "significance"),
                              y = max_y * 1.05,
                              size = 5,
                              fontface = "bold",
                              inherit.aes = FALSE)
          }
        }

        # For numerical variables, add correlation annotations in legend area
        if(is_numerical && !is.null(stats_results) && nrow(stats_results) > 0) {
          # Add correlation info as subtitle
          n_sig <- sum(stats_results$significance != "ns", na.rm = TRUE)
          subtitle_text <- paste0(nrow(stats_results), " taxa tested, ",
                                 n_sig, " significant correlations")

          p <- p + labs(subtitle = subtitle_text)
        }

      } else {
        # Simple boxplot (no metadata grouping)
        p <- ggplot(data, aes(x = .data[[taxonomy_level]], y = pNpS_transformed)) +
          geom_boxplot(aes(fill = .data[[taxonomy_level]]), alpha = 0.7, outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
          labs(title = paste("pN/pS by", taxonomy_level),
               x = taxonomy_level,
               y = y_label) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none")
      }

      # Convert to plotly with boxmode="group" for proper side-by-side display
      # This is CRITICAL to prevent overlapping boxes when multiple groups exist
      ggplotly(p, tooltip = c("x", "y", "fill")) %>%
        plotly::layout(boxmode = "group")
    })

    # Statistical results table
    output$stats_table <- DT::renderDataTable({
      req(rv$stats_results)

      stats <- rv$stats_results
      show_significant_only <- input$show_significant_only

      if(nrow(stats) == 0) {
        return(NULL)
      }

      # Filter to significant taxa only if requested
      if(show_significant_only) {
        stats <- stats %>% filter(significance != "ns")

        if(nrow(stats) == 0) {
          return(NULL)
        }
      }

      # Format columns based on test type
      display_cols <- c("taxon", "test_type", "n", "statistic", "p.value", "p.adj", "significance")

      if("rho" %in% names(stats)) {
        display_cols <- c(display_cols[1:3], "rho", display_cols[4:length(display_cols)])
      }

      display_data <- stats[, intersect(display_cols, names(stats))]

      # Rename columns
      col_names <- c(
        taxon = "Taxon",
        test_type = "Test",
        n = "N",
        rho = "Correlation",
        statistic = "Statistic",
        p.value = "P-value",
        p.adj = "P-adj (BH)",
        significance = "Sig"
      )

      for(old in names(col_names)) {
        if(old %in% names(display_data)) {
          names(display_data)[names(display_data) == old] <- col_names[old]
        }
      }

      # Create datatable
      dt <- DT::datatable(
        display_data,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          order = list(list(which(names(display_data) == "P-value") - 1, 'asc'))
        ),
        rownames = FALSE,
        extensions = 'Buttons'
      )

      # Format Correlation (rho) column for Spearman results
      if("Correlation" %in% names(display_data)) {
        dt <- dt %>% DT::formatRound(columns = "Correlation", digits = 4)
      }

      # Format Statistic column
      if("Statistic" %in% names(display_data)) {
        dt <- dt %>% DT::formatRound(columns = "Statistic", digits = 2)
      }

      # Format P-value column
      if("P-value" %in% names(display_data)) {
        dt <- dt %>% DT::formatRound(columns = "P-value", digits = 4)
      }

      # Format P-adj column only if it contains non-NA values
      if("P-adj (BH)" %in% names(display_data) && any(!is.na(display_data[["P-adj (BH)"]]))) {
        dt <- dt %>% DT::formatRound(columns = "P-adj (BH)", digits = 4)
      }

      # Color-code significance
      if("Sig" %in% names(display_data)) {
        dt <- dt %>% DT::formatStyle(
          "Sig",
          backgroundColor = DT::styleEqual(
            c("***", "**", "*", "ns"),
            c("#d4edda", "#fff3cd", "#ffe5d0", "#f8f9fa")
          )
        )
      }

      dt
    })

    # Summary text
    output$summary_text <- renderPrint({
      req(rv$filtered_data)

      data <- rv$filtered_data

      cat("pN/pS Data Summary\n")
      cat("==================\n\n")
      cat("Total records:", nrow(data), "\n")
      cat("Unique samples:", length(unique(data$sample_id)), "\n")
      cat("Unique genomes:", length(unique(data$clean_genome)), "\n")
      cat("Taxa analyzed:", length(unique(data[[input$tax_level]])), "\n\n")
      cat("pN/pS statistics:\n")
      cat("  Min:", round(min(data$pNpS, na.rm = TRUE), 4), "\n")
      cat("  Median:", round(median(data$pNpS, na.rm = TRUE), 4), "\n")
      cat("  Mean:", round(mean(data$pNpS, na.rm = TRUE), 4), "\n")
      cat("  Max:", round(max(data$pNpS, na.rm = TRUE), 4), "\n")
    })

    # ========== Download Modal Handler ==========
    comparison_ggplot <- reactive({
      req(rv$filtered_data, input$tax_level)

      data <- rv$filtered_data
      taxonomy_level <- input$tax_level
      metadata_var <- input$metadata_var
      transformation <- input$transformation
      stats_results <- rv$stats_results
      show_significant_only <- input$show_significant_only

      y_label <- if(transformation == "log") "pN/pS (log10)" else "pN/pS"

      # Filter to significant taxa only if requested
      if(show_significant_only && !is.null(stats_results) && nrow(stats_results) > 0) {
        significant_taxa <- stats_results %>%
          filter(significance != "ns") %>%
          pull(taxon)

        if(length(significant_taxa) > 0) {
          data <- data %>% filter(.data[[taxonomy_level]] %in% significant_taxa)
        }
      }

      # Order taxa by significance if stats available
      if(!is.null(stats_results) && nrow(stats_results) > 0) {
        taxa_order <- stats_results %>%
          arrange(p.value) %>%
          pull(taxon)
        data[[taxonomy_level]] <- factor(data[[taxonomy_level]], levels = taxa_order)
      }

      req(nrow(data) > 0)

      if(metadata_var != "none" && metadata_var %in% names(data)) {
        # Categorical: boxplot
        data[[metadata_var]] <- as.factor(data[[metadata_var]])
        data[[metadata_var]] <- forcats::fct_drop(data[[metadata_var]])

        p <- ggplot(data, aes_string(x = taxonomy_level, y = "pNpS_transformed", fill = metadata_var)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge2(width = 0.65)) +
          geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65),
                    alpha = 0.35, size = 0.9, show.legend = FALSE) +
          labs(title = paste("pN/pS by", taxonomy_level, "and", metadata_var),
               x = taxonomy_level, y = y_label, fill = metadata_var) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
          scale_fill_brewer(palette = "Set1")
      } else {
        p <- ggplot(data, aes_string(x = taxonomy_level, y = "pNpS_transformed", fill = taxonomy_level)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
          labs(title = paste("pN/pS by", taxonomy_level), x = taxonomy_level, y = y_label) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      }

      p
    })

    setup_download_modal_handler(input, output, session, "comparison", comparison_ggplot, "pnps_comparison")

  })
}

# ==============================================================================
# Helper Functions
# ==============================================================================

run_pnps_stats <- function(data, taxonomy_level, metadata_var) {

  # Check metadata variable type (robust detection like Nucleotide Diversity module)
  is_numerical <- suppressWarnings({
    var_data <- data[[metadata_var]]
    var_data <- var_data[!is.na(var_data)]
    if(length(var_data) == 0) return(FALSE)

    # Try to convert to numeric
    numeric_data <- as.numeric(as.character(var_data))
    # If conversion successful (not all NA), it's numerical
    !all(is.na(numeric_data))
  })

  cat("  Variable type:", ifelse(is_numerical, "numerical (Spearman)", "categorical"), "\n")

  # Get unique taxa
  taxa <- unique(data[[taxonomy_level]])

  results_list <- list()

  for(taxon in taxa) {
    taxon_data <- data %>%
      filter(.data[[taxonomy_level]] == taxon,
             !is.na(pNpS_transformed),
             !is.na(.data[[metadata_var]]))

    if(nrow(taxon_data) < 3) next

    if(is_numerical) {
      # Spearman correlation for numerical variables
      meta_values <- as.numeric(as.character(taxon_data[[metadata_var]]))
      pnps_values <- taxon_data$pNpS_transformed

      # Skip if not enough valid numeric values
      valid_pairs <- !is.na(meta_values) & !is.na(pnps_values)
      if(sum(valid_pairs) < 3) next

      cor_result <- tryCatch({
        cor.test(meta_values[valid_pairs], pnps_values[valid_pairs],
                method = "spearman", exact = FALSE)
      }, error = function(e) {
        cat("    ERROR in Spearman for", taxon, ":", e$message, "\n")
        list(estimate = NA, p.value = NA, statistic = NA)
      })

      results_list[[taxon]] <- data.frame(
        taxon = taxon,
        test_type = "Spearman",
        n = nrow(taxon_data),
        rho = cor_result$estimate,
        statistic = cor_result$statistic,
        p.value = cor_result$p.value,
        stringsAsFactors = FALSE
      )

    } else {
      # Categorical variable
      groups <- unique(taxon_data[[metadata_var]])
      n_groups <- length(groups)

      if(n_groups < 2) next

      if(n_groups == 2) {
        # Wilcoxon test - use formula with proper syntax
        formula_str <- as.formula(paste("pNpS_transformed ~", metadata_var))
        test_result <- tryCatch({
          wilcox.test(formula_str, data = taxon_data)
        }, error = function(e) {
          cat("  ERROR in Wilcoxon for", taxon, ":", e$message, "\n")
          list(statistic = NA, p.value = NA)
        })

        results_list[[taxon]] <- data.frame(
          taxon = taxon,
          test_type = "Wilcoxon",
          n = nrow(taxon_data),
          statistic = as.numeric(test_result$statistic),
          p.value = as.numeric(test_result$p.value),
          stringsAsFactors = FALSE
        )

      } else {
        # Kruskal-Wallis test - use formula with proper syntax
        formula_str <- as.formula(paste("pNpS_transformed ~", metadata_var))
        test_result <- tryCatch({
          kruskal.test(formula_str, data = taxon_data)
        }, error = function(e) {
          cat("  ERROR in Kruskal-Wallis for", taxon, ":", e$message, "\n")
          list(statistic = NA, p.value = NA)
        })

        results_list[[taxon]] <- data.frame(
          taxon = taxon,
          test_type = "Kruskal-Wallis",
          n = nrow(taxon_data),
          statistic = as.numeric(test_result$statistic),
          p.value = as.numeric(test_result$p.value),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if(length(results_list) == 0) {
    return(NULL)
  }

  # Combine results
  results <- bind_rows(results_list)

  # Apply BH correction ONLY for Kruskal-Wallis (multi-level categorical)
  # NO correction for: Wilcoxon (binary), Spearman (numerical)
  if(unique(results$test_type)[1] == "Kruskal-Wallis") {
    # Multiple group comparisons - apply BH correction
    results$p.adj <- p.adjust(results$p.value, method = "BH")
    # Use adjusted p-value for significance
    sig_p <- results$p.adj
  } else {
    # Wilcoxon or Spearman - NO correction applied
    results$p.adj <- NA  # Not applicable
    # Use raw p-value for significance
    sig_p <- results$p.value
  }

  # Add significance notation based on appropriate p-value
  results$significance <- sapply(sig_p, function(p) {
    if(is.na(p)) return("ns")
    if(p < 0.001) return("***")
    if(p < 0.01) return("**")
    if(p < 0.05) return("*")
    return("ns")
  })

  # Order by p-value
  results <- results %>% arrange(p.value)

  return(results)
}
