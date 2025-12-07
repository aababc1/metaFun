# ==============================================================================
# Nucleotide Diversity vs pN/pS Correlation Module
# ==============================================================================
# This module creates scatter plots showing the relationship between
# nucleotide diversity (X-axis) and pN/pS ratio (Y-axis) averaged by
# taxonomy groups. Same taxa across metadata groups are connected by
# dashed lines for easy comparison.
# ==============================================================================

library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(RColorBrewer)

# ==============================================================================
# UI Function
# ==============================================================================

diversity_pnps_correlation_UI <- function(id) {
  ns <- NS(id)

  tagList(
    # Full height container - matching other pN/pS tabs
    card(
      full_screen = TRUE,
      class = "main-card",
      fill = TRUE,
      card_header(
        style = "background-color: #97C6F1; color: white;",
        tags$span("Nucleotide Diversity vs pN/pS Correlation Analysis")
      ),
        card_body(
          layout_sidebar(
            sidebar = sidebar(
              title = "Analysis Settings",
              width = 280,

              # Info box matching other pN/pS tabs
              div(
                style = "background-color: #e7f3ff; padding: 8px; margin: 10px 0; border-radius: 4px; border-left: 3px solid #2196F3;",
                tags$small(
                  icon("info-circle"),
                  " Visualize correlation between Nucleotide Diversity and pN/pS by taxonomy level."
                )
              ),

              # Analysis message
              div(
                id = ns("analysis_message"),
                style = "color: red; font-weight: bold; margin-bottom: 10px;",
                "Please load integrated data first."
              ),

              accordion(
                accordion_panel(
                  title = "Analysis Parameters",
                  icon = bsicons::bs_icon("gear"),
                  open = TRUE,

                  selectInput(ns("tax_level"),
                             "Taxonomic Level:",
                             choices = c("Domain" = "domain",
                                        "Phylum" = "phylum",
                                        "Class" = "class",
                                        "Order" = "order",
                                        "Family" = "family",
                                        "Genus" = "genus",
                                        "Species" = "species"),
                             selected = "genus"),

                  selectInput(ns("group_var"),
                             "Group by Metadata:",
                             choices = NULL,
                             selected = NULL),

                  sliderInput(ns("top_n_taxa"),
                             "Number of Taxa to Display:",
                             min = 5, max = 50, value = 15, step = 1),

                  selectInput(ns("sort_by"),
                             "Sort Taxa by:",
                             choices = c("pN/pS p-value (significant first)" = "pnps_pvalue",
                                        "Sample Count" = "n_samples",
                                        "Mean Nucl. Diversity" = "mean_nucl_div",
                                        "Mean pN/pS" = "mean_pnps",
                                        "Alphabetical" = "alpha"),
                             selected = "pnps_pvalue")
                ),

                accordion_panel(
                  title = "Quality Filters",
                  icon = bsicons::bs_icon("funnel"),
                  open = FALSE,

                  sliderInput(ns("min_coverage"),
                             "Minimum Coverage:",
                             min = 1, max = 50, value = 5, step = 1),

                  sliderInput(ns("min_breadth"),
                             "Minimum Breadth:",
                             min = 0.1, max = 1.0, value = 0.5, step = 0.1),

                  sliderInput(ns("min_genomes"),
                             "Min Genomes per Taxa-Group:",
                             min = 1, max = 20, value = 1, step = 1),

                  helpText("Note: For Species level, use min_genomes = 1"),

                  # Filter status display (unified with Genome-wide pN/pS)
                  div(
                    style = "background-color: #f8f9fa; padding: 8px; margin: 10px 0; border-radius: 4px; border: 1px solid #dee2e6; font-size: 12px;",
                    tags$b(icon("filter"), " Applied Filters:"),
                    br(),
                    uiOutput(ns("filter_info"))
                  )
                ),

                accordion_panel(
                  title = "Plot Appearance",
                  icon = bsicons::bs_icon("palette"),
                  open = FALSE,

                  selectInput(ns("color_palette"),
                             "Color Palette:",
                             choices = c("Category20" = "category20",
                                        "Set1" = "Set1",
                                        "Set2" = "Set2",
                                        "Set3" = "Set3",
                                        "Paired" = "Paired",
                                        "Dark2" = "Dark2"),
                             selected = "category20"),

                  sliderInput(ns("point_size"),
                             "Point Size:",
                             min = 2, max = 15, value = 6, step = 1),

                  checkboxInput(ns("show_connecting_lines"),
                               "Connect Same Taxa (Dashed Lines)",
                               value = TRUE),

                  checkboxInput(ns("show_error_bars"),
                               "Show Error Bars (SE)",
                               value = FALSE)
                )
              ),

              br(),
              actionButton(ns("update_analysis"), "Update Analysis",
                          class = "btn-primary btn-block",
                          icon = icon("sync"))
            ),

            # Main content area - vertical layout
            div(
              style = "display: flex; flex-direction: column; gap: 15px; padding: 10px;",

              # Scatter Plot Card - larger height
              card(
                height = "550px",
                card_header(
                  div(
                    style = "display: flex; justify-content: space-between; align-items: center;",
                    span("Scatter Plot: Nucleotide Diversity vs pN/pS"),
                    create_plot_download_btn(ns, "scatter")
                  )
                ),
                card_body(
                  style = "padding: 5px;",
                  plotlyOutput(ns("scatter_plot"), height = "500px")
                )
              ),

              # Diversity Metrics Comparison Table
              card(
                height = "350px",
                card_header("Per-Taxon Wilcoxon Test (Groups Comparison)"),
                card_body(
                  style = "overflow-y: auto;",
                  DTOutput(ns("comparison_table"))
                )
              ),

              # Data Table (collapsible)
              accordion(
                accordion_panel(
                  title = "Aggregated Data Table",
                  icon = bsicons::bs_icon("table"),
                  div(
                    style = "max-height: 400px; overflow-y: auto;",
                    DTOutput(ns("data_table"))
                  )
                )
              )
            )
          )
        )
      )
    )
}

# ==============================================================================
# Server Function
# ==============================================================================

diversity_pnps_correlation_Server <- function(id, diversity_data, pnps_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    rv <- reactiveValues(
      merged_data = NULL,
      genome_means = NULL,
      aggregated_data = NULL,
      group_comparison = NULL,
      filter_stats = list(div_before = 0, div_after = 0, pnps_before = 0, merged = 0)
    )

    metadata_vars <- reactive({
      div_data <- diversity_data()
      if(is.null(div_data)) return(NULL)

      meta_cols <- c("disease_group", "host_sex", "country", "continent",
                    "age_group", "AJCC_stage")
      available <- intersect(meta_cols, names(div_data))

      categorical <- character(0)
      for(col in available) {
        if(col %in% names(div_data)) {
          vals <- unique(div_data[[col]])
          vals <- vals[!is.na(vals)]
          if(length(vals) >= 2 && length(vals) <= 20) {
            categorical <- c(categorical, col)
          }
        }
      }
      return(categorical)
    })

    observe({
      meta_vars <- metadata_vars()
      if(!is.null(meta_vars) && length(meta_vars) > 0) {
        updateSelectInput(session, "group_var",
                         choices = meta_vars,
                         selected = meta_vars[1])
        shinyjs::hide("analysis_message")
      } else {
        shinyjs::show("analysis_message")
      }
    })

    # Filter info display (unified with Genome-wide pN/pS)
    output$filter_info <- renderUI({
      div(
        tags$span(paste0("Coverage >= ", input$min_coverage)),
        br(),
        tags$span(paste0("Breadth >= ", input$min_breadth)),
        br(),
        tags$span(paste0("Min Genomes >= ", input$min_genomes)),
        br(),
        if(rv$filter_stats$merged > 0) {
          tags$span(
            style = "color: #6c757d;",
            paste0("Diversity: ", rv$filter_stats$div_before, " -> ", rv$filter_stats$div_after),
            br(),
            paste0("After inner_join: ", rv$filter_stats$merged)
          )
        } else {
          tags$span(style = "color: #6c757d;", "Run analysis to see filter stats")
        }
      )
    })

    observeEvent(input$update_analysis, {
      req(diversity_data(), pnps_data())

      div_data <- diversity_data()
      pnps_raw <- pnps_data()

      if(is.list(pnps_raw) && "genome_wide_pnps" %in% names(pnps_raw)) {
        pnps_genome <- pnps_raw$genome_wide_pnps
      } else {
        pnps_genome <- pnps_raw
      }

      if(is.null(div_data) || is.null(pnps_genome)) {
        showNotification("Data not available", type = "error")
        return()
      }

      tax_level <- input$tax_level
      group_var <- input$group_var

      cat("=== Diversity vs pN/pS Analysis ===\n")
      cat("Tax level:", tax_level, "\n")
      cat("Group var:", group_var, "\n")

      # Track filter statistics
      rv$filter_stats$div_before <- nrow(div_data)
      rv$filter_stats$pnps_before <- nrow(pnps_genome)

      # Quality filters
      div_filtered <- div_data %>%
        filter(coverage >= input$min_coverage,
               breadth >= input$min_breadth,
               !is.na(nucl_diversity),
               !is.na(.data[[tax_level]]),
               .data[[tax_level]] != "",
               !is.na(.data[[group_var]]))

      pnps_filtered <- pnps_genome %>%
        filter(!is.na(pNpS), pNpS > 0)

      rv$filter_stats$div_after <- nrow(div_filtered)

      cat("Filtered diversity rows:", nrow(div_filtered), "\n")
      cat("Filtered pnps rows:", nrow(pnps_filtered), "\n")

      # Merge
      div_select <- div_filtered %>%
        select(sample_id, clean_genome, nucl_diversity,
               all_of(c(tax_level, group_var)))

      pnps_select <- pnps_filtered %>%
        select(sample_id, clean_genome, pNpS)

      merged <- inner_join(div_select, pnps_select, by = c("sample_id", "clean_genome"))
      rv$filter_stats$merged <- nrow(merged)
      cat("Merged rows:", nrow(merged), "\n")

      if(nrow(merged) == 0) {
        showNotification("No matching data", type = "error")
        return()
      }

      rv$merged_data <- merged

      # Step 1: Mean per genome
      genome_means <- merged %>%
        group_by(clean_genome, .data[[tax_level]], .data[[group_var]]) %>%
        summarise(
          genome_mean_nucl_div = mean(nucl_diversity, na.rm = TRUE),
          genome_mean_pnps = mean(pNpS, na.rm = TRUE),
          n_samples_per_genome = n(),
          .groups = "drop"
        )

      rv$genome_means <- genome_means
      cat("Genome means:", nrow(genome_means), "\n")

      # Step 2: Mean per taxonomy + metadata group
      taxa_means <- genome_means %>%
        group_by(.data[[tax_level]], .data[[group_var]]) %>%
        summarise(
          mean_nucl_div = mean(genome_mean_nucl_div, na.rm = TRUE),
          mean_pnps = mean(genome_mean_pnps, na.rm = TRUE),
          se_nucl_div = if(n() > 1) sd(genome_mean_nucl_div, na.rm = TRUE) / sqrt(n()) else 0,
          se_pnps = if(n() > 1) sd(genome_mean_pnps, na.rm = TRUE) / sqrt(n()) else 0,
          n_genomes = n(),
          total_samples = sum(n_samples_per_genome),
          .groups = "drop"
        ) %>%
        filter(n_genomes >= input$min_genomes)

      cat("Taxa means (after min_genomes filter):", nrow(taxa_means), "\n")

      if(nrow(taxa_means) == 0) {
        showNotification("No taxa meet the minimum genome requirement. Try lowering 'Min Genomes per Taxa-Group'.", type = "warning")
        return()
      }

      # Step 3: Run Wilcoxon test on ALL taxa FIRST (for p-value based sorting)
      groups <- unique(merged[[group_var]])
      cat("Groups found:", paste(groups, collapse = ", "), "\n")

      all_wilcox_df <- NULL

      if(length(groups) == 2) {
        group1 <- groups[1]
        group2 <- groups[2]

        # Find ALL taxa present in both groups
        all_taxa_both <- merged %>%
          group_by(.data[[tax_level]]) %>%
          filter(length(unique(.data[[group_var]])) == 2) %>%
          pull(.data[[tax_level]]) %>%
          unique()

        cat("All taxa in both groups:", length(all_taxa_both), "\n")

        if(length(all_taxa_both) > 0) {
          wilcox_results <- lapply(all_taxa_both, function(taxon) {
            taxon_data <- merged %>% filter(.data[[tax_level]] == taxon)

            grp1_nucl <- taxon_data %>% filter(.data[[group_var]] == group1) %>% pull(nucl_diversity)
            grp2_nucl <- taxon_data %>% filter(.data[[group_var]] == group2) %>% pull(nucl_diversity)
            grp1_pnps <- taxon_data %>% filter(.data[[group_var]] == group1) %>% pull(pNpS)
            grp2_pnps <- taxon_data %>% filter(.data[[group_var]] == group2) %>% pull(pNpS)

            wilcox_nucl <- tryCatch({
              if(length(grp1_nucl) >= 2 && length(grp2_nucl) >= 2) {
                wilcox.test(grp1_nucl, grp2_nucl)
              } else {
                list(p.value = NA)
              }
            }, error = function(e) list(p.value = NA))

            wilcox_pnps <- tryCatch({
              if(length(grp1_pnps) >= 2 && length(grp2_pnps) >= 2) {
                wilcox.test(grp1_pnps, grp2_pnps)
              } else {
                list(p.value = NA)
              }
            }, error = function(e) list(p.value = NA))

            data.frame(
              Taxon = taxon,
              n_grp1 = length(grp1_nucl),
              n_grp2 = length(grp2_nucl),
              mean_nucl_div_grp1 = mean(grp1_nucl, na.rm = TRUE),
              mean_nucl_div_grp2 = mean(grp2_nucl, na.rm = TRUE),
              diff_nucl_div = mean(grp1_nucl, na.rm = TRUE) - mean(grp2_nucl, na.rm = TRUE),
              p_nucl_div = wilcox_nucl$p.value,
              mean_pnps_grp1 = mean(grp1_pnps, na.rm = TRUE),
              mean_pnps_grp2 = mean(grp2_pnps, na.rm = TRUE),
              diff_pnps = mean(grp1_pnps, na.rm = TRUE) - mean(grp2_pnps, na.rm = TRUE),
              p_pnps = wilcox_pnps$p.value,
              stringsAsFactors = FALSE
            )
          })

          all_wilcox_df <- bind_rows(wilcox_results)

          # Significance markers
          all_wilcox_df$sig_nucl_div <- ifelse(is.na(all_wilcox_df$p_nucl_div), "",
                                              ifelse(all_wilcox_df$p_nucl_div < 0.001, "***",
                                                    ifelse(all_wilcox_df$p_nucl_div < 0.01, "**",
                                                          ifelse(all_wilcox_df$p_nucl_div < 0.05, "*", "ns"))))
          all_wilcox_df$sig_pnps <- ifelse(is.na(all_wilcox_df$p_pnps), "",
                                          ifelse(all_wilcox_df$p_pnps < 0.001, "***",
                                                ifelse(all_wilcox_df$p_pnps < 0.01, "**",
                                                      ifelse(all_wilcox_df$p_pnps < 0.05, "*", "ns"))))

          # Rename columns with group names
          colnames(all_wilcox_df)[colnames(all_wilcox_df) == "n_grp1"] <- paste0("n_", group1)
          colnames(all_wilcox_df)[colnames(all_wilcox_df) == "n_grp2"] <- paste0("n_", group2)
          colnames(all_wilcox_df)[colnames(all_wilcox_df) == "mean_nucl_div_grp1"] <- paste0("NuclDiv_", group1)
          colnames(all_wilcox_df)[colnames(all_wilcox_df) == "mean_nucl_div_grp2"] <- paste0("NuclDiv_", group2)
          colnames(all_wilcox_df)[colnames(all_wilcox_df) == "mean_pnps_grp1"] <- paste0("pNpS_", group1)
          colnames(all_wilcox_df)[colnames(all_wilcox_df) == "mean_pnps_grp2"] <- paste0("pNpS_", group2)

          cat("Wilcoxon test completed for all taxa:", nrow(all_wilcox_df), "taxa\n")
        }
      }

      # Step 4: Select top N taxa (now can sort by p-value)
      taxa_overall <- taxa_means %>%
        group_by(.data[[tax_level]]) %>%
        summarise(
          total_samples_all = sum(total_samples),
          overall_mean_nucl_div = mean(mean_nucl_div, na.rm = TRUE),
          overall_mean_pnps = mean(mean_pnps, na.rm = TRUE),
          .groups = "drop"
        )

      # Join with p-values if available
      if(!is.null(all_wilcox_df)) {
        taxa_overall <- taxa_overall %>%
          left_join(all_wilcox_df %>% select(Taxon, p_pnps, p_nucl_div),
                   by = setNames("Taxon", tax_level))
      }

      # Sort by selected criteria
      if(input$sort_by == "pnps_pvalue" && !is.null(all_wilcox_df)) {
        # Sort by p-value (smallest first = most significant)
        taxa_overall <- taxa_overall %>%
          arrange(p_pnps, desc(total_samples_all))
        cat("Sorting by pN/pS p-value (significant first)\n")
      } else if(input$sort_by == "n_samples") {
        taxa_overall <- taxa_overall %>% arrange(desc(total_samples_all))
      } else if(input$sort_by == "mean_nucl_div") {
        taxa_overall <- taxa_overall %>% arrange(desc(overall_mean_nucl_div))
      } else if(input$sort_by == "mean_pnps") {
        taxa_overall <- taxa_overall %>% arrange(desc(overall_mean_pnps))
      } else {
        taxa_overall <- taxa_overall %>% arrange(.data[[tax_level]])
      }

      top_taxa <- head(taxa_overall[[tax_level]], input$top_n_taxa)

      aggregated <- taxa_means %>%
        filter(.data[[tax_level]] %in% top_taxa)

      aggregated[[tax_level]] <- factor(aggregated[[tax_level]], levels = top_taxa)

      rv$aggregated_data <- aggregated
      cat("Final aggregated:", nrow(aggregated), "points,", length(top_taxa), "taxa\n")

      # Step 5: Filter Wilcoxon results to top_taxa only
      if(!is.null(all_wilcox_df)) {
        rv$group_comparison <- all_wilcox_df %>%
          filter(Taxon %in% top_taxa)
        cat("Wilcoxon comparison table (top taxa):", nrow(rv$group_comparison), "\n")
      } else {
        rv$group_comparison <- NULL
      }

      showNotification(
        paste("Analysis complete:", nrow(aggregated), "points from", length(top_taxa), "taxa"),
        type = "message"
      )
    })

    # Scatter plot - using plotly directly for separate legends
    output$scatter_plot <- renderPlotly({
      req(rv$aggregated_data)

      data <- rv$aggregated_data
      tax_level <- input$tax_level
      group_var <- input$group_var

      if(nrow(data) == 0) {
        return(plotly_empty() %>% layout(title = "No data available"))
      }

      taxa_levels <- levels(data[[tax_level]])
      n_taxa <- length(taxa_levels)
      group_levels <- unique(as.character(data[[group_var]]))
      n_groups <- length(group_levels)

      # Color palette for taxa
      if(input$color_palette == "category20") {
        colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                   "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                   "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
        color_map <- setNames(rep(colors, length.out = n_taxa), taxa_levels)
      } else {
        pal_colors <- brewer.pal(max(3, min(n_taxa, 8)), input$color_palette)
        color_map <- setNames(colorRampPalette(pal_colors)(n_taxa), taxa_levels)
      }

      # Shape mapping for metadata groups (plotly symbols)
      # circle, triangle-up, square, diamond, x, cross
      shape_map <- setNames(c("circle", "triangle-up", "square", "diamond", "x", "cross")[1:n_groups], group_levels)

      # Axis limits
      x_range <- range(data$mean_nucl_div, na.rm = TRUE)
      y_range <- range(data$mean_pnps, na.rm = TRUE)
      x_margin <- diff(x_range) * 0.1
      y_margin <- diff(y_range) * 0.1

      # Initialize plot
      p <- plot_ly()

      # Get significant taxa from Wilcoxon results (using raw p-values)
      sig_taxa <- c()
      if(!is.null(rv$group_comparison)) {
        sig_taxa <- rv$group_comparison %>%
          filter(p_nucl_div < 0.05 | p_pnps < 0.05) %>%
          pull(Taxon)
      }

      # Add connecting lines first (if enabled)
      # Solid line = p < 0.05, Dashed line = ns
      if(input$show_connecting_lines && n_groups > 1) {
        taxa_in_multiple <- data %>%
          group_by(.data[[tax_level]]) %>%
          filter(n() > 1) %>%
          ungroup()

        if(nrow(taxa_in_multiple) > 0) {
          for(taxon in unique(taxa_in_multiple[[tax_level]])) {
            taxon_data <- taxa_in_multiple %>%
              filter(.data[[tax_level]] == taxon) %>%
              arrange(.data[[group_var]])

            # Determine line style based on significance
            is_sig <- as.character(taxon) %in% sig_taxa
            line_dash <- ifelse(is_sig, "solid", "dash")
            line_width <- ifelse(is_sig, 2, 1)
            line_opacity <- ifelse(is_sig, 0.8, 0.4)

            p <- p %>% add_trace(
              data = taxon_data,
              x = ~mean_nucl_div,
              y = ~mean_pnps,
              type = "scatter",
              mode = "lines",
              line = list(color = color_map[as.character(taxon)], dash = line_dash, width = line_width),
              opacity = line_opacity,
              showlegend = FALSE,
              hoverinfo = "skip"
            )
          }
        }
      }

      # Add points - one trace per metadata group (for shape legend)
      for(grp in group_levels) {
        grp_data <- data %>% filter(.data[[group_var]] == grp)

        p <- p %>% add_trace(
          data = grp_data,
          x = ~mean_nucl_div,
          y = ~mean_pnps,
          type = "scatter",
          mode = "markers",
          marker = list(
            symbol = shape_map[grp],
            size = input$point_size * 2,
            color = color_map[as.character(grp_data[[tax_level]])],
            line = list(color = "white", width = 0.5)
          ),
          text = ~paste0(
            "Taxon: ", .data[[tax_level]], "<br>",
            "Group: ", .data[[group_var]], "<br>",
            "Nucl.Div: ", signif(mean_nucl_div, 4), "<br>",
            "pN/pS: ", signif(mean_pnps, 4), "<br>",
            "N genomes: ", n_genomes
          ),
          hoverinfo = "text",
          name = grp,
          legendgroup = "shape",
          showlegend = TRUE
        )
      }

      # Layout
      p <- p %>% layout(
        title = list(
          text = paste("Nucleotide Diversity vs pN/pS by", tax_level),
          font = list(size = 16)
        ),
        xaxis = list(
          title = "Mean Nucleotide Diversity",
          range = c(x_range[1] - x_margin, x_range[2] + x_margin)
        ),
        yaxis = list(
          title = "Mean pN/pS",
          range = c(y_range[1] - y_margin, y_range[2] + y_margin)
        ),
        legend = list(
          title = list(text = paste0("<b>", group_var, "</b>")),
          orientation = "v",
          x = 1.05,
          y = 1,
          font = list(size = 11)
        ),
        margin = list(r = 220, t = 50, b = 50, l = 60)
      )

      # Add annotations for additional legends
      legend_y_start <- 0.92
      legend_y_step <- 0.05

      # Add "Line" legend section
      line_legend_y <- legend_y_start - (n_groups + 1) * legend_y_step
      p <- p %>% add_annotations(
        x = 1.05, y = line_legend_y,
        xref = "paper", yref = "paper",
        text = "<b>Connection</b>",
        showarrow = FALSE,
        font = list(size = 11),
        xanchor = "left"
      )
      p <- p %>% add_annotations(
        x = 1.05, y = line_legend_y - legend_y_step,
        xref = "paper", yref = "paper",
        text = "── p < 0.05 (sig)",
        showarrow = FALSE,
        font = list(size = 9),
        xanchor = "left"
      )
      p <- p %>% add_annotations(
        x = 1.05, y = line_legend_y - 2 * legend_y_step,
        xref = "paper", yref = "paper",
        text = "- - - p ≥ 0.05 (ns)",
        showarrow = FALSE,
        font = list(size = 9),
        xanchor = "left"
      )

      # Add color legend title
      color_legend_y <- line_legend_y - 3.5 * legend_y_step
      p <- p %>% add_annotations(
        x = 1.05, y = color_legend_y,
        xref = "paper", yref = "paper",
        text = paste0("<b>", tax_level, "</b>"),
        showarrow = FALSE,
        font = list(size = 11),
        xanchor = "left"
      )

      # Add color swatches for each taxon
      for(i in seq_along(taxa_levels)) {
        taxon <- taxa_levels[i]
        y_pos <- color_legend_y - i * legend_y_step

        if(y_pos > -0.3) {
          p <- p %>% add_annotations(
            x = 1.05, y = y_pos,
            xref = "paper", yref = "paper",
            text = paste0("<span style='color:", color_map[taxon], ";'>●</span> ", taxon),
            showarrow = FALSE,
            font = list(size = 9),
            xanchor = "left"
          )
        }
      }

      p
    })

    # Per-taxon Wilcoxon comparison table
    output$comparison_table <- renderDT({
      req(rv$group_comparison)

      data <- rv$group_comparison

      # Format numeric columns
      numeric_cols <- c("diff_nucl_div", "diff_pnps", "p_nucl_div", "p_pnps")
      for(col in numeric_cols) {
        if(col %in% names(data)) data[[col]] <- signif(data[[col]], 3)
      }

      mean_cols <- grep("NuclDiv_|pNpS_", names(data), value = TRUE)
      for(col in mean_cols) {
        if(col %in% names(data)) data[[col]] <- signif(data[[col]], 4)
      }

      datatable(data,
                options = list(
                  pageLength = 15,
                  scrollX = TRUE,
                  scrollY = "250px",
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel'),
                  columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ),
                rownames = FALSE,
                extensions = 'Buttons',
                caption = "Per-taxon Wilcoxon test (raw p-values). Solid line = p < 0.05, dashed = ns.") %>%
        formatStyle('sig_nucl_div',
                   backgroundColor = styleEqual(c("***", "**", "*", "ns"),
                                               c('#90EE90', '#87CEEB', '#FFD700', 'white'))) %>%
        formatStyle('sig_pnps',
                   backgroundColor = styleEqual(c("***", "**", "*", "ns"),
                                               c('#90EE90', '#87CEEB', '#FFD700', 'white')))
    })

    # Data table
    output$data_table <- renderDT({
      req(rv$aggregated_data)

      tax_level <- input$tax_level
      group_var <- input$group_var

      data <- rv$aggregated_data %>%
        arrange(desc(n_genomes)) %>%
        mutate(
          mean_nucl_div = signif(mean_nucl_div, 4),
          mean_pnps = signif(mean_pnps, 4),
          se_nucl_div = signif(se_nucl_div, 3),
          se_pnps = signif(se_pnps, 3)
        )

      colnames(data) <- c(tax_level, group_var, "Mean Nucl.Div", "Mean pN/pS",
                         "SE Nucl.Div", "SE pN/pS", "N Genomes", "Total Samples")

      datatable(data,
                options = list(
                  pageLength = 15,
                  scrollX = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel')
                ),
                rownames = FALSE,
                extensions = 'Buttons',
                filter = 'top')
    })

    # ========== Download Modal Handler ==========
    scatter_ggplot <- reactive({
      req(rv$aggregated_data)

      data <- rv$aggregated_data
      tax_level <- input$tax_level
      group_var <- input$group_var

      req(nrow(data) > 0)

      # Get taxa levels and groups
      taxa_levels <- levels(data[[tax_level]])
      if(is.null(taxa_levels)) taxa_levels <- unique(as.character(data[[tax_level]]))
      n_taxa <- length(taxa_levels)
      group_levels <- unique(as.character(data[[group_var]]))

      # Color palette matching interactive plot
      if(input$color_palette == "category20") {
        colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                   "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                   "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
        color_values <- setNames(rep(colors, length.out = n_taxa), taxa_levels)
      } else {
        pal_colors <- brewer.pal(max(3, min(n_taxa, 8)), input$color_palette)
        color_values <- setNames(colorRampPalette(pal_colors)(n_taxa), taxa_levels)
      }

      # Get significant taxa from Wilcoxon results
      sig_taxa <- c()
      if(!is.null(rv$group_comparison)) {
        sig_taxa <- rv$group_comparison %>%
          filter(p_nucl_div < 0.05 | p_pnps < 0.05) %>%
          pull(Taxon)
      }

      # Base plot
      p <- ggplot(data, aes(x = mean_nucl_div, y = mean_pnps))

      # Add connecting lines first (if enabled) with solid/dashed based on significance
      if (input$show_connecting_lines && length(group_levels) > 1) {
        taxa_in_multiple <- data %>%
          group_by(.data[[tax_level]]) %>%
          filter(n() > 1) %>%
          ungroup()

        if (nrow(taxa_in_multiple) > 0) {
          # Separate significant and non-significant taxa
          sig_taxa_data <- taxa_in_multiple %>%
            filter(as.character(.data[[tax_level]]) %in% sig_taxa)
          nonsig_taxa_data <- taxa_in_multiple %>%
            filter(!as.character(.data[[tax_level]]) %in% sig_taxa)

          # Add solid lines for significant taxa
          if(nrow(sig_taxa_data) > 0) {
            p <- p + geom_line(data = sig_taxa_data,
                              aes(group = .data[[tax_level]], color = .data[[tax_level]]),
                              linetype = "solid", linewidth = 1, alpha = 0.7)
          }
          # Add dashed lines for non-significant taxa
          if(nrow(nonsig_taxa_data) > 0) {
            p <- p + geom_line(data = nonsig_taxa_data,
                              aes(group = .data[[tax_level]], color = .data[[tax_level]]),
                              linetype = "dashed", linewidth = 0.8, alpha = 0.4)
          }
        }
      }

      # Add error bars if enabled
      if (input$show_error_bars) {
        p <- p +
          geom_errorbar(aes(ymin = mean_pnps - se_pnps, ymax = mean_pnps + se_pnps,
                           color = .data[[tax_level]]),
                       width = 0, alpha = 0.5) +
          geom_errorbarh(aes(xmin = mean_nucl_div - se_nucl_div, xmax = mean_nucl_div + se_nucl_div,
                            color = .data[[tax_level]]),
                        height = 0, alpha = 0.5)
      }

      # Add points
      p <- p +
        geom_point(aes(color = .data[[tax_level]], shape = .data[[group_var]]),
                  size = input$point_size, alpha = 0.85) +
        scale_color_manual(values = color_values, name = tax_level) +
        scale_shape_manual(values = c(16, 17, 15, 18, 8, 3)[1:length(group_levels)],
                          name = group_var) +
        labs(title = paste("Nucleotide Diversity vs pN/pS by", tax_level),
             subtitle = ifelse(input$show_connecting_lines,
                              "Solid line = p < 0.05, Dashed line = ns",
                              NULL),
             x = "Mean Nucleotide Diversity",
             y = "Mean pN/pS") +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray40"),
          legend.position = "right",
          legend.box = "vertical",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid.minor = element_blank()
        ) +
        guides(
          color = guide_legend(order = 1, ncol = 1),
          shape = guide_legend(order = 2)
        )

      p
    })

    setup_download_modal_handler(input, output, session, "scatter", scatter_ggplot, "diversity_pnps_correlation")

  })
}
