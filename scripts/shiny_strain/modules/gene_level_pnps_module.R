# gene_level_pnps_module.R
# Shiny module for gene-level pN/pS analysis with COG category integration

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(DT)
library(viridis)
library(parallel)  # For parallel processing

# ==============================================================================
# UI Function
# ==============================================================================

gene_level_pnps_UI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),
    # Custom header tabs (clickable) and hide the default tab bar
    tags$head(tags$style(HTML(paste0(
      "#", ns("main_tabs"), " .nav,",
      "#", ns("main_tabs"), " .nav.nav-tabs,",
      "#", ns("main_tabs"), " .nav.nav-pills,",
      "#", ns("main_tabs"), " .tabbable > .nav { display: none !important; }",
      "#", ns("header_tabs"), " .nav-link{ cursor:pointer; color:#ffffff; font-weight:600; }",
      "#", ns("header_tabs"), " .nav-link:hover{ text-decoration: underline; }",
      "#", ns("header_tabs"), " .nav-link.active{ text-decoration: underline; }",
      "#", ns("header_tabs"), " { margin-top: 6px; }"
    )))),

    card(
      full_screen = TRUE,
      class = "main-card",
      fill = TRUE,
      card_header(
        style = "background-color: #97C6F1; color: white;",
        # Tabs placed in the blue header
        tags$ul(id = ns("header_tabs"), style = "list-style:none; padding-left:0; margin:6px 0 0 0; display:flex; gap:16px;",
          tags$li(class = "nav-item", actionLink(ns("tab_visualization"), label = "Visualization", class = "nav-link active")),
          tags$li(class = "nav-item", actionLink(ns("tab_significance"), label = "Statistical Significance", class = "nav-link")),
          tags$li(class = "nav-item", actionLink(ns("tab_statistics"), label = "Summary Statistics", class = "nav-link")),
          tags$li(class = "nav-item", actionLink(ns("tab_data"), label = "Gene Details", class = "nav-link"))
        )
      ),
      card_body(
        layout_sidebar(
          sidebar = sidebar(
            title = "Gene-Level Analysis Settings",

            # Instructional text
            div(
              style = "background-color: #e7f3ff; padding: 10px; margin: 10px 0; border-radius: 4px; border-left: 3px solid #2196F3;",
              tags$small(
                icon("info-circle"),
                " Analyze gene-level pN/pS by functional categories (COG) and taxonomy."
              )
            ),

            # Red instruction text
            div(
              style = "background-color: #ffebee; padding: 10px; margin: 10px 0; border-radius: 4px; border-left: 3px solid #f44336; color: #c62828;",
              tags$small(
                icon("exclamation-triangle"),
                tags$b(" To change visualization settings, adjust parameters below and click 'Update Analysis'.")
              )
            ),

            # Data loading button
            actionButton(ns("load_data"),
                        "Load Gene-Level Data",
                        class = "btn-success btn-block",
                        icon = icon("database")),

            br(),

            # Taxonomic level selection - moved here for easier access
            selectInput(ns("tax_level"),
                       "Taxonomic Level:",
                       choices = c("Phylum" = "phylum",
                                  "Class" = "class",
                                  "Order" = "order",
                                  "Family" = "family",
                                  "Genus" = "genus",
                                  "Species" = "species"),
                       selected = "phylum"),

            # Metadata variable selection (dynamic) - moved here for easier access
            uiOutput(ns("metadata_var_ui")),

            # Update Analysis button moved up
            actionButton(ns("run_analysis"),
                        "Update Analysis",
                        class = "btn-primary btn-block",
                        icon = icon("play")),

            br(),

            accordion(
              accordion_panel(
                title = "Analysis Parameters",
                icon = bsicons::bs_icon("gear"),

                # Analysis type - fixed to COG (KEGG disabled)
                selectInput(ns("analysis_type"),
                           "Analysis Type:",
                           choices = c("COG Categories" = "cog"),
                           selected = "cog"),

                # Visualization type selection
                radioButtons(ns("viz_type"),
                            "Visualization Type:",
                            choices = c("Heatmap" = "heatmap",
                                       "Bar Chart" = "bar"),
                            selected = "heatmap"),

                # COG category selection (conditional)
                conditionalPanel(
                  condition = "input.analysis_type == 'cog'",
                  ns = ns,
                  selectInput(ns("cog_display"),
                             "COG Display:",
                             choices = c("Individual Categories (A, B, C, ...)" = "individual",
                                        "Functional Groups" = "grouped"),
                             selected = "individual"),

                  selectInput(ns("selected_cog"),
                             "Select COG Category:",
                             choices = c("All" = "all"),
                             selected = "all",
                             multiple = TRUE)
                ),

                # KEGG pathway selection (disabled - hidden)
                # conditionalPanel(
                #   condition = "input.analysis_type == 'kegg'",
                #   ns = ns,
                #   selectInput(ns("selected_kegg"),
                #              "Select KEGG Pathway:",
                #              choices = c("All" = "all"),
                #              selected = "all",
                #              multiple = TRUE)
                # ),

                # Number of top taxa to show (dynamic max via uiOutput)
                uiOutput(ns("n_top_taxa_ui")),

                # Number of top functional categories to show
                sliderInput(ns("n_top_categories"),
                           "Max functional categories:",
                           min = 10, max = 200, value = 100, step = 5),

                # Transformation
                radioButtons(ns("transformation"),
                            "pN/pS Transformation:",
                            choices = c("Raw" = "raw", "Log10" = "log"),
                            selected = "raw"),

                # Filters
                sliderInput(ns("min_coverage"),
                           "Min Coverage:",
                           min = 1, max = 100, value = 5, step = 1),

                sliderInput(ns("min_breadth"),
                           "Min Breadth:",
                           min = 0, max = 1, value = 0.5, step = 0.05),

                # CPU cores for parallel processing
                sliderInput(ns("n_cores"),
                           "CPU cores for statistics:",
                           min = 1, max = parallel::detectCores(),
                           value = min(16, parallel::detectCores()), step = 1),

                helpText("More cores = faster significance matrix calculation"),

                # Taxa sorting option
                selectInput(ns("taxa_sort_by"),
                           "Sort Taxa by:",
                           choices = c("Significance Count (most * first)" = "sig_count",
                                      "Gene Count (abundance)" = "gene_count",
                                      "Alphabetical" = "alpha"),
                           selected = "sig_count")
              )
            )
          ),

          # Main content with hidden tabs (controlled by header)
          tabsetPanel(
            id = ns("main_tabs"),
            type = "hidden",

            # Tab 1: Visualization
            tabPanelBody("visualization",
              card(
                card_header(
                  div(
                    style = "display: flex; justify-content: space-between; align-items: center;",
                    span("Gene-Level pN/pS by Functional Category and Taxonomy"),
                    create_plot_download_btn(ns, "gene_plot")
                  )
                ),
                card_body(
                  style = "overflow-y: auto; max-height: 2500px;",
                  plotlyOutput(ns("gene_plot"), height = "auto")
                )
              )
            ),

            # Tab 2: Summary Statistics
            tabPanelBody("statistics",
              card(
                card_header("Summary by Functional Category and Taxonomy"),
                card_body(
                  DTOutput(ns("cog_summary_table"))
                )
              )
            ),

            # Tab 3: Statistical Significance
            tabPanelBody("significance",
              card(
                card_header(
                  div(
                    style = "display: flex; justify-content: space-between; align-items: center;",
                    span("Statistical Significance Matrix (p-values by COG Category × Taxonomy)"),
                    create_plot_download_btn(ns, "significance_heatmap")
                  )
                ),
                card_body(
                  div(
                    style = "background-color: #e7f3ff; padding: 10px; margin-bottom: 15px; border-radius: 4px; border-left: 3px solid #2196F3;",
                    tags$small(
                      icon("info-circle"),
                      " This matrix shows statistical significance of pN/pS differences between metadata groups.",
                      br(),
                      "• Binary metadata (2 groups): Wilcoxon test",
                      br(),
                      "• Categorical metadata (>2 groups): Kruskal-Wallis test",
                      br(),
                      "• Numeric metadata: Spearman correlation",
                      br(),
                      "Significance levels: *** p<0.001, ** p<0.01, * p<0.05",
                      br(), br(),
                      tags$b("Tip: "), "Drag to select a region to zoom in. Double-click to reset view."
                    )
                  ),
                  plotlyOutput(ns("significance_heatmap"), height = "auto")
                )
              )
            ),

            # Tab 4: Gene Details
            tabPanelBody("data",
              card(
                card_header("Individual Gene Data"),
                card_body(
                  DTOutput(ns("gene_details_table"))
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

gene_level_pnps_Server <- function(id, gtdb_data = reactive(NULL), sample_metadata = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ===========================================================================
    # Tab Switching Logic (Header Tabs)
    # ===========================================================================

    # Switch to Visualization tab
    observeEvent(input$tab_visualization, {
      updateTabsetPanel(session, "main_tabs", selected = "visualization")
      shinyjs::runjs(sprintf("$('#%s .nav-link').removeClass('active'); $('#%s').addClass('active');",
                            ns("header_tabs"), ns("tab_visualization")))
    })

    # Switch to Summary Statistics tab
    observeEvent(input$tab_statistics, {
      updateTabsetPanel(session, "main_tabs", selected = "statistics")
      shinyjs::runjs(sprintf("$('#%s .nav-link').removeClass('active'); $('#%s').addClass('active');",
                            ns("header_tabs"), ns("tab_statistics")))
    })

    # Switch to Statistical Significance tab
    observeEvent(input$tab_significance, {
      updateTabsetPanel(session, "main_tabs", selected = "significance")
      shinyjs::runjs(sprintf("$('#%s .nav-link').removeClass('active'); $('#%s').addClass('active');",
                            ns("header_tabs"), ns("tab_significance")))
    })

    # Switch to Gene Details tab
    observeEvent(input$tab_data, {
      updateTabsetPanel(session, "main_tabs", selected = "data")
      shinyjs::runjs(sprintf("$('#%s .nav-link').removeClass('active'); $('#%s').addClass('active');",
                            ns("header_tabs"), ns("tab_data")))
    })

    # ===========================================================================
    # Reactive Values
    # ===========================================================================

    rv <- reactiveValues(
      gene_data = NULL,
      filtered_data = NULL,  # Aggregated summary data
      raw_data = NULL,        # Filtered raw gene-level data
      data_loaded = FALSE,
      metadata_available = NULL  # Store available metadata columns
    )

    # Render metadata variable dropdown dynamically (most robust approach)
    output$metadata_var_ui <- renderUI({
      cat("\n=== Rendering Metadata Dropdown ===\n")

      # Default dropdown before data is loaded
      if (is.null(rv$gene_data)) {
        cat("Gene data not loaded yet - showing default dropdown\n")
        return(selectInput(ns("metadata_var"),
                          "Metadata Variable:",
                          choices = c("None (load data first)" = "none"),
                          selected = "none"))
      }

      # Detect metadata variables
      exclude_cols <- c("gene", "sample_id", "clean_genome", "pNpS", "pN", "pS",
                       "coverage", "breadth", "gene_length", "breadth_minCov",
                       "COG_category", "COG_category_name", "COG_simplified",
                       "KEGG_Pathway", "KEGG_ko",
                       "domain", "phylum", "class", "order", "family", "genus", "species",
                       "gtdb_taxonomy", "accession_used_in_analysis",
                       "genome", "scaffold", "SNV_N_count", "SNV_S_count",
                       "N_sites", "S_sites", "Description", "Preferred_name", "bioproject_accession")

      all_cols <- names(rv$gene_data)
      meta_cols <- setdiff(all_cols, exclude_cols)

      cat("Detected", length(meta_cols), "metadata columns\n")

      if (length(meta_cols) > 0) {
        cat("Metadata variables:", paste(meta_cols, collapse = ", "), "\n")
        rv$metadata_available <- meta_cols

        # Create dropdown with metadata options
        selectInput(ns("metadata_var"),
                   "Metadata Variable:",
                   choices = c("None" = "none", setNames(meta_cols, meta_cols)),
                   selected = "none")
      } else {
        cat("No metadata found - showing default dropdown\n")
        rv$metadata_available <- NULL

        selectInput(ns("metadata_var"),
                   "Metadata Variable:",
                   choices = c("None (no metadata available)" = "none"),
                   selected = "none")
      }
    })

    # Dynamic top taxa slider based on actual taxa count
    output$n_top_taxa_ui <- renderUI({
      if (is.null(rv$gene_data)) {
        return(sliderInput(ns("n_top_taxa"),
                          "Top taxa to display:",
                          min = 5, max = 100, value = 15, step = 5))
      }

      tax_col <- input$tax_level
      if (is.null(tax_col) || !tax_col %in% names(rv$gene_data)) {
        return(sliderInput(ns("n_top_taxa"),
                          "Top taxa to display:",
                          min = 5, max = 100, value = 15, step = 5))
      }

      # Count distinct taxa
      n_distinct_taxa <- length(unique(rv$gene_data[[tax_col]]))
      max_taxa <- max(5, n_distinct_taxa)
      current_val <- min(15, max_taxa)

      sliderInput(ns("n_top_taxa"),
                 paste0("Top taxa to display (max: ", n_distinct_taxa, "):"),
                 min = 5, max = max_taxa, value = current_val, step = 1)
    })

    # Load gene-level data with COG annotations (OPTIMIZED)
    observeEvent(input$load_data, {
      # Check if GTDB data and metadata are available
      gtdb <- gtdb_data()
      metadata <- sample_metadata()

      if (is.null(gtdb) || nrow(gtdb) == 0) {
        showNotification(
          "⚠️ GTDB data not loaded! Please click 'Load All Data' in Master Controller first.",
          type = "warning",
          duration = 10
        )
        cat("\n❌ ERROR: Attempted to load gene data without GTDB!\n")
        cat("User must click 'Load All Data' in Master Controller first.\n\n")
        return()
      }

      if (is.null(metadata) || nrow(metadata) == 0) {
        showNotification(
          "⚠️ Sample metadata not loaded! Please click 'Load All Data' in Master Controller first.",
          type = "warning",
          duration = 10
        )
        cat("\n❌ ERROR: Attempted to load gene data without sample metadata!\n")
        cat("User must click 'Load All Data' in Master Controller first.\n\n")
        return()
      }

      withProgress(message = "Loading gene-level data...", {
        tryCatch({
          # Source optimized preprocessing
          source("modules/gene_level_pnps_preprocessing_optimized.R")

          # Load data using optimized path
          setProgress(0.3, detail = "Loading filtered gene data...")
          rv$gene_data <- load_gene_pnps_optimized(use_filtered = TRUE)

          # Fallback message if no optimized files found
          if (nrow(rv$gene_data) > 5000000) {
            showNotification(
              "TIP: Run scripts/create_filtered_gene_pnps.R for faster loading!",
              type = "message",
              duration = 8
            )
          }

          # Add GTDB taxonomy data
          setProgress(0.5, detail = "Adding taxonomic information...")
          gtdb <- gtdb_data()

          cat("\n=== GTDB Join Diagnostics ===\n")
          cat("GTDB data is NULL?", is.null(gtdb), "\n")
          if (!is.null(gtdb)) {
            cat("GTDB rows:", nrow(gtdb), "\n")
            cat("GTDB columns:", paste(names(gtdb), collapse = ", "), "\n")
            cat("Has gtdb_taxonomy column?", "gtdb_taxonomy" %in% names(gtdb), "\n")
            cat("Has accession column?", "accession" %in% names(gtdb), "\n")
            if ("accession" %in% names(gtdb)) {
              cat("Sample GTDB accessions:", paste(head(gtdb$accession, 3), collapse = ", "), "\n")
            }
          }

          # Show sample gene data
          cat("\nGene data sample:\n")
          cat("Gene data rows:", nrow(rv$gene_data), "\n")
          cat("Gene data columns:", paste(names(rv$gene_data), collapse = ", "), "\n")
          cat("Sample clean_genome values:", paste(head(unique(rv$gene_data$clean_genome), 3), collapse = ", "), "\n")
          cat("=============================\n\n")

          # Check if taxonomy columns already exist in the data
          has_taxonomy <- all(c("phylum", "class", "order", "family", "genus", "species") %in% names(rv$gene_data)) &&
                          sum(!is.na(rv$gene_data$phylum)) > 0

          if (has_taxonomy) {
            cat("Taxonomy columns already present in data - skipping GTDB join\n")
            n_with_tax <- sum(!is.na(rv$gene_data$phylum))
            cat("Existing taxonomy coverage:", n_with_tax, "genes with phylum\n")
          } else if (!is.null(gtdb) && nrow(gtdb) > 0 && "gtdb_taxonomy" %in% names(gtdb)) {
            tryCatch({
              cat("Processing GTDB taxonomy...\n")

              # Get unique genomes to minimize processing
              unique_genomes <- unique(rv$gene_data$clean_genome)
              cat("Found", length(unique_genomes), "unique genomes in gene data\n")

              # Try direct match first
              gtdb_filtered <- gtdb %>%
                filter(accession %in% unique_genomes) %>%
                select(accession, gtdb_taxonomy)

              cat("Direct match:", nrow(gtdb_filtered), "GTDB entries\n")

              # If no matches, try with prefix removal (RS_/GB_)
              if (nrow(gtdb_filtered) == 0) {
                cat("Trying with RS_/GB_ prefix removal...\n")
                gtdb_clean <- gtdb %>%
                  mutate(accession_clean = sub("^(RS_|GB_)", "", accession))

                gtdb_filtered <- gtdb_clean %>%
                  filter(accession_clean %in% unique_genomes) %>%
                  select(accession = accession_clean, gtdb_taxonomy)

                cat("After prefix removal:", nrow(gtdb_filtered), "GTDB entries\n")
              }

              # If still no matches, try matching gene data WITH prefixes
              if (nrow(gtdb_filtered) == 0) {
                cat("Trying to add RS_/GB_ prefixes to gene data...\n")
                gene_with_prefix <- c(
                  paste0("RS_", unique_genomes),
                  paste0("GB_", unique_genomes)
                )

                gtdb_filtered <- gtdb %>%
                  filter(accession %in% gene_with_prefix) %>%
                  select(accession, gtdb_taxonomy) %>%
                  mutate(accession = sub("^(RS_|GB_)", "", accession))

                cat("After adding prefixes:", nrow(gtdb_filtered), "GTDB entries\n")
              }

              cat("Final matched entries:", nrow(gtdb_filtered), "\n")
              if (nrow(gtdb_filtered) > 0) {
                cat("Sample matched accessions:", paste(head(gtdb_filtered$accession, 3), collapse = ", "), "\n")
              }

              if (nrow(gtdb_filtered) > 0) {
                # Parse taxonomy string properly by splitting on semicolon
                tax_split <- strsplit(gtdb_filtered$gtdb_taxonomy, ";", fixed = TRUE)

                gtdb_filtered <- gtdb_filtered %>%
                  mutate(
                    domain = sapply(tax_split, function(x) if (length(x) >= 1) gsub("^d__", "", x[1]) else NA_character_),
                    phylum = sapply(tax_split, function(x) if (length(x) >= 2) gsub("^p__", "", x[2]) else NA_character_),
                    class = sapply(tax_split, function(x) if (length(x) >= 3) gsub("^c__", "", x[3]) else NA_character_),
                    order = sapply(tax_split, function(x) if (length(x) >= 4) gsub("^o__", "", x[4]) else NA_character_),
                    family = sapply(tax_split, function(x) if (length(x) >= 5) gsub("^f__", "", x[5]) else NA_character_),
                    genus = sapply(tax_split, function(x) if (length(x) >= 6) gsub("^g__", "", x[6]) else NA_character_),
                    species = sapply(tax_split, function(x) if (length(x) >= 7) gsub("^s__", "", x[7]) else NA_character_)
                  )

                # Join with gene data
                rv$gene_data <- rv$gene_data %>%
                  left_join(gtdb_filtered, by = c("clean_genome" = "accession"))

                n_with_tax <- sum(!is.na(rv$gene_data$phylum))
                cat("Added taxonomy to", n_with_tax, "genes\n")
              } else {
                cat("WARNING: No matching GTDB entries found!\n")
                rv$gene_data <- rv$gene_data %>%
                  mutate(domain = NA, phylum = NA, class = NA, order = NA,
                         family = NA, genus = NA, species = NA, gtdb_taxonomy = NA)
              }
            }, error = function(e) {
              cat("ERROR in GTDB join:", e$message, "\n")
              showNotification(
                paste("GTDB join error:", e$message),
                type = "error",
                duration = 10
              )
              rv$gene_data <- rv$gene_data %>%
                mutate(domain = NA, phylum = NA, class = NA, order = NA,
                       family = NA, genus = NA, species = NA, gtdb_taxonomy = NA)
            })
          } else {
            cat("GTDB data not available or empty\n")
            showNotification(
              "GTDB data not available. Load data in Master Controller first.",
              type = "warning",
              duration = 5
            )
            # Add empty taxonomy columns
            rv$gene_data <- rv$gene_data %>%
              mutate(domain = NA, phylum = NA, class = NA, order = NA,
                     family = NA, genus = NA, species = NA, gtdb_taxonomy = NA)
          }

          # Update COG category choices - show individual categories with descriptions
          unique_cogs <- rv$gene_data %>%
            filter(!is.na(COG_category), COG_category != "-", nchar(COG_category) == 1) %>%
            distinct(COG_category, COG_category_name, COG_functional_group) %>%
            arrange(COG_functional_group, COG_category)

          cat("Found", nrow(unique_cogs), "individual COG categories\n")

          cog_choices <- setNames(
            c("all", unique_cogs$COG_category),
            c("All Categories", paste0(unique_cogs$COG_category, ": ", unique_cogs$COG_category_name))
          )

          updateSelectInput(session, "selected_cog", choices = cog_choices)

          # Update KEGG pathway choices
          kegg_pathways <- rv$gene_data %>%
            filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
            pull(KEGG_Pathway) %>%
            unique() %>%
            sort()

          updateSelectInput(session, "selected_kegg",
                           choices = c("All" = "all", kegg_pathways))

          # Join sample metadata (metadata dropdown will update automatically via observe())
          setProgress(0.7, detail = "Adding sample metadata...")
          metadata <- sample_metadata()

          # Check if metadata columns already exist in the data (e.g., disease_group, country, etc.)
          # These are the expected metadata columns from sample metadata
          expected_metadata_cols <- c("disease_group", "country", "continent", "host_sex", "age_group")
          has_metadata <- any(expected_metadata_cols %in% names(rv$gene_data)) &&
                          any(sapply(expected_metadata_cols, function(col) {
                            col %in% names(rv$gene_data) && sum(!is.na(rv$gene_data[[col]])) > 0
                          }))

          if (has_metadata) {
            cat("Metadata columns already present in data - skipping metadata join\n")
            existing_meta <- intersect(expected_metadata_cols, names(rv$gene_data))
            cat("Existing metadata columns:", paste(existing_meta, collapse = ", "), "\n")
          } else if (!is.null(metadata) && nrow(metadata) > 0) {
            tryCatch({
              cat("Processing sample metadata...\n")
              cat("Metadata rows:", nrow(metadata), "\n")
              cat("Metadata columns:", paste(names(metadata), collapse = ", "), "\n")

              # Check for sample_id column
              if ("sample_id" %in% names(rv$gene_data) && "accession_used_in_analysis" %in% names(metadata)) {
                # Join metadata to gene data on sample_id
                cat("Joining metadata on sample_id...\n")
                cat("Before join - gene data columns:", length(names(rv$gene_data)), "\n")

                rv$gene_data <- rv$gene_data %>%
                  left_join(metadata, by = c("sample_id" = "accession_used_in_analysis"))

                cat("After join - gene data columns:", length(names(rv$gene_data)), "\n")
                cat("New columns added:", paste(setdiff(names(rv$gene_data), c("sample_id", "clean_genome", "genome", "scaffold", "gene", "coverage", "breadth", "breadth_minCov", "SNV_N_count", "SNV_S_count", "N_sites", "S_sites", "pN", "pS", "pNpS", "COG_category", "Description", "Preferred_name", "KEGG_ko", "KEGG_Pathway", "COG_category_name", "COG_simplified", "domain", "phylum", "class", "order", "family", "genus", "species", "gtdb_taxonomy")), collapse = ", "), "\n")

                # Check if metadata columns were added
                if ("disease_group" %in% names(rv$gene_data)) {
                  n_with_metadata <- sum(!is.na(rv$gene_data$disease_group))
                  cat("Successfully added metadata! Genes with disease_group:", n_with_metadata, "\n")
                } else {
                  cat("WARNING: Metadata columns not found after join!\n")
                }
              } else {
                cat("WARNING: Cannot join metadata - missing join columns\n")
                cat("Gene data has sample_id?", "sample_id" %in% names(rv$gene_data), "\n")
                cat("Metadata has accession_used_in_analysis?", "accession_used_in_analysis" %in% names(metadata), "\n")
              }
            }, error = function(e) {
              cat("ERROR in metadata join:", e$message, "\n")
              showNotification(
                paste("Metadata join error:", e$message),
                type = "warning",
                duration = 5
              )
            })
          } else {
            cat("Sample metadata not available\n")
          }

          # Set data_loaded flag (metadata dropdown will update automatically via observe())
          rv$data_loaded <- TRUE
          cat("\nData loading complete. Metadata dropdown will update automatically.\n")

          # Summary notification with taxonomy status
          n_with_tax <- sum(!is.na(rv$gene_data$phylum))
          total_genes <- nrow(rv$gene_data)
          tax_percent <- round(100 * n_with_tax / total_genes, 1)

          if (n_with_tax == 0) {
            showNotification(
              paste0("⚠️ Loaded ", format(total_genes, big.mark = ","), " genes but NO TAXONOMY!\n",
                     "Analysis by taxonomic level will NOT work. ",
                     "Make sure GTDB data was loaded in Master Controller."),
              type = "error",
              duration = 15
            )
          } else if (tax_percent < 50) {
            showNotification(
              paste0("⚠️ Loaded ", format(total_genes, big.mark = ","), " genes; ",
                     format(n_with_tax, big.mark = ","), " with taxonomy (", tax_percent, "%)\n",
                     "Low taxonomy coverage may affect analysis."),
              type = "warning",
              duration = 10
            )
          } else {
            showNotification(
              paste0("✓ Loaded ", format(total_genes, big.mark = ","), " genes; ",
                     format(n_with_tax, big.mark = ","), " with taxonomy (", tax_percent, "%)"),
              type = "message",
              duration = 5
            )
          }

        }, error = function(e) {
          showNotification(
            paste("Error loading data:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    })

    # Update COG choices when display mode changes
    observeEvent(input$cog_display, {
      req(rv$data_loaded, rv$gene_data)

      if (input$cog_display == "grouped" && "COG_functional_group" %in% names(rv$gene_data)) {
        # Show functional groups
        cog_choices <- c("All" = "all", sort(unique(rv$gene_data$COG_functional_group)))
        updateSelectInput(session, "selected_cog", choices = cog_choices, selected = "all")

      } else {
        # Show individual categories (default)
        unique_cogs <- rv$gene_data %>%
          filter(!is.na(COG_category), COG_category != "-", nchar(COG_category) == 1) %>%
          distinct(COG_category, COG_category_name, COG_functional_group) %>%
          arrange(COG_functional_group, COG_category)

        cog_choices <- setNames(
          c("all", unique_cogs$COG_category),
          c("All Categories", paste0(unique_cogs$COG_category, ": ", unique_cogs$COG_category_name))
        )

        updateSelectInput(session, "selected_cog", choices = cog_choices, selected = "all")
      }
    })

    # Filter data
    observeEvent(input$run_analysis, {
      req(rv$gene_data)
      req(input$n_top_taxa)  # Ensure dynamic slider is rendered
      req(input$analysis_type)
      req(input$cog_display)
      req(input$tax_level)

      withProgress(message = "Filtering data...", {
        data <- rv$gene_data

        cat("\n=== Starting data filtering ===\n")
        cat("Initial rows:", nrow(data), "\n")

        # Apply coverage and breadth filters
        data <- data %>%
          filter(
            coverage >= input$min_coverage,
            breadth_minCov >= input$min_breadth,
            !is.na(pNpS),
            is.finite(pNpS)
          )

        cat("After quality filters:", nrow(data), "\n")

        # Determine functional category column
        if (input$analysis_type == "cog") {
          # Decide which column to use for grouping
          if (input$cog_display == "grouped" && "COG_functional_group" %in% names(data)) {
            func_col <- "COG_functional_group"
            cat("Using functional groups\n")
          } else {
            # Individual categories (default) - ensure only single-letter codes
            func_col <- "COG_category"
            cat("Using individual COG categories\n")

            # Filter to only single-letter COG categories
            data <- data %>% filter(nchar(COG_category) == 1)
            cat("After filtering to single-letter COG codes:", nrow(data), "\n")
          }

          # Remove NA/missing COG categories
          data <- data %>% filter(!is.na(.data[[func_col]]), .data[[func_col]] != "", .data[[func_col]] != "-")
          cat("After removing NA COG categories:", nrow(data), "\n")

          # Create combined labels (Code: Description) for better Y-axis labels
          if (input$cog_display == "individual" && "COG_category_name" %in% names(data)) {
            data <- data %>%
              mutate(COG_label = paste0(.data[[func_col]], ": ", COG_category_name))
            func_col <- "COG_label"
            cat("Created individual COG labels with descriptions\n")
            cat("Unique COG labels:", length(unique(data$COG_label)), "\n")
            cat("Sample labels:", paste(head(sort(unique(data$COG_label)), 10), collapse = "; "), "\n")
          }

          # Apply COG filter
          if (!"all" %in% input$selected_cog && !is.null(input$selected_cog) && length(input$selected_cog) > 0) {
            cat("Filtering by selected categories:", paste(input$selected_cog, collapse = ", "), "\n")
            # For individual with labels, filter by original category code
            if (func_col == "COG_label") {
              data <- data %>% filter(COG_category %in% input$selected_cog)
            } else if (func_col == "COG_functional_group") {
              data <- data %>% filter(COG_functional_group %in% input$selected_cog)
            } else {
              data <- data %>% filter(.data[[func_col]] %in% input$selected_cog)
            }
            cat("After COG filter:", nrow(data), "\n")
          } else {
            cat("No COG filter applied (showing all)\n")
          }
        } else {
          # KEGG analysis
          func_col <- "KEGG_Pathway"
          cat("Using KEGG pathway column\n")

          # Filter out missing KEGG data first
          data <- data %>% filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-", KEGG_Pathway != "")
          cat("After removing missing KEGG:", nrow(data), "\n")

          # Apply KEGG filter
          if (!"all" %in% input$selected_kegg && !is.null(input$selected_kegg) && length(input$selected_kegg) > 0) {
            cat("Filtering by KEGG pathways:", paste(input$selected_kegg, collapse = ", "), "\n")
            data <- data %>% filter(.data[[func_col]] %in% input$selected_kegg)
            cat("After KEGG filter:", nrow(data), "\n")
          } else {
            cat("No KEGG filter applied (all selected)\n")
          }
        }

        # Get taxonomic level
        tax_col <- input$tax_level
        cat("Taxonomic level:", tax_col, "\n")

        # Remove NA taxonomic assignments
        data <- data %>% filter(!is.na(.data[[tax_col]]))
        cat("After removing NA taxonomy:", nrow(data), "\n")

        # Check if any data remains
        if (nrow(data) == 0) {
          showNotification(
            paste("No genes with", tax_col, "taxonomy information!",
                  "Make sure to load GTDB data first in the Master Controller."),
            type = "error",
            duration = 10
          )
          rv$filtered_data <- data.frame()
          rv$raw_data <- data.frame()
          return(NULL)
        }

        # Optional: Limit to top N functional categories if slider is set below total
        # (Only apply limit if user has reduced the slider - otherwise show all)
        all_categories <- unique(data[[func_col]])
        n_total_categories <- length(all_categories)

        cat("\n=== Category Limiting Check ===\n")
        cat("Total unique categories in data:", n_total_categories, "\n")
        cat("Slider value (n_top_categories):", input$n_top_categories, "\n")
        cat("Will limit?", input$n_top_categories < n_total_categories, "\n")

        # DISABLED CATEGORY LIMITING - ALWAYS SHOW ALL CATEGORIES
        # This forces all categories to display regardless of slider value
        cat("**FORCING ALL CATEGORIES TO SHOW (limiting disabled)**\n")
        cat("Showing all", n_total_categories, "functional categories\n")

        # COMMENTED OUT THE LIMITING LOGIC:
        # if (("all" %in% input$selected_cog || "all" %in% input$selected_kegg ||
        #     is.null(input$selected_cog) || is.null(input$selected_kegg)) &&
        #     input$n_top_categories < n_total_categories) {
        #   cat("User selected top", input$n_top_categories, "out of", n_total_categories, "categories\n")
        #   top_categories <- data %>%
        #     group_by(.data[[func_col]]) %>%
        #     summarise(total_genes = n(), .groups = "drop") %>%
        #     arrange(desc(total_genes)) %>%
        #     head(input$n_top_categories) %>%
        #     pull(.data[[func_col]])
        #   data <- data %>% filter(.data[[func_col]] %in% top_categories)
        #   cat("After limiting categories:", nrow(data), "genes in", length(top_categories), "categories\n")
        #   cat("Limited categories:", paste(head(top_categories, 10), collapse = ", "), "\n")
        # }

        # Calculate mean pN/pS by functional category, taxonomy, and metadata (if selected)
        metadata_var <- input$metadata_var
        setProgress(0.5, detail = "Calculating means by taxonomy...")
        cat("Calculating means for", length(unique(data[[func_col]])), "functional categories and",
            length(unique(data[[tax_col]])), "taxa\n")

        # Add metadata variable to grouping if selected
        if (metadata_var != "none" && metadata_var %in% names(data)) {
          cat("Including metadata variable:", metadata_var, "\n")
          # Remove NA metadata values
          data <- data %>% filter(!is.na(.data[[metadata_var]]))
          cat("After removing NA metadata:", nrow(data), "\n")

          summary_data <- data %>%
            group_by(.data[[func_col]], .data[[tax_col]], .data[[metadata_var]]) %>%
            summarise(
              mean_pNpS = mean(pNpS, na.rm = TRUE),
              median_pNpS = median(pNpS, na.rm = TRUE),
              n_genes = n(),
              .groups = "drop"
            )
        } else {
          summary_data <- data %>%
            group_by(.data[[func_col]], .data[[tax_col]]) %>%
            summarise(
              mean_pNpS = mean(pNpS, na.rm = TRUE),
              median_pNpS = median(pNpS, na.rm = TRUE),
              n_genes = n(),
              .groups = "drop"
            )
        }

        # Identify top N taxa with user-selected sorting
        taxa_stats <- data %>%
          group_by(.data[[tax_col]]) %>%
          summarise(
            total_genes = n(),
            mean_pnps = mean(pNpS, na.rm = TRUE),
            .groups = "drop"
          )

        # Sort based on user selection
        # Note: For Visualization tab, "sig_count" falls back to "gene_count"
        # since p-values are only computed in Statistical Analysis tab
        sort_by <- input$taxa_sort_by
        if (is.null(sort_by) || sort_by == "sig_count") sort_by <- "gene_count"

        if (sort_by == "gene_count") {
          taxa_stats <- taxa_stats %>% arrange(desc(total_genes))
        } else {
          taxa_stats <- taxa_stats %>% arrange(.data[[tax_col]])
        }

        top_taxa <- head(taxa_stats[[tax_col]], input$n_top_taxa)

        # Create "Other" category for remaining taxa
        # Remove GTDB prefix (e.g., "p__", "c__", "o__", etc.) for cleaner display
        summary_data <- summary_data %>%
          mutate(
            taxa_display = if_else(.data[[tax_col]] %in% top_taxa,
                                   gsub("^[a-z]__", "", .data[[tax_col]]),  # Remove GTDB prefix
                                   "Other")
          )

        cat("\n=== Before 'Other' Re-aggregation ===\n")
        cat("Unique categories BEFORE re-aggregation:", length(unique(summary_data[[func_col]])), "\n")
        cat("Categories:", paste(head(sort(unique(summary_data[[func_col]])), 10), collapse = ", "), "...\n")
        cat("Total rows:", nrow(summary_data), "\n")

        # Recalculate means for "Other" group (preserve metadata grouping if selected)
        if (metadata_var != "none" && metadata_var %in% names(summary_data)) {
          # With metadata: group by func_col, taxa_display, AND metadata_var
          summary_data <- summary_data %>%
            group_by(.data[[func_col]], taxa_display, .data[[metadata_var]]) %>%
            summarise(
              mean_pNpS = mean(mean_pNpS, na.rm = TRUE),
              n_genes = sum(n_genes),
              .groups = "drop"
            )
        } else {
          # Without metadata: group by func_col and taxa_display only
          summary_data <- summary_data %>%
            group_by(.data[[func_col]], taxa_display) %>%
            summarise(
              mean_pNpS = mean(mean_pNpS, na.rm = TRUE),
              n_genes = sum(n_genes),
              .groups = "drop"
            )
        }

        cat("\n=== After 'Other' Re-aggregation ===\n")
        cat("Unique categories AFTER re-aggregation:", length(unique(summary_data[[func_col]])), "\n")
        cat("Categories:", paste(head(sort(unique(summary_data[[func_col]])), 10), collapse = ", "), "...\n")
        cat("Total rows:", nrow(summary_data), "\n")

        # Apply transformation
        summary_data <- summary_data %>%
          mutate(
            pNpS_transformed = if (input$transformation == "log") {
              log10(mean_pNpS + 1e-10)
            } else {
              mean_pNpS
            }
          )

        # Store functional column name for plotting
        summary_data$func_col_name <- func_col
        summary_data$tax_col_name <- tax_col

        cat("\n=== Summary ===\n")
        cat("Summary data rows:", nrow(summary_data), "\n")
        cat("Functional categories:", length(unique(summary_data[[func_col]])), "\n")
        cat("Taxa displayed:", length(unique(summary_data$taxa_display)), "\n")
        cat("Column names:", paste(names(summary_data), collapse = ", "), "\n")

        rv$filtered_data <- summary_data
        rv$raw_data <- data  # Store raw data for detailed table

        showNotification(
          paste("Analyzed", nrow(data), "genes across",
                length(unique(summary_data$taxa_display)), "taxa"),
          type = "message"
        )
      })
    })

    # Visualization
    output$gene_plot <- renderPlotly({
      req(rv$filtered_data)

      tryCatch({
        data <- rv$filtered_data

        # Validate data
        if (nrow(data) == 0) {
          return(plotly::plot_ly() %>%
            plotly::layout(
              title = "No data to display",
              xaxis = list(title = ""),
              yaxis = list(title = ""),
              annotations = list(
                text = "No data passed the filters. Try adjusting filter settings or selecting different categories.",
                xref = "paper", yref = "paper",
                x = 0.5, y = 0.5, showarrow = FALSE,
                font = list(size = 14, color = "red")
              )
            ))
        }

        # Get functional column name
        func_col <- unique(data$func_col_name)[1]
        if (is.null(func_col) || is.na(func_col)) {
          stop("Functional column name not found in data")
        }

        tax_col <- input$tax_level
        metadata_var <- input$metadata_var
        y_label <- if (input$transformation == "log") "Mean pN/pS (log10)" else "Mean pN/pS"

        # Create title based on analysis type and metadata
        title <- if (input$analysis_type == "cog") {
          base_title <- paste0("Mean Gene pN/pS by COG Category and ", tools::toTitleCase(tax_col))
          if (metadata_var != "none") paste0(base_title, " (grouped by ", metadata_var, ")") else base_title
        } else {
          base_title <- paste0("Mean Gene pN/pS by KEGG Pathway and ", tools::toTitleCase(tax_col))
          if (metadata_var != "none") paste0(base_title, " (grouped by ", metadata_var, ")") else base_title
        }

        # Calculate number of categories and taxa for plot sizing
        n_categories <- length(unique(data[[func_col]]))
        n_taxa <- length(unique(data$taxa_display))

        # ========================================================================
        # VISUALIZATION TYPE: HEATMAP or BAR CHART
        # ========================================================================

        if (input$viz_type == "heatmap") {
          # ===== HEATMAP VISUALIZATION =====

          # Prepare data for heatmap
          # For heatmap with metadata, create separate heatmaps per metadata group
          if (metadata_var != "none" && metadata_var %in% names(data)) {
            # With metadata: create faceted heatmaps
            n_groups <- length(unique(data[[metadata_var]]))

            # Order functional categories
            # For individual COG categories, order by functional group then alphabetically
            if (input$cog_display == "individual" && "COG_category" %in% names(data)) {
              # Extract COG_category from COG_label if needed
              if (func_col == "COG_label") {
                data <- data %>%
                  mutate(COG_category_for_sort = substr(COG_label, 1, 1))
              } else {
                data <- data %>%
                  mutate(COG_category_for_sort = COG_category)
              }

              # Order by functional group, then by mean pN/pS within group
              func_order <- data %>%
                group_by(.data[[func_col]], COG_functional_group, COG_category_for_sort) %>%
                summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
                arrange(COG_functional_group, desc(avg_pnps)) %>%
                pull(.data[[func_col]])
            } else {
              # For functional groups, order by mean pN/pS
              func_order <- data %>%
                group_by(.data[[func_col]]) %>%
                summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
                arrange(desc(avg_pnps)) %>%
                pull(.data[[func_col]])
            }

            # Order taxa by overall mean pN/pS (descending), "Other" last
            taxa_order <- data %>%
              group_by(taxa_display) %>%
              summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              arrange(desc(avg_pnps)) %>%
              pull(taxa_display)

            if ("Other" %in% taxa_order) {
              taxa_order <- c(setdiff(taxa_order, "Other"), "Other")
            }

            # Calculate mean pN/pS for each CATEGORY across all taxa (separate column on right)
            # This creates a "Mean" column showing category averages
            category_means <- data %>%
              group_by(.data[[func_col]], .data[[metadata_var]]) %>%
              summarise(mean_across_taxa = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              mutate(taxa_display = "Mean")  # Special taxa for category means

            # Add category means to data for visualization
            data_with_mean <- bind_rows(
              data,
              category_means %>% rename(pNpS_transformed = mean_across_taxa)
            )

            # Set factor levels for ordering
            data_with_mean[[func_col]] <- factor(data_with_mean[[func_col]],
                                                  levels = rev(func_order))  # Categories descending by mean

            # Update taxa_order to include "Mean" on the RIGHT (last position)
            taxa_order_with_mean <- c(taxa_order, "Mean")
            data_with_mean$taxa_display <- factor(data_with_mean$taxa_display, levels = taxa_order_with_mean)

            cat("Added Mean column with", nrow(category_means), "category averages\n")

            p <- ggplot(data_with_mean, aes_string(x = "taxa_display", y = func_col, fill = "pNpS_transformed")) +
              geom_tile(color = "white", size = 0.5) +
              facet_wrap(as.formula(paste("~", metadata_var)),
                        scales = "free_x",  # Only X-axis free, Y-axis shared (labels only on left)
                        ncol = min(n_groups, 2)) +
              scale_fill_gradient2(
                low = "blue", mid = "white", high = "red",  # coolwarm colormap (like seaborn)
                midpoint = if (input$transformation == "log") 0 else 1,  # log10(1)=0, or raw 1.0
                name = y_label,
                limits = c(min(data_with_mean$pNpS_transformed, na.rm = TRUE),
                          max(data_with_mean$pNpS_transformed, na.rm = TRUE))
              ) +
              labs(
                title = title,
                x = tools::toTitleCase(tax_col),
                y = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway"
              ) +
              theme_minimal(base_size = 13) +
              theme(
                # Ensure consistent text sizes across all facets
                axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
                axis.text.y = element_text(size = 10),
                axis.title = element_text(size = 12, face = "bold"),
                legend.position = "right",
                legend.title = element_text(size = 11, face = "bold"),
                legend.text = element_text(size = 10),
                strip.text = element_text(size = 12, face = "bold", color = "white"),
                strip.background = element_rect(fill = "#2C3E50", color = NA),
                panel.spacing = unit(1.5, "lines"),
                plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                plot.margin = margin(10, 10, 20, 10),
                # Prevent text scaling based on facet size
                strip.clip = "off"
              ) +
              # Force consistent axis text size across facets
              coord_cartesian(clip = "off")

            plot_height <- max(800, n_categories * 35 + 250)

          } else {
            # Without metadata: single heatmap

            # Order functional categories
            # For individual COG categories, order by functional group then alphabetically
            if (input$cog_display == "individual" && "COG_category" %in% names(data)) {
              # Extract COG_category from COG_label if needed
              if (func_col == "COG_label") {
                data <- data %>%
                  mutate(COG_category_for_sort = substr(COG_label, 1, 1))
              } else {
                data <- data %>%
                  mutate(COG_category_for_sort = COG_category)
              }

              # Order by functional group, then by mean pN/pS within group
              func_order <- data %>%
                group_by(.data[[func_col]], COG_functional_group, COG_category_for_sort) %>%
                summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
                arrange(COG_functional_group, desc(avg_pnps)) %>%
                pull(.data[[func_col]])
            } else {
              # For functional groups, order by mean pN/pS
              func_order <- data %>%
                group_by(.data[[func_col]]) %>%
                summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
                arrange(desc(avg_pnps)) %>%
                pull(.data[[func_col]])
            }

            # Order taxa by mean pN/pS (descending), "Other" last
            taxa_order <- data %>%
              group_by(taxa_display) %>%
              summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              arrange(desc(avg_pnps)) %>%
              pull(taxa_display)

            if ("Other" %in% taxa_order) {
              taxa_order <- c(setdiff(taxa_order, "Other"), "Other")
            }

            # Calculate mean pN/pS for each CATEGORY across all taxa (separate column on right)
            category_means <- data %>%
              group_by(.data[[func_col]]) %>%
              summarise(mean_across_taxa = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              mutate(taxa_display = "Mean")  # Special taxa for category means

            # Add category means to data for visualization
            data_with_mean <- bind_rows(
              data,
              category_means %>% rename(pNpS_transformed = mean_across_taxa)
            )

            # Set factor levels for ordering
            data_with_mean[[func_col]] <- factor(data_with_mean[[func_col]],
                                                  levels = rev(func_order))  # Categories descending by mean

            # Update taxa_order to include "Mean" on the RIGHT (last position)
            taxa_order_with_mean <- c(taxa_order, "Mean")
            data_with_mean$taxa_display <- factor(data_with_mean$taxa_display, levels = taxa_order_with_mean)

            cat("Added Mean column with", nrow(category_means), "category averages\n")

            p <- ggplot(data_with_mean, aes_string(x = "taxa_display", y = func_col, fill = "pNpS_transformed")) +
              geom_tile(color = "white", size = 0.5) +
              scale_fill_gradient2(
                low = "blue", mid = "white", high = "red",  # coolwarm colormap (like seaborn)
                midpoint = if (input$transformation == "log") 0 else 1,  # log10(1)=0, or raw 1.0
                name = y_label,
                limits = c(min(data_with_mean$pNpS_transformed, na.rm = TRUE),
                          max(data_with_mean$pNpS_transformed, na.rm = TRUE))
              ) +
              labs(
                title = title,
                x = tools::toTitleCase(tax_col),
                y = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway"
              ) +
              theme_minimal(base_size = 14) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
                axis.text.y = element_text(size = 11, hjust = 1),
                axis.title = element_text(size = 13, face = "bold"),
                legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10),
                legend.key.height = unit(1.2, "cm"),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.margin = margin(15, 15, 15, 15),
                panel.grid = element_blank()
              )

            plot_height <- max(900, n_categories * 40 + 300)
          }

          cat("Heatmap plot dimensions: height =", plot_height, "px, categories =", n_categories, ", taxa =", n_taxa, "\n")

          ggplotly(p, tooltip = c("x", "y", "fill")) %>%
            plotly::layout(
              height = plot_height,
              margin = list(l = 180, r = 200, t = 120, b = 150),  # Increased for more taxa/categories
              xaxis = list(tickfont = list(size = 11)),
              yaxis = list(tickfont = list(size = 10))
            )

        } else if (input$viz_type == "bar") {
          # ===== BAR CHART VISUALIZATION =====

          # Order taxa: "Other" last, rest by mean pN/pS
          taxa_in_data <- unique(data$taxa_display)
          taxa_order <- data %>%
            group_by(taxa_display) %>%
            summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(avg_pnps)) %>%
            pull(taxa_display)

          # Put "Other" at the end if it exists
          if ("Other" %in% taxa_order) {
            taxa_order <- c(setdiff(taxa_order, "Other"), "Other")
          }

          data$taxa_display <- factor(data$taxa_display, levels = taxa_order)

          # Create grouped bar chart with optional metadata faceting
          if (metadata_var != "none" && metadata_var %in% names(data)) {
          # With metadata: use facets to show metadata groups
          # Count number of metadata groups to adjust layout
          n_groups <- length(unique(data[[metadata_var]]))

          p <- ggplot(data, aes_string(x = func_col, y = "pNpS_transformed", fill = "taxa_display")) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.8) +
            facet_wrap(as.formula(paste("~", metadata_var)),
                      scales = "free_x",
                      ncol = min(n_groups, 2)) +  # Max 2 columns for more width per panel
            labs(
              title = title,
              x = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway",
              y = y_label,
              fill = tools::toTitleCase(tax_col)
            ) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 15)),
              axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 15)),
              legend.position = "right",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 11, face = "bold"),
              legend.key.size = unit(0.8, "cm"),
              strip.text = element_text(size = 12, face = "bold", color = "white"),
              strip.background = element_rect(fill = "#2C3E50", color = NA),
              panel.spacing = unit(1.5, "lines"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.margin = margin(10, 10, 20, 10)
            )
        } else {
          # Without metadata: standard grouped bar chart
          p <- ggplot(data, aes_string(x = func_col, y = "pNpS_transformed", fill = "taxa_display")) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.8) +
            labs(
              title = title,
              x = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway",
              y = y_label,
              fill = tools::toTitleCase(tax_col)
            ) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
              axis.text.y = element_text(size = 11),
              axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 15)),
              axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 15)),
              legend.position = "right",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12, face = "bold"),
              legend.key.size = unit(0.8, "cm"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.margin = margin(10, 10, 20, 10)
            )
        }

        # Use distinct colors for taxa
        n_taxa <- length(unique(data$taxa_display))
        if (n_taxa <= 10) {
          p <- p + scale_fill_brewer(palette = "Paired")
        } else if (n_taxa <= 20) {
          p <- p + scale_fill_viridis_d(option = "turbo", end = 0.95)
        } else {
          p <- p + scale_fill_manual(values = rainbow(n_taxa))
        }

        # Adjust height based on whether facets are used and amount of data
        plot_height <- if (metadata_var != "none" && metadata_var %in% names(data)) {
          # Faceted plot: increase height based on number of groups and categories
          n_groups <- length(unique(data[[metadata_var]]))
          n_rows <- ceiling(n_groups / 2)  # Using 2 columns max
          base_height <- 600 + (n_categories * 30)  # More height for more categories
          max(base_height, n_rows * 500)  # At least 500px per row
        } else {
          # Non-faceted: base height + extra for number of categories
          800 + (n_categories * 35)  # Scale with number of categories
        }

        cat("Plot dimensions: height =", plot_height, "px, categories =", n_categories, ", taxa =", n_taxa, "\n")

        ggplotly(p, tooltip = c("x", "y", "fill")) %>%
          plotly::layout(
            barmode = "group",
            legend = list(
              orientation = "v",
              x = 1.02,
              y = 0.5,
              font = list(size = 11),
              tracegroupgap = 0
            ),
            height = plot_height,
            margin = list(l = 80, r = 150, t = 100, b = 150),  # More margin for labels
            xaxis = list(tickfont = list(size = 11)),
            yaxis = list(tickfont = list(size = 11))
          )

        }  # End of bar chart block

      }, error = function(e) {
        cat("ERROR in gene_plot:", e$message, "\n")
        return(plotly::plot_ly() %>%
          plotly::layout(
            title = "Visualization Error",
            xaxis = list(title = ""),
            yaxis = list(title = ""),
            annotations = list(
              text = paste("Error:", e$message, "\nCheck console for details."),
              xref = "paper", yref = "paper",
              x = 0.5, y = 0.5, showarrow = FALSE,
              font = list(size = 12, color = "red")
            )
          ))
      })
    })

    # Summary statistics by functional category and taxonomy
    output$cog_summary_table <- DT::renderDataTable({
      req(rv$filtered_data)

      data <- rv$filtered_data
      func_col <- unique(data$func_col_name)[1]
      tax_col <- input$tax_level

      # Create summary table
      summary_stats <- data %>%
        select(
          Category = !!sym(func_col),
          Taxa = taxa_display,
          Mean_pNpS = mean_pNpS,
          N_genes = n_genes
        ) %>%
        arrange(Category, desc(Mean_pNpS))

      # Add category names for individual COG view
      if (input$analysis_type == "cog" && input$cog_display == "individual" && exists("rv$raw_data")) {
        if (!is.null(rv$raw_data)) {
          cog_names <- rv$raw_data %>%
            filter(nchar(COG_category) == 1) %>%
            distinct(COG_category, COG_category_name) %>%
            rename(Category = COG_category)

          summary_stats <- summary_stats %>%
            left_join(cog_names, by = "Category") %>%
            select(Category, COG_category_name, everything())
        }
      }

      DT::datatable(
        summary_stats,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          order = list(list(2, 'desc'))
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        caption = paste("Summary by",
                       if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway",
                       "and", tools::toTitleCase(tax_col))
      ) %>%
        DT::formatRound(columns = c("Mean_pNpS"), digits = 4)
    })

    # Gene details table
    output$gene_details_table <- DT::renderDataTable({
      req(rv$raw_data)

      tax_col <- input$tax_level

      # Select columns based on analysis type
      if (input$analysis_type == "cog") {
        gene_details <- rv$raw_data %>%
          select(
            gene,
            sample_id,
            clean_genome,
            !!sym(tax_col),
            pNpS,
            coverage,
            breadth_minCov,
            COG_category,
            COG_category_name,
            Description,
            Preferred_name
          ) %>%
          mutate(!!sym(tax_col) := gsub("^[a-z]__", "", .data[[tax_col]])) %>%  # Clean GTDB prefix
          arrange(desc(pNpS))
      } else {
        gene_details <- rv$raw_data %>%
          select(
            gene,
            sample_id,
            clean_genome,
            !!sym(tax_col),
            pNpS,
            coverage,
            breadth_minCov,
            KEGG_Pathway,
            KEGG_ko,
            Description,
            Preferred_name
          ) %>%
          mutate(!!sym(tax_col) := gsub("^[a-z]__", "", .data[[tax_col]])) %>%  # Clean GTDB prefix
          arrange(desc(pNpS))
      }

      DT::datatable(
        gene_details,
        options = list(
          pageLength = 50,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          order = list(list(4, 'desc'))  # Order by pNpS column
        ),
        rownames = FALSE,
        filter = "top",
        extensions = 'Buttons',
        caption = paste("Individual gene-level data (filtered by quality thresholds)")
      ) %>%
        DT::formatRound(columns = c("pNpS", "coverage", "breadth_minCov"), digits = 4)
    })

    # Data summary
    output$data_summary <- renderPrint({
      req(rv$raw_data, rv$filtered_data)

      raw_data <- rv$raw_data
      summary_data <- rv$filtered_data
      tax_col <- input$tax_level
      func_col <- unique(summary_data$func_col_name)[1]

      cat("Gene-Level pN/pS Data Summary\n")
      cat("================================\n\n")

      cat("Analysis Type:", if (input$analysis_type == "cog") "COG Categories" else "KEGG Pathways", "\n")
      cat("Taxonomic Level:", tools::toTitleCase(tax_col), "\n\n")

      cat("Raw Data (after quality filters):\n")
      cat("  Total genes:", nrow(raw_data), "\n")
      cat("  Unique samples:", length(unique(raw_data$sample_id)), "\n")
      cat("  Unique genomes:", length(unique(raw_data$clean_genome)), "\n")
      cat("  Unique", tax_col, ":", length(unique(raw_data[[tax_col]])), "\n\n")

      if (input$analysis_type == "cog") {
        cat("  Genes with COG annotation:",
            sum(!is.na(raw_data$COG_category) & raw_data$COG_category != "-"), "\n\n")
      } else {
        cat("  Genes with KEGG annotation:",
            sum(!is.na(raw_data$KEGG_Pathway) & raw_data$KEGG_Pathway != "-"), "\n\n")
      }

      cat("pN/pS Statistics (raw gene-level):\n")
      cat("  Min:", round(min(raw_data$pNpS, na.rm = TRUE), 4), "\n")
      cat("  Q1:", round(quantile(raw_data$pNpS, 0.25, na.rm = TRUE), 4), "\n")
      cat("  Median:", round(median(raw_data$pNpS, na.rm = TRUE), 4), "\n")
      cat("  Mean:", round(mean(raw_data$pNpS, na.rm = TRUE), 4), "\n")
      cat("  Q3:", round(quantile(raw_data$pNpS, 0.25, na.rm = TRUE), 4), "\n")
      cat("  Max:", round(max(raw_data$pNpS, na.rm = TRUE), 4), "\n\n")

      cat("Aggregated Summary:\n")
      cat("  Functional categories:", length(unique(summary_data[[func_col]])), "\n")
      cat("  Taxa displayed:", length(unique(summary_data$taxa_display)), "\n")
      cat("  (Top", input$n_top_taxa, "taxa + 'Other'):\n")

      # List top taxa
      top_taxa_list <- unique(summary_data$taxa_display)
      top_taxa_list <- top_taxa_list[top_taxa_list != "Other"]
      for (i in seq_along(top_taxa_list)) {
        cat("    ", i, ". ", top_taxa_list[i], "\n", sep = "")
      }
      if ("Other" %in% unique(summary_data$taxa_display)) {
        cat("    ... Other taxa grouped\n")
      }
    })

    # ===========================================================================
    # Statistical Significance Heatmap
    # ===========================================================================

    output$significance_heatmap <- renderPlotly({
      req(rv$raw_data)

      tryCatch({
        # Check metadata variable is selected
        metadata_var <- input$metadata_var
        if (is.null(metadata_var) || metadata_var == "none") {
          return(plotly::plot_ly() %>%
            plotly::layout(
              title = "No Metadata Selected",
              annotations = list(
                text = "Please select a metadata variable to perform statistical tests",
                xref = "paper", yref = "paper",
                x = 0.5, y = 0.5, showarrow = FALSE,
                font = list(size = 14, color = "red")
              )
            ))
        }

        cat("\n=== Computing Statistical Significance Matrix ===\n")

        data <- rv$raw_data
        tax_col <- input$tax_level

        # Get functional column
        if (input$analysis_type == "cog") {
          if (input$cog_display == "grouped" && "COG_functional_group" %in% names(data)) {
            func_col <- "COG_functional_group"
          } else {
            func_col <- "COG_category"
            # Filter to single-letter COG codes
            data <- data %>% filter(nchar(COG_category) == 1)
          }
        } else {
          func_col <- "KEGG_Pathway"
        }

        # Filter data
        data <- data %>%
          filter(
            !is.na(.data[[func_col]]),
            .data[[func_col]] != "-",
            .data[[func_col]] != "",
            !is.na(.data[[tax_col]]),
            !is.na(.data[[metadata_var]]),
            !is.na(pNpS),
            is.finite(pNpS)
          )

        cat("Data rows for testing:", nrow(data), "\n")

        # Determine metadata type
        meta_values <- unique(data[[metadata_var]])
        is_numeric <- is.numeric(data[[metadata_var]])
        n_groups <- length(meta_values)

        cat("Metadata variable:", metadata_var, "\n")
        cat("Metadata type:", if(is_numeric) "numeric" else "categorical", "\n")
        cat("Number of groups:", n_groups, "\n")

        # Get unique categories and ALL taxa first
        categories <- sort(unique(data[[func_col]]))
        all_taxa <- unique(data[[tax_col]])

        cat("Categories:", length(categories), "\n")
        cat("All taxa:", length(all_taxa), "\n")

        # Get number of cores from input
        n_cores <- input$n_cores
        if (is.null(n_cores)) n_cores <- 1
        n_cores <- min(n_cores, parallel::detectCores())

        cat("Using", n_cores, "CPU cores for parallel processing\n")

        # ============================================================
        # Step 1: Calculate p-values for ALL taxa (for significance-based sorting)
        # ============================================================
        combinations_all <- expand.grid(cat_name = categories, tax_name = all_taxa,
                                       stringsAsFactors = FALSE)
        n_combinations_all <- nrow(combinations_all)

        cat("Total combinations to test (all taxa):", n_combinations_all, "\n")

        # Define test function
        run_test <- function(i) {
          cat_name <- combinations_all$cat_name[i]
          tax_name <- combinations_all$tax_name[i]

          test_data <- data %>%
            filter(.data[[func_col]] == cat_name, .data[[tax_col]] == tax_name)

          if (nrow(test_data) < 3) {
            return(list(cat = cat_name, tax = tax_name, pval = NA))
          }

          pval <- tryCatch({
            if (is_numeric) {
              cor_test <- cor.test(test_data$pNpS, test_data[[metadata_var]],
                                  method = "spearman", exact = FALSE)
              cor_test$p.value
            } else if (n_groups == 2) {
              groups <- split(test_data$pNpS, test_data[[metadata_var]])
              if (length(groups) == 2 && all(sapply(groups, length) >= 2)) {
                wilcox.test(groups[[1]], groups[[2]])$p.value
              } else {
                NA
              }
            } else {
              kruskal.test(pNpS ~ .data[[metadata_var]], data = test_data)$p.value
            }
          }, error = function(e) NA)

          list(cat = cat_name, tax = tax_name, pval = pval)
        }

        # Run tests in parallel (or sequential if only 1 core)
        if (n_cores > 1 && .Platform$OS.type != "windows") {
          results_all <- parallel::mclapply(1:n_combinations_all, run_test, mc.cores = n_cores)
        } else {
          results_all <- lapply(1:n_combinations_all, run_test)
        }

        # Convert results to data frame for easier manipulation
        results_df <- do.call(rbind, lapply(results_all, function(x) {
          data.frame(cat = x$cat, tax = x$tax, pval = x$pval, stringsAsFactors = FALSE)
        }))

        cat("P-values computed for all taxa\n")

        # ============================================================
        # Step 2: Calculate significance count per taxon for sorting
        # ============================================================
        taxa_sig_count <- results_df %>%
          group_by(tax) %>%
          summarise(
            sig_count = sum(pval < 0.05, na.rm = TRUE),
            gene_count = n(),
            .groups = "drop"
          )

        # Join with gene count from original data
        taxa_gene_count <- data %>%
          group_by(.data[[tax_col]]) %>%
          summarise(total_genes = n(), .groups = "drop")

        taxa_stats <- taxa_sig_count %>%
          left_join(taxa_gene_count, by = c("tax" = tax_col))

        # Sort taxa based on user selection
        sort_by <- input$taxa_sort_by
        if (is.null(sort_by)) sort_by <- "sig_count"

        if (sort_by == "sig_count") {
          taxa_stats <- taxa_stats %>% arrange(desc(sig_count), desc(total_genes))
          cat("Sorting by significance count (most * first)\n")
        } else if (sort_by == "gene_count") {
          taxa_stats <- taxa_stats %>% arrange(desc(total_genes))
          cat("Sorting by gene count\n")
        } else {
          taxa_stats <- taxa_stats %>% arrange(tax)
          cat("Sorting alphabetically\n")
        }

        # Apply top N taxa filter
        n_top <- input$n_top_taxa
        if (is.null(n_top)) n_top <- 15

        taxa <- head(taxa_stats$tax, n_top)

        cat("Top", n_top, "taxa selected\n")
        cat("Taxa:", length(taxa), "\n")

        # Debug: Show taxa with their significance counts (already sorted)
        debug_df <- head(taxa_stats, n_top)
        cat("Taxa ordering (by sig_count):\n")
        for (i in 1:min(5, nrow(debug_df))) {
          cat(sprintf("  %d. %s: %d significant tests\n",
                      i, debug_df$tax[i], debug_df$sig_count[i]))
        }
        if (nrow(debug_df) > 5) cat("  ...\n")

        # ============================================================
        # Step 3: Filter results to top taxa only and create matrix
        # ============================================================
        results_filtered <- results_df %>% filter(tax %in% taxa)

        # Create matrix with proper ordering
        pvalue_matrix <- matrix(NA, nrow = length(categories), ncol = length(taxa),
                               dimnames = list(categories, taxa))

        for (i in 1:nrow(results_filtered)) {
          pvalue_matrix[results_filtered$cat[i], results_filtered$tax[i]] <- results_filtered$pval[i]
        }

        cat("P-value matrix created for top taxa\n")

        # Convert matrix to long format for plotting
        pvalue_df <- as.data.frame(pvalue_matrix) %>%
          tibble::rownames_to_column("Category") %>%
          tidyr::pivot_longer(-Category, names_to = "Taxon", values_to = "p_value")

        # Debug: Check pvalue_df structure before factor conversion
        cat("DEBUG pvalue_df: nrow=", nrow(pvalue_df), ", ncol=", ncol(pvalue_df), "\n")
        cat("DEBUG pvalue_df columns:", paste(names(pvalue_df), collapse=", "), "\n")
        cat("DEBUG unique Taxon values:", paste(head(unique(pvalue_df$Taxon), 5), collapse=", "), "\n")
        cat("DEBUG taxa vector:", paste(head(taxa, 5), collapse=", "), "\n")

        # IMPORTANT: Set Taxon as factor with proper ordering (by significance count)
        # Only use taxa that actually exist in pvalue_df
        valid_taxa <- intersect(taxa, unique(pvalue_df$Taxon))
        if (length(valid_taxa) == 0) {
          stop("No valid taxa found after filtering")
        }
        pvalue_df$Taxon <- factor(pvalue_df$Taxon, levels = valid_taxa)
        # Update taxa to match valid ones for scale_x_discrete
        taxa <- valid_taxa

        cat("Taxon order (by sig count):", paste(head(taxa, 5), collapse = ", "), "...\n")

        # Add significance stars
        pvalue_df <- pvalue_df %>%
          mutate(
            significance = case_when(
              is.na(p_value) ~ "",
              p_value < 0.001 ~ "***",
              p_value < 0.01 ~ "**",
              p_value < 0.05 ~ "*",
              TRUE ~ ""
            ),
            log_p = -log10(p_value + 1e-300),  # Transform for color scale
            display_text = case_when(
              is.na(p_value) ~ "N/A",
              p_value < 0.001 ~ paste0("***\np=", formatC(p_value, format = "e", digits = 2)),
              p_value < 0.01 ~ paste0("**\np=", formatC(p_value, format = "f", digits = 4)),
              p_value < 0.05 ~ paste0("*\np=", formatC(p_value, format = "f", digits = 3)),
              TRUE ~ paste0("ns\np=", formatC(p_value, format = "f", digits = 3))
            )
          )

        # Order categories by functional group if individual COG
        if (input$cog_display == "individual" && "COG_category" %in% names(data) &&
            "COG_functional_group" %in% names(data)) {
          # Get functional group ordering
          cat_order_df <- data %>%
            filter(nchar(COG_category) == 1) %>%
            distinct(COG_category, COG_functional_group) %>%
            arrange(COG_functional_group, COG_category)

          cat_order <- cat_order_df$COG_category
          # Ensure all categories in pvalue_df are in cat_order
          if (length(cat_order) > 0) {
            # Add any missing categories to the end
            missing_cats <- setdiff(unique(pvalue_df$Category), cat_order)
            cat_order <- c(cat_order, missing_cats)
            pvalue_df$Category <- factor(pvalue_df$Category, levels = cat_order)
          }
        } else {
          # Default: alphabetical order
          pvalue_df$Category <- factor(pvalue_df$Category, levels = sort(unique(pvalue_df$Category)))
        }

        # Create heatmap
        test_type <- if (is_numeric) {
          "Spearman Correlation"
        } else if (n_groups == 2) {
          "Wilcoxon Test"
        } else {
          "Kruskal-Wallis Test"
        }

        title <- paste0("Statistical Significance Matrix (",
                       test_type, " by ", metadata_var, ")")

        # Debug: Final check before ggplot
        cat("DEBUG before ggplot:\n")
        cat("  pvalue_df nrow:", nrow(pvalue_df), "\n")
        cat("  pvalue_df columns:", paste(names(pvalue_df), collapse=", "), "\n")
        cat("  taxa vector length:", length(taxa), "\n")
        cat("  taxa values:", paste(taxa, collapse=", "), "\n")
        cat("  Taxon levels:", paste(levels(pvalue_df$Taxon), collapse=", "), "\n")

        # Ensure data is not empty
        if (nrow(pvalue_df) == 0) {
          stop("No data available for plotting")
        }

        # Plot with plotly
        # Note: Factor levels for Taxon are already set, so ggplot will use them for ordering
        # Remove scale_x_discrete(limits=) to avoid "Must use existing variables" error
        p <- ggplot(pvalue_df, aes(x = Taxon, y = Category, fill = log_p)) +
          geom_tile(color = "white", size = 0.5) +
          geom_text(aes(label = significance), size = 5, fontface = "bold", color = "black") +
          scale_fill_gradient2(
            low = "blue", mid = "white", high = "red",
            midpoint = -log10(0.05),  # Midpoint at p=0.05
            name = "-log10(p)",
            na.value = "grey60"  # More visible grey for NA/insufficient data
          ) +
          scale_x_discrete(drop = FALSE) +  # Use factor levels, don't drop unused
          scale_y_discrete(drop = FALSE) +  # Same for categories
          labs(
            title = title,
            x = tools::toTitleCase(tax_col),
            y = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway"
          ) +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(size = 11, face = "bold"),
            legend.position = "right",
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid = element_blank()
          )

        plot_height <- max(600, length(categories) * 25 + 200)

        ggplotly(p, tooltip = "display_text") %>%
          plotly::layout(
            height = plot_height,
            margin = list(l = 150, r = 150, t = 100, b = 120)
          )

      }, error = function(e) {
        cat("ERROR in significance_heatmap:", e$message, "\n")
        return(plotly::plot_ly() %>%
          plotly::layout(
            title = "Statistical Testing Error",
            annotations = list(
              text = paste("Error:", e$message, "\nCheck console for details."),
              xref = "paper", yref = "paper",
              x = 0.5, y = 0.5, showarrow = FALSE,
              font = list(size = 12, color = "red")
            )
          ))
      })
    })

    # ===========================================================================
    # Taxon Comparison Plot (Average pN/pS by Metadata Group)
    # ===========================================================================

    output$taxon_comparison_plot <- renderPlotly({
      req(rv$raw_data)

      tryCatch({
        metadata_var <- input$metadata_var
        if (is.null(metadata_var) || metadata_var == "none") {
          return(plotly::plot_ly() %>%
            plotly::layout(
              title = "No Metadata Selected",
              annotations = list(
                text = "Please select a metadata variable to compare groups",
                xref = "paper", yref = "paper",
                x = 0.5, y = 0.5, showarrow = FALSE,
                font = list(size = 14, color = "red")
              )
            ))
        }

        cat("\n=== Creating Taxon Comparison Plot ===\n")

        data <- rv$raw_data
        tax_col <- input$tax_level

        # Filter data
        data <- data %>%
          filter(
            !is.na(.data[[tax_col]]),
            !is.na(.data[[metadata_var]]),
            !is.na(pNpS),
            is.finite(pNpS)
          )

        cat("Data rows:", nrow(data), "\n")

        # Get top taxa with user-selected sorting
        taxa_stats <- data %>%
          group_by(.data[[tax_col]]) %>%
          summarise(
            gene_count = n(),
            mean_pnps = mean(pNpS, na.rm = TRUE),
            .groups = "drop"
          )

        # Sort based on user selection
        # Note: "sig_count" falls back to "gene_count" here since p-values aren't available
        sort_by <- input$taxa_sort_by
        if (is.null(sort_by) || sort_by == "sig_count") sort_by <- "gene_count"

        if (sort_by == "gene_count") {
          taxa_stats <- taxa_stats %>% arrange(desc(gene_count))
        } else {
          taxa_stats <- taxa_stats %>% arrange(.data[[tax_col]])
        }

        top_taxa <- head(taxa_stats[[tax_col]], input$n_top_taxa)

        # Filter to top taxa
        data <- data %>%
          filter(.data[[tax_col]] %in% top_taxa)

        # Calculate mean and SE for each taxon + metadata group
        summary_data <- data %>%
          group_by(.data[[tax_col]], .data[[metadata_var]]) %>%
          summarise(
            mean_pnps = mean(pNpS, na.rm = TRUE),
            se_pnps = sd(pNpS, na.rm = TRUE) / sqrt(n()),
            n = n(),
            .groups = "drop"
          ) %>%
          mutate(
            se_pnps = ifelse(is.na(se_pnps), 0, se_pnps),
            ymin = mean_pnps - se_pnps,
            ymax = mean_pnps + se_pnps
          )

        # Order taxa by overall mean pN/pS
        taxa_order <- summary_data %>%
          group_by(.data[[tax_col]]) %>%
          summarise(overall_mean = mean(mean_pnps, na.rm = TRUE), .groups = "drop") %>%
          arrange(desc(overall_mean)) %>%
          pull(.data[[tax_col]])

        summary_data[[tax_col]] <- factor(summary_data[[tax_col]], levels = taxa_order)

        cat("Taxa in plot:", length(taxa_order), "\n")
        cat("Groups:", paste(unique(summary_data[[metadata_var]]), collapse = ", "), "\n")

        # Create grouped bar chart
        p <- ggplot(summary_data, aes(x = .data[[tax_col]], y = mean_pnps,
                                      fill = .data[[metadata_var]])) +
          geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
          geom_errorbar(aes(ymin = ymin, ymax = ymax),
                       position = position_dodge(width = 0.8), width = 0.3, color = "black") +
          labs(
            title = paste("Average pN/pS by", tools::toTitleCase(tax_col), "and", metadata_var),
            x = tools::toTitleCase(tax_col),
            y = "Mean pN/pS (± SE)",
            fill = metadata_var
          ) +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 11, face = "bold"),
            legend.position = "top",
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
          ) +
          scale_fill_brewer(palette = "Set1")

        ggplotly(p, tooltip = c("x", "y", "fill")) %>%
          plotly::layout(
            margin = list(l = 80, r = 50, t = 80, b = 150)
          )

      }, error = function(e) {
        cat("ERROR in taxon_comparison_plot:", e$message, "\n")
        return(plotly::plot_ly() %>%
          plotly::layout(
            title = "Error Creating Plot",
            annotations = list(
              text = paste("Error:", e$message),
              xref = "paper", yref = "paper",
              x = 0.5, y = 0.5, showarrow = FALSE,
              font = list(size = 12, color = "red")
            )
          ))
      })
    })

    # Make significance heatmap calculate immediately on Update Analysis
    # (not lazily when tab is viewed)
    outputOptions(output, "significance_heatmap", suspendWhenHidden = FALSE)

    # ========== Download Modal Handlers ==========

    # Gene Plot ggplot - matches main visualization logic exactly
    gene_plot_ggplot <- reactive({
      req(rv$filtered_data)
      data <- rv$filtered_data
      req(nrow(data) > 0)

      func_col <- unique(data$func_col_name)[1]
      if (is.null(func_col) || is.na(func_col)) {
        stop("Functional column name not found in data")
      }

      tax_col <- input$tax_level
      metadata_var <- input$metadata_var
      y_label <- if (input$transformation == "log") "Mean pN/pS (log10)" else "Mean pN/pS"

      # Create title based on analysis type and metadata
      title <- if (input$analysis_type == "cog") {
        base_title <- paste0("Mean Gene pN/pS by COG Category and ", tools::toTitleCase(tax_col))
        if (metadata_var != "none") paste0(base_title, " (grouped by ", metadata_var, ")") else base_title
      } else {
        base_title <- paste0("Mean Gene pN/pS by KEGG Pathway and ", tools::toTitleCase(tax_col))
        if (metadata_var != "none") paste0(base_title, " (grouped by ", metadata_var, ")") else base_title
      }

      # Calculate number of categories and taxa
      n_categories <- length(unique(data[[func_col]]))
      n_taxa <- length(unique(data$taxa_display))

      # ========================================================================
      # VISUALIZATION TYPE: HEATMAP or BAR CHART (matching main plot logic)
      # ========================================================================

      if (input$viz_type == "heatmap") {
        # ===== HEATMAP VISUALIZATION =====

        if (metadata_var != "none" && metadata_var %in% names(data)) {
          # With metadata: create faceted heatmaps
          n_groups <- length(unique(data[[metadata_var]]))

          # Order functional categories
          if (input$cog_display == "individual" && "COG_category" %in% names(data)) {
            # Extract COG_category from COG_label if needed
            if (func_col == "COG_label") {
              data <- data %>%
                mutate(COG_category_for_sort = substr(COG_label, 1, 1))
            } else {
              data <- data %>%
                mutate(COG_category_for_sort = COG_category)
            }

            # Order by functional group, then by mean pN/pS within group
            func_order <- data %>%
              group_by(.data[[func_col]], COG_functional_group, COG_category_for_sort) %>%
              summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              arrange(COG_functional_group, desc(avg_pnps)) %>%
              pull(.data[[func_col]])
          } else {
            # For functional groups, order by mean pN/pS
            func_order <- data %>%
              group_by(.data[[func_col]]) %>%
              summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              arrange(desc(avg_pnps)) %>%
              pull(.data[[func_col]])
          }

          # Order taxa by overall mean pN/pS (descending), "Other" last
          taxa_order <- data %>%
            group_by(taxa_display) %>%
            summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(avg_pnps)) %>%
            pull(taxa_display)

          if ("Other" %in% taxa_order) {
            taxa_order <- c(setdiff(taxa_order, "Other"), "Other")
          }

          # Calculate mean pN/pS for each CATEGORY across all taxa (separate column on right)
          category_means <- data %>%
            group_by(.data[[func_col]], .data[[metadata_var]]) %>%
            summarise(mean_across_taxa = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
            mutate(taxa_display = "Mean")

          # Add category means to data for visualization
          data_with_mean <- bind_rows(
            data,
            category_means %>% rename(pNpS_transformed = mean_across_taxa)
          )

          # Set factor levels for ordering
          data_with_mean[[func_col]] <- factor(data_with_mean[[func_col]], levels = rev(func_order))

          # Update taxa_order to include "Mean" on the RIGHT (last position)
          taxa_order_with_mean <- c(taxa_order, "Mean")
          data_with_mean$taxa_display <- factor(data_with_mean$taxa_display, levels = taxa_order_with_mean)

          p <- ggplot(data_with_mean, aes_string(x = "taxa_display", y = func_col, fill = "pNpS_transformed")) +
            geom_tile(color = "white", size = 0.5) +
            facet_wrap(as.formula(paste("~", metadata_var)),
                      scales = "free_x",
                      ncol = min(n_groups, 2)) +
            scale_fill_gradient2(
              low = "blue", mid = "white", high = "red",
              midpoint = if (input$transformation == "log") 0 else 1,
              name = y_label,
              limits = c(min(data_with_mean$pNpS_transformed, na.rm = TRUE),
                        max(data_with_mean$pNpS_transformed, na.rm = TRUE))
            ) +
            labs(
              title = title,
              x = tools::toTitleCase(tax_col),
              y = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway"
            ) +
            theme_minimal(base_size = 13) +
            theme(
              # Ensure consistent text sizes across all facets
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"),
              legend.position = "right",
              legend.title = element_text(size = 11, face = "bold"),
              legend.text = element_text(size = 10),
              strip.text = element_text(size = 12, face = "bold", color = "white"),
              strip.background = element_rect(fill = "#2C3E50", color = NA),
              panel.spacing = unit(1.5, "lines"),
              plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
              plot.margin = margin(10, 10, 20, 10),
              # Prevent text scaling based on facet size
              strip.clip = "off"
            ) +
            # Force consistent axis text size across facets
            coord_cartesian(clip = "off")

        } else {
          # Without metadata: single heatmap

          # Order functional categories
          if (input$cog_display == "individual" && "COG_category" %in% names(data)) {
            if (func_col == "COG_label") {
              data <- data %>%
                mutate(COG_category_for_sort = substr(COG_label, 1, 1))
            } else {
              data <- data %>%
                mutate(COG_category_for_sort = COG_category)
            }

            func_order <- data %>%
              group_by(.data[[func_col]], COG_functional_group, COG_category_for_sort) %>%
              summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              arrange(COG_functional_group, desc(avg_pnps)) %>%
              pull(.data[[func_col]])
          } else {
            func_order <- data %>%
              group_by(.data[[func_col]]) %>%
              summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
              arrange(desc(avg_pnps)) %>%
              pull(.data[[func_col]])
          }

          # Order taxa by mean pN/pS (descending), "Other" last
          taxa_order <- data %>%
            group_by(taxa_display) %>%
            summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(avg_pnps)) %>%
            pull(taxa_display)

          if ("Other" %in% taxa_order) {
            taxa_order <- c(setdiff(taxa_order, "Other"), "Other")
          }

          # Calculate mean pN/pS for each CATEGORY across all taxa (separate column on right)
          category_means <- data %>%
            group_by(.data[[func_col]]) %>%
            summarise(mean_across_taxa = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
            mutate(taxa_display = "Mean")

          # Add category means to data for visualization
          data_with_mean <- bind_rows(
            data,
            category_means %>% rename(pNpS_transformed = mean_across_taxa)
          )

          # Set factor levels for ordering
          data_with_mean[[func_col]] <- factor(data_with_mean[[func_col]], levels = rev(func_order))

          # Update taxa_order to include "Mean" on the RIGHT (last position)
          taxa_order_with_mean <- c(taxa_order, "Mean")
          data_with_mean$taxa_display <- factor(data_with_mean$taxa_display, levels = taxa_order_with_mean)

          p <- ggplot(data_with_mean, aes_string(x = "taxa_display", y = func_col, fill = "pNpS_transformed")) +
            geom_tile(color = "white", size = 0.5) +
            scale_fill_gradient2(
              low = "blue", mid = "white", high = "red",
              midpoint = if (input$transformation == "log") 0 else 1,
              name = y_label,
              limits = c(min(data_with_mean$pNpS_transformed, na.rm = TRUE),
                        max(data_with_mean$pNpS_transformed, na.rm = TRUE))
            ) +
            labs(
              title = title,
              x = tools::toTitleCase(tax_col),
              y = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway"
            ) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
              axis.text.y = element_text(size = 11, hjust = 1),
              axis.title = element_text(size = 13, face = "bold"),
              legend.position = "right",
              legend.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 10),
              legend.key.height = unit(1.2, "cm"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.margin = margin(15, 15, 15, 15),
              panel.grid = element_blank()
            )
        }

        return(p)

      } else if (input$viz_type == "bar") {
        # ===== BAR CHART VISUALIZATION =====

        # Order taxa: "Other" last, rest by mean pN/pS
        taxa_order <- data %>%
          group_by(taxa_display) %>%
          summarise(avg_pnps = mean(pNpS_transformed, na.rm = TRUE), .groups = "drop") %>%
          arrange(desc(avg_pnps)) %>%
          pull(taxa_display)

        if ("Other" %in% taxa_order) {
          taxa_order <- c(setdiff(taxa_order, "Other"), "Other")
        }

        data$taxa_display <- factor(data$taxa_display, levels = taxa_order)

        # Create grouped bar chart with optional metadata faceting
        if (metadata_var != "none" && metadata_var %in% names(data)) {
          n_groups <- length(unique(data[[metadata_var]]))

          p <- ggplot(data, aes_string(x = func_col, y = "pNpS_transformed", fill = "taxa_display")) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.8) +
            facet_wrap(as.formula(paste("~", metadata_var)),
                      scales = "free_x",
                      ncol = min(n_groups, 2)) +
            labs(
              title = title,
              x = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway",
              y = y_label,
              fill = tools::toTitleCase(tax_col)
            ) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 15)),
              axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 15)),
              legend.position = "right",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 11, face = "bold"),
              legend.key.size = unit(0.8, "cm"),
              strip.text = element_text(size = 12, face = "bold", color = "white"),
              strip.background = element_rect(fill = "#2C3E50", color = NA),
              panel.spacing = unit(1.5, "lines"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.margin = margin(10, 10, 20, 10)
            )
        } else {
          p <- ggplot(data, aes_string(x = func_col, y = "pNpS_transformed", fill = "taxa_display")) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.8) +
            labs(
              title = title,
              x = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway",
              y = y_label,
              fill = tools::toTitleCase(tax_col)
            ) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
              axis.text.y = element_text(size = 11),
              axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 15)),
              axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 15)),
              legend.position = "right",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12, face = "bold"),
              legend.key.size = unit(0.8, "cm"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.margin = margin(10, 10, 20, 10)
            )
        }

        # Use distinct colors for taxa
        n_taxa_plot <- length(unique(data$taxa_display))
        if (n_taxa_plot <= 10) {
          p <- p + scale_fill_brewer(palette = "Paired")
        } else if (n_taxa_plot <= 20) {
          p <- p + scale_fill_viridis_d(option = "turbo", end = 0.95)
        } else {
          p <- p + scale_fill_manual(values = rainbow(n_taxa_plot))
        }

        return(p)
      }
    })

    # Significance Heatmap ggplot - matches screen plot for download
    significance_heatmap_ggplot <- reactive({
      req(rv$raw_data)

      tryCatch({
        metadata_var <- input$metadata_var
        if (is.null(metadata_var) || metadata_var == "none") {
          return(ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "No metadata variable selected") +
            theme_void())
        }

        data <- rv$raw_data
        tax_col <- input$tax_level

        # Get functional column
        if (input$analysis_type == "cog") {
          if (input$cog_display == "grouped" && "COG_functional_group" %in% names(data)) {
            func_col <- "COG_functional_group"
          } else {
            func_col <- "COG_category"
            data <- data %>% filter(nchar(COG_category) == 1)
          }
        } else {
          func_col <- "KEGG_Pathway"
        }

        # Filter data
        data <- data %>%
          filter(
            !is.na(.data[[func_col]]), .data[[func_col]] != "-", .data[[func_col]] != "",
            !is.na(.data[[tax_col]]), !is.na(.data[[metadata_var]]),
            !is.na(pNpS), is.finite(pNpS)
          )

        # Determine test type
        meta_values <- unique(data[[metadata_var]])
        is_numeric <- is.numeric(data[[metadata_var]])
        n_groups <- length(meta_values)

        categories <- sort(unique(data[[func_col]]))
        taxa <- sort(unique(data[[tax_col]]))

        # Create p-value matrix
        pvalue_matrix <- matrix(NA, nrow = length(categories), ncol = length(taxa),
                               dimnames = list(categories, taxa))

        for (cat_name in categories) {
          for (tax_name in taxa) {
            test_data <- data %>%
              filter(.data[[func_col]] == cat_name, .data[[tax_col]] == tax_name)

            if (nrow(test_data) >= 3) {
              pval <- tryCatch({
                if (is_numeric) {
                  cor.test(test_data$pNpS, test_data[[metadata_var]], method = "spearman", exact = FALSE)$p.value
                } else if (n_groups == 2) {
                  groups <- split(test_data$pNpS, test_data[[metadata_var]])
                  if (length(groups) == 2 && all(sapply(groups, length) >= 2)) {
                    wilcox.test(groups[[1]], groups[[2]])$p.value
                  } else NA
                } else {
                  kruskal.test(pNpS ~ .data[[metadata_var]], data = test_data)$p.value
                }
              }, error = function(e) NA)
              pvalue_matrix[cat_name, tax_name] <- pval
            }
          }
        }

        # Convert to long format
        pvalue_df <- as.data.frame(pvalue_matrix) %>%
          tibble::rownames_to_column("Category") %>%
          tidyr::pivot_longer(-Category, names_to = "Taxon", values_to = "p_value") %>%
          mutate(
            significance = case_when(
              is.na(p_value) ~ "", p_value < 0.001 ~ "***",
              p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ ""
            ),
            log_p = -log10(p_value + 1e-300)
          )

        # Create test type label
        test_type <- if (is_numeric) "Spearman Correlation" else if (n_groups == 2) "Wilcoxon Test" else "Kruskal-Wallis Test"
        title <- paste0("Statistical Significance Matrix (", test_type, " by ", metadata_var, ")")

        p <- ggplot(pvalue_df, aes(x = Taxon, y = Category, fill = log_p)) +
          geom_tile(color = "white", size = 0.5) +
          geom_text(aes(label = significance), size = 5, fontface = "bold", color = "black") +
          scale_fill_gradient2(
            low = "blue", mid = "white", high = "red",
            midpoint = -log10(0.05), name = "-log10(p)",
            na.value = "grey60"
          ) +
          labs(title = title, x = tools::toTitleCase(tax_col),
               y = if (input$analysis_type == "cog") "COG Category" else "KEGG Pathway") +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(size = 11, face = "bold"),
            legend.position = "right",
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid = element_blank()
          )

        p

      }, error = function(e) {
        cat("Error in significance_heatmap_ggplot:", e$message, "\n")
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
          theme_void()
      })
    })

    setup_download_modal_handler(input, output, session, "gene_plot", gene_plot_ggplot, "gene_pnps_heatmap")
    setup_download_modal_handler(input, output, session, "significance_heatmap", significance_heatmap_ggplot, "gene_pnps_significance")

  })
}
