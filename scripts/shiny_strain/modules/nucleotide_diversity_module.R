# =============================================================================
# Nucleotide Diversity Analysis Module - Completely Revised Version
# =============================================================================

# Load required libraries for statistical annotations
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  # Install ggpubr if not available
  install.packages("ggpubr", repos = "https://cran.rstudio.com/")
}
library(ggpubr)


# Load statistical helper functions
tryCatch({
  if (exists("app_dir")) {
    source(file.path(app_dir, "modules/nucleotide_diversity_stats_helpers.R"), local = TRUE)
  }
}, error = function(e) {
  warning("Could not load statistical helper functions: ", e$message)
})

# Metadata detection function - For Nucleotide Diversity (for diversity_data, original metadata only)
get_diversity_metadata_info <- function(data) {
  # Only use actual original metadata column names
  # Actual columns from CRC_Control113_PRJNA447983_metadata.txt file
  original_metadata_cols <- c(
    "bioproject_accession", "accession_used_in_analysis", "country", "continent",
    "host_age", "host_body_mass_index", "host_sex", "disease_group",
    "AJCC_stage", "age_group"
  )

  # Find only original metadata columns that actually exist in data
  available_cols <- colnames(data)
  meta_cols <- intersect(original_metadata_cols, available_cols)

  # Classify metadata types - AUTO-DETECT numerical vs categorical
  categorical_vars <- character(0)
  numerical_vars <- character(0)

  for (var in meta_cols) {
    if (var %in% colnames(data)) {
      var_data <- data[[var]]
      var_data <- var_data[!is.na(var_data)]

      if (length(var_data) > 0) {
        # Auto-detect if variable is numerical by attempting conversion
        is_numeric_var <- suppressWarnings({
          numeric_data <- as.numeric(as.character(var_data))
          # Check if conversion was successful (not all NA)
          !all(is.na(numeric_data))
        })

        if (is_numeric_var) {
          # Successfully converted to numeric - it's a numerical variable
          numerical_vars <- c(numerical_vars, var)
          cat("DEBUG: Variable", var, "detected as NUMERICAL\n")
        } else {
          # Cannot convert to numeric - it's categorical
          categorical_vars <- c(categorical_vars, var)
          cat("DEBUG: Variable", var, "detected as CATEGORICAL\n")
        }
      }
    }
  }

  cat("DEBUG: Final classification - Numerical:", paste(numerical_vars, collapse = ", "), "\n")
  cat("DEBUG: Final classification - Categorical:", paste(categorical_vars, collapse = ", "), "\n")

  return(list(categorical = categorical_vars, numerical = numerical_vars))
}

# UI Module
nucleotide_diversity_UI <- function(id) {
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
        # Tabs placed in the blue header (title removed per request)
        tags$ul(id = ns("header_tabs"), style = "list-style:none; padding-left:0; margin:6px 0 0 0; display:flex; gap:16px;",
          tags$li(class = "nav-item", actionLink(ns("tab_overview"), label = "Overview", class = "nav-link active")),
          tags$li(class = "nav-item", actionLink(ns("tab_comparison"), label = "Group Comparison", class = "nav-link")),
          tags$li(class = "nav-item", actionLink(ns("tab_taxonomy"), label = "Taxonomy Analysis", class = "nav-link")),
          # COMMENTED OUT: Quality Control tab disabled
          # tags$li(class = "nav-item", actionLink(ns("tab_quality"), label = "Quality Control", class = "nav-link")),
          tags$li(class = "nav-item", actionLink(ns("tab_data"), label = "Data Table", class = "nav-link"))
        )
      ),
      card_body(
    layout_sidebar(
      sidebar = sidebar(
            title = "Analysis Settings",
        
        # Load message
        div(
          id = ns("analysis_message"),
          style = "color: red; font-weight: bold; margin-bottom: 10px;",
          "Please load integrated data first."
        ),
        
        accordion(
          accordion_panel(
            title = "Analysis Parameters",
            icon = bsicons::bs_icon("gear"),
            
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
            
            # Dynamic metadata variable selection
            selectInput(ns("metadata_var"),
                       "Metadata Variable:",
                       choices = NULL,
                       selected = NULL),
            
            # Metadata type information
            div(
              style = "color: blue; font-size: 0.9em; margin: 10px 0; padding: 8px; background: #f5f5f5; border-radius: 4px;",
              htmlOutput(ns("metadata_type_info"))
            ),
            
            # Numerical variable options
            conditionalPanel(
              condition = "output.is_numerical_var == true",
              ns = ns,
              selectInput(ns("numerical_analysis_type"),
                         "Analysis Type for Numerical Variable:",
                             choices = c("Original" = "original"),
                             selected = "original")
            ),
            
              
            # Quality filter
            sliderInput(ns("min_coverage"),
                       "Minimum Coverage:",
                       min = 1, max = 50, value = 5, step = 1),
            
            sliderInput(ns("min_breadth"),
                       "Minimum Breadth:",
                       min = 0.1, max = 1.0, value = 0.5, step = 0.1),
            
            # Dynamic sample size slider (will be rendered based on data)
            uiOutput(ns("min_samples_slider_ui")),

            # Show significant taxa only checkbox
            checkboxInput(ns("show_significant_only"),
                         "Show only significant taxa",
                         value = FALSE),

  
            # Statistical options - Auto-select based on variable type
            div(
              style = "color: blue; font-size: 0.9em; margin: 10px 0; padding: 8px; background: #f5f5f5; border-radius: 4px;",
              htmlOutput(ns("statistical_test_info"))
            ),
            
            # Manual test selection (for advanced users)
            conditionalPanel(
              condition = "output.show_manual_test_selection == true",
              ns = ns,
              selectInput(ns("manual_test_type"),
                         "Manual Test Selection (Advanced):",
                         choices = c("Auto-select based on data" = "auto",
                                   "Welch t-test" = "welch",
                                   "Student t-test" = "student_t",
                                   "Mann-Whitney U" = "mann_whitney",
                                   "Kruskal-Wallis" = "kruskal",
                                   "ANOVA" = "anova"),
                         selected = "auto")
            ),
            
            numericInput(ns("alpha_level"),
                        "Significance Level:",
                        value = 0.05, min = 0.001, max = 0.1, step = 0.01),
            
            
            actionButton(ns("update_analysis"), "Update Analysis", 
                        class = "btn-primary btn-block")
          ),
          
              # Plot Size Controls
              accordion_panel(
                "Plot Size",
                icon = bsicons::bs_icon("arrows-angle-expand"),
                numericInput(
                  ns("plot_width"),
                  label = span("Plot Width (px)", bs_icon("info-circle")),
                  value = 800,
                  min = 100,
                  max = 2000,
                  step = 10
                ) %>% tooltip(HTML("Adjust the width of the plot in pixels.")),
                
                numericInput(
                  ns("plot_height"),
                  label = span("Plot Height (px)", bs_icon("info-circle")),
                  value = 500,
                  min = 100,
                  max = 2000,
                  step = 10
                ) %>% tooltip(HTML("Adjust the height of the plot in pixels.")),
                
                sliderInput(
                  ns("plot_margin"),
                  label = span("Plot Margin (cm)", bs_icon("info-circle")),
                  min = 0,
                  max = 2,
                  value = 0.5,
                  step = 0.1
                ) %>% tooltip(HTML("Adjust the margin around the plot in centimeters."))
              ),
              
              # Plot Theme and Color
              accordion_panel(
                "Plot Theme and Color",
                icon = bsicons::bs_icon("palette"),
                selectInput(
                  inputId = ns("plot_theme"),
                  label = span("Plot Theme", bs_icon("info-circle")),
                  choices = names(theme_options),
                  selected = "Black & White"
                ) %>% tooltip(HTML("Select the overall theme for the plot.")),
                
                checkboxInput(
                  inputId = ns("show_grid"),
                  label = span("Show Grid", bs_icon("info-circle")),
                  value = TRUE
                ) %>% tooltip(HTML("Toggle the visibility of grid lines in the plot.")),
                
                selectInput(
                  inputId = ns("font_family"),
                  label = span("Font Family", bs_icon("info-circle")),
                  choices = c("Arial", "Times New Roman", "Courier", "Helvetica", "Georgia"),
                  selected = "Arial"
                ) %>% tooltip(HTML("Select the font family for all text elements in the plot.")),
                
                selectInput(
                  inputId = ns("color_palette"),
                  label = span("Categorical Color Palette", bs_icon("info-circle")),
                  choices = names(available_palettes$categorical),
                  selected = "pal_npg"
                ) %>% tooltip(HTML("Color palette for categorical variables.")),
                
                selectInput(
                  inputId = ns("gradient_palette"),
                  label = span("Numerical Color Palette", bs_icon("info-circle")),
                  choices = names(available_palettes$gradient)
                ) %>% tooltip(HTML("Select a color gradient for numerical variables."))
              ),
              
              # Point Style
              accordion_panel(
                "Point Style",
                icon = bsicons::bs_icon("circle-fill"),
                sliderInput(ns("point_size"), 
                  label = span("Point Size", bs_icon("info-circle")),
                  min = 0.25, max = 10, value = 1, step = 0.25
                ) %>% tooltip(HTML("Adjust the size of points in the plot.")),
                
                sliderInput(ns("point_alpha"), 
                  label = span("Point Alpha", bs_icon("info-circle")),
                  min = 0, max = 1, value = 0.3, step = 0.05
                ) %>% tooltip(HTML("Adjust the transparency of points in the plot."))
              ),
              
          create_text_controls(ns),
          create_legend_controls(ns)
        )
      ),
      
      # Main content
      div(
        # Message before data is loaded
        conditionalPanel(
          condition = "output.data_loaded == false",
          ns = ns,
          div(
            style = "color: red; font-weight: bold; font-size: 20px; text-align: center; margin: 50px;",
            "Load integrated data first"
          )
        ),
        
        # Content after data is loaded
        conditionalPanel(
          condition = "output.data_loaded == true",
          ns = ns,

              navset_hidden(
                id = ns("main_tabs"),
                selected = "Overview",
                nav_panel_hidden(value = "Overview",
                  # Metadata info value boxes (moved inside Overview)
              layout_column_wrap(
                fill = FALSE,
                value_box(
                  title = "Available Categorical Variables",
                  value = textOutput(ns("num_categorical")),
                  showcase = bs_icon("bar-chart"),
                  theme = "mint"
                ),
                value_box(
                  title = "Available Numerical Variables", 
                  value = textOutput(ns("num_numerical")),
                      showcase = bs_icon("123"),
                  theme = "mint"
                )
              ),
              br(),
              
                  # Summary cards with better layout and order (4 cards like Available Variables)
                  layout_column_wrap(
                    fill = FALSE,
                    value_box(
                      title = tags$span("Total Genomes",
                        tags$small(style="display:block; font-weight:normal; font-size:11px; opacity:0.8;",
                                   "All genome records in dataset")),
                      value = textOutput(ns("overview_total_genomes_text")),
                      showcase = bs_icon("123"),
                      theme = "mint"
                    ),
                    value_box(
                      title = "Total Samples",
                      value = textOutput(ns("overview_total_samples_text")),
                      showcase = bs_icon("stack"),
                      theme = "mint"
                    ),
                    value_box(
                      title = tags$span("Filtered Genomes",
                        tags$small(style="display:block; font-weight:normal; font-size:11px; opacity:0.8;",
                                   "After coverage/breadth filters")),
                      value = textOutput(ns("overview_filtered_genomes_text")),
                      showcase = bs_icon("diagram-2"),
                      theme = "mint"
                    ),
                    value_box(
                      title = "Unique Taxa",
                      value = textOutput(ns("overview_unique_taxa_text")),
                      showcase = bs_icon("arrow-left-right"),
                      theme = "mint"
                    )
                  ),
                  br(),
                  # Plots with better layout (no titles)
                  fluidRow(
                    column(
                      width = 6,
                      div(
                        style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                        create_plot_download_btn(ns, "overview_diversity")
                      ),
                      plotlyOutput(ns("overview_diversity_distribution_plot"), height = "400px")
                    ),
                    column(
                      width = 6,
                      div(
                        style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                        create_plot_download_btn(ns, "overview_quality")
                      ),
                      plotlyOutput(ns("overview_quality_plot"), height = "400px")
                    )
                  )
                ),
                
                nav_panel_hidden(value = "Group Comparison",
                  # Main comparison plot
                  div(
                    style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                    create_plot_download_btn(ns, "metadata_comparison")
                  ),
                  plotlyOutput(ns("metadata_comparison_plot"), height = "600px"),
                  br(),

                  # Statistical results in a single panel
                  div(
                    style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
                    h4("Statistical Test Results", style = "color: #495057; margin-bottom: 15px;"),
                    div(
                      style = "display: flex; justify-content: space-between; align-items: flex-start; gap: 20px;",
                      div(
                        style = "flex: 1;",
                        DT::dataTableOutput(ns("stats_results_table"))
                      ),
                      div(
                        style = "flex: 0 0 300px;",
                        h5("Significance Legend", style = "color: #6c757d; margin-bottom: 10px;"),
                        div(
                          style = "background-color: white; padding: 15px; border-radius: 6px; border: 1px solid #dee2e6;",
                          htmlOutput(ns("significance_legend"))
                        )
                      )
                    )
                  )
                ),

                nav_panel_hidden(value = "Taxonomy Analysis",
                  fluidRow(
                    column(
                      width = 12,
                      div(
                        style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
                        h4("Diversity Metrics Across Taxa", style = "margin: 0;"),
                        create_plot_download_btn(ns, "taxonomy_boxplot")
                      ),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 4px;",
                        HTML(paste0(
                          "<strong>Taxonomy-based Analysis:</strong> This visualization shows diversity metrics across all taxa at the selected taxonomic rank (",
                          "<span id='", ns("current_tax_level_display"), "'></span>",
                          "). Each boxplot represents the distribution of the selected diversity metric for all genomes belonging to that taxon, grouped by metadata categories.",
                          "<br><strong>Statistical Testing:</strong> Mann-Whitney U test for 2-group comparisons, Kruskal-Wallis test for 3+ groups. ",
                          "Significant Kruskal-Wallis results trigger post-hoc Dunn's test to identify specific group differences."
                        ))
                      ),
                      # Individual taxa range slider for SNVs per kbp plot
                      uiOutput(ns("snv_taxa_range_ui")),

                      plotlyOutput(ns("taxonomy_boxplot"), height = "600px")
                    )
                  ),
                  br(),
                  # SNVs per kbp - Statistical Test Results Table
                  fluidRow(
                    column(
                      width = 12,
                      h4("SNVs per kbp - Statistical Test Results"),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 4px;",
                        "Statistical test results for SNVs per kbp comparisons across taxa. Shows p-values, effect sizes, and group comparisons."
                      ),
                      div(style = "height: 300px; overflow-y: auto;",
                        DT::dataTableOutput(ns("snvs_stats_table"))
                      )
                    )
                  ),
                  br(),
                  # Add separate nucleotide diversity plot
                  fluidRow(
                    column(
                      width = 12,
                      div(
                        style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
                        h4("Median Nucleotide Diversity by Taxon", style = "margin: 0;"),
                        create_plot_download_btn(ns, "median_nucl")
                      ),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 4px;",
                        "This plot shows the median nucleotide diversity for each taxon, grouped by the selected metadata variable."
                      ),
                      # Individual taxa range slider for Median Nucleotide Diversity plot
                      uiOutput(ns("nucl_taxa_range_ui")),

                      plotlyOutput(ns("median_nucl_plot"), height = "500px")
                    )
                  ),
                  br(),
                  # Nucleotide Diversity - Statistical Test Results Table
                  fluidRow(
                    column(
                      width = 12,
                      h4("Nucleotide Diversity - Statistical Test Results"),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 4px;",
                        "Statistical test results for nucleotide diversity comparisons across taxa. Shows p-values, effect sizes, and group comparisons."
                      ),
                      div(style = "height: 300px; overflow-y: auto;",
                        DT::dataTableOutput(ns("nucl_stats_table"))
                      )
                    )
                  ),
                  br(),
                  # SNVs vs Nucleotide Diversity Correlation Analysis
                  fluidRow(
                    column(
                      width = 12,
                      div(
                        style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
                        h4("SNVs per kbp vs Nucleotide Diversity Correlation", style = "margin: 0;"),
                        create_plot_download_btn(ns, "snv_nucl_correlation")
                      ),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 4px;",
                        HTML(paste0(
                          "<strong>Correlation Analysis:</strong> This plot shows the relationship between SNVs per kbp and Nucleotide Diversity for each taxon. ",
                          "Each taxon is represented by a different color with a fitted regression line. ",
                          "Spearman correlation coefficients and p-values are calculated for each taxon independently."
                        ))
                      ),
                      plotlyOutput(ns("snv_nucl_correlation_plot"), height = "600px")
                    )
                  ),
                  br(),
                  # Correlation Statistics Table
                  fluidRow(
                    column(
                      width = 12,
                      h4("Correlation Statistics by Taxon"),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 4px;",
                        "Spearman correlation coefficients between SNVs per kbp and Nucleotide Diversity for each taxon. Significant correlations (p < 0.05) are highlighted."
                      ),
                      div(style = "height: 400px; overflow-y: auto;",
                        DT::dataTableOutput(ns("snv_nucl_correlation_table"))
                      )
                    )
                  ),
                  br(),
                  fluidRow(
                    column(12,
                      h4("Taxa Summary"),
                      div(
                        style = "color: #666; font-size: 0.9em; margin-bottom: 10px;",
                        "Summary statistics for each taxon including sample sizes, median diversity, and interquartile ranges."
                      ),
                      DT::dataTableOutput(ns("taxonomy_summary_table"))
                    )
                  )
                ),

                # COMMENTED OUT: Quality Control panel disabled
                # nav_panel_hidden(value = "Quality Control",
                #   fluidRow(
                #     column(
                #       width = 6,
                #       div(
                #         style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                #         create_plot_download_btn(ns, "coverage_breadth")
                #       ),
                #       plotlyOutput(ns("coverage_breadth_plot"), height = "400px")
                #     ),
                #     column(
                #       width = 6,
                #       div(
                #         style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                #         create_plot_download_btn(ns, "genome_length")
                #       ),
                #       plotlyOutput(ns("genome_length_plot"), height = "400px")
                #     )
                #   ),
                #   br(),
                #   div(
                #     style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                #     create_plot_download_btn(ns, "quality_metrics")
                #   ),
                #   plotlyOutput(ns("quality_metrics_plot"), height = "400px")
                # ),

                nav_panel_hidden(value = "Data Table",
                  DT::dataTableOutput(ns("diversity_data_table"))
                )
              )
            )
          )
        )
      )
    )
  )
  
}

# Server Module
nucleotide_diversity_Server <- function(id, integrated_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive value for clicked data
    clicked_data_rv <- reactiveValues(data = NULL)
    
    # Function for error handling
    get_current_inputs <- function(...) {
      vars <- c(...)  
      res <- list()
      for (v in vars) {
        res[[v]] <- input[[v]] 
      }
      return(res)
    }
    
    # Data load status
    output$data_loaded <- reactive({
      !is.null(integrated_data())
    })
    outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
    
    # Metadata information
    metadata_info <- reactive({
      req(integrated_data())
      
      tryCatch({
        data <- integrated_data()
        
        # Check if diversity_data exists
        if("diversity_data" %in% names(data) && !is.null(data$diversity_data)) {
          diversity_data <- data$diversity_data
          message("Detecting metadata in diversity_data, rows: ", nrow(diversity_data))
          get_diversity_metadata_info(diversity_data)
        } else {
          message("No diversity_data found, trying main data")
          get_diversity_metadata_info(data)
        }
      }, error = function(e) {
        showNotification(paste("Metadata detection error:", e$message), type = "warning")
        return(list(categorical = character(0), numerical = character(0)))
      })
    })
    
    # Update value boxes
    output$num_categorical <- renderText({
      meta_info <- metadata_info()
      if(is.null(meta_info)) return("0")
      return(as.character(length(meta_info$categorical)))
    })
    
    output$num_numerical <- renderText({
      meta_info <- metadata_info()
      if(is.null(meta_info)) return("0")
      return(as.character(length(meta_info$numerical)))
    })
    
    # Type of currently selected variable
    current_var_type <- reactive({
      req(input$metadata_var, metadata_info())
      meta_info <- metadata_info()

      cat("DEBUG current_var_type: Checking variable", input$metadata_var, "\n")
      cat("DEBUG current_var_type: Categorical vars:", paste(meta_info$categorical, collapse = ", "), "\n")
      cat("DEBUG current_var_type: Numerical vars:", paste(meta_info$numerical, collapse = ", "), "\n")

      if (input$metadata_var %in% meta_info$categorical) {
        cat("DEBUG current_var_type:", input$metadata_var, "is CATEGORICAL\n")
        return("categorical")
      } else if (input$metadata_var %in% meta_info$numerical) {
        cat("DEBUG current_var_type:", input$metadata_var, "is NUMERICAL\n")
        return("numerical")
      }
      cat("DEBUG current_var_type:", input$metadata_var, "is UNKNOWN\n")
      return("unknown")
    })
    
    # Output whether variable is numerical
    output$is_numerical_var <- reactive({
      tryCatch({
        current_var_type() == "numerical"
      }, error = function(e) {
        return(FALSE)
      })
    })
    outputOptions(output, "is_numerical_var", suspendWhenHidden = FALSE)
    
    output$is_categorical_var <- reactive({
      tryCatch({
        current_var_type() == "categorical"
      }, error = function(e) {
        return(FALSE)
      })
    })
    outputOptions(output, "is_categorical_var", suspendWhenHidden = FALSE)
    
    # Display metadata type information
    output$metadata_type_info <- renderText({
      if(is.null(input$metadata_var) || input$metadata_var == "") {
        return(HTML("<span style='color: #6c757d;'>No variable selected</span>"))
      }
      
      tryCatch({
        var_type <- current_var_type()
        type_info <- switch(var_type,
          "categorical" = list(text = "Categorical Variable", color = "#28a745"),
          "numerical" = list(text = "Numerical Variable", color = "#007bff"),
          list(text = "Unknown Type", color = "#6c757d")
        )
        
        HTML(paste0("<span style='color: ", type_info$color, "; font-weight: bold;'>", 
                    type_info$text, "</span>"))
      }, error = function(e) {
        return(HTML("<span style='color: #dc3545;'>Error detecting type</span>"))
      })
    })
    
    # Metadata analysis for automatic statistical test selection
    metadata_analysis <- reactive({
      req(integrated_data(), input$metadata_var)
      
      tryCatch({
        data <- integrated_data()
        
        # Analyze metadata from diversity_data
        if("diversity_data" %in% names(data) && !is.null(data$diversity_data)) {
          diversity_data <- data$diversity_data
        } else {
          diversity_data <- data
        }
        
        if(!input$metadata_var %in% colnames(diversity_data)) {
          return(list(
            var_type = "unknown",
            unique_count = 0,
            groups = character(0),
            recommended_test = "unknown"
          ))
        }
        
        var_data <- diversity_data[[input$metadata_var]]
        var_data <- var_data[!is.na(var_data)]

        unique_count <- length(unique(var_data))
        groups <- unique(var_data)

        # Variable type determination - check numeric first, then categorical
        # Convert to numeric to check if it's truly numerical
        is_numeric_var <- suppressWarnings({
          numeric_data <- as.numeric(as.character(var_data))
          !all(is.na(numeric_data))
        })

        if(is_numeric_var) {
          # Numeric variable
          if(unique_count <= 2) {
            var_type <- "binary"
            recommended_test <- "Mann-Whitney U"
          } else {
            var_type <- "numerical"
            recommended_test <- "Spearman Correlation"
          }
        } else {
          # Categorical variable
          if(unique_count <= 2) {
            var_type <- "binary"
            recommended_test <- "Mann-Whitney U"
          } else if(unique_count <= 10) {
            var_type <- "categorical"
            recommended_test <- "Pairwise Wilcoxon"
          } else {
            var_type <- "categorical"
            recommended_test <- "Kruskal-Wallis"
          }
        }
        
        return(list(
          var_type = var_type,
          unique_count = unique_count,
          groups = groups,
          recommended_test = recommended_test
        ))
        
      }, error = function(e) {
        return(list(
          var_type = "unknown",
          unique_count = 0,
          groups = character(0),
          recommended_test = "unknown"
        ))
      })
    })
    
    # Display statistical test information (auto-selection)
    output$statistical_test_info <- renderText({
      if(is.null(input$metadata_var) || input$metadata_var == "") {
        return(HTML("<span style='color: #6c757d;'>No variable selected</span>"))
      }
      
      tryCatch({
        analysis <- metadata_analysis()
        
        if(analysis$var_type == "unknown") {
          return(HTML("<span style='color: #dc3545;'>Variable not found in data</span>"))
        }
        
        # Auto-selected test information
        test_info <- HTML(paste0(
          "<strong>Auto-Selected Test:</strong> ", analysis$recommended_test, "<br>",
          "<strong>Variable Type:</strong> ", analysis$var_type, " (", analysis$unique_count, " unique values)<br>",
          if(analysis$var_type != "numerical") {
            paste0("<strong>Groups:</strong> ", paste(head(analysis$groups, 10), collapse = ", "),
                  if(length(analysis$groups) > 10) "..." else "", "<br>")
          } else {
            paste0("<strong>Range:</strong> ", min(analysis$groups, na.rm = TRUE), " to ",
                  max(analysis$groups, na.rm = TRUE), "<br>")
          },
          "<strong>Reason:</strong> ",
          if(analysis$var_type == "binary") {
            "Binary variable → Mann-Whitney U test"
          } else if(analysis$var_type == "categorical") {
            "Categorical variable → Kruskal-Wallis or Pairwise Wilcoxon tests (BH correction)"
          } else if(analysis$var_type == "numerical") {
            "Numerical variable → Spearman correlation (tests association with diversity)"
          }
        ))
        
        test_info
      }, error = function(e) {
        return(HTML("<span style='color: #dc3545;'>Error analyzing variable</span>"))
      })
    })
    
    # Whether to show manual test selection
    output$show_manual_test_selection <- reactive({
      # Manual selection option for advanced users (hidden by default)
      FALSE
    })
    outputOptions(output, "show_manual_test_selection", suspendWhenHidden = FALSE)
    
    # UI control - When data is available
    observe({
      req(integrated_data())

      shinyjs::hide("analysis_message")

      meta_info <- metadata_info()
      all_vars <- c(meta_info$categorical, meta_info$numerical)

      if (length(all_vars) > 0) {
        # Update metadata selection
        choices <- setNames(all_vars, tools::toTitleCase(gsub("_", " ", all_vars)))

        current_selection <- input$metadata_var
        if (is.null(current_selection) || !current_selection %in% all_vars) {
          # Select first categorical if available, otherwise first variable
          default_selection <- if(length(meta_info$categorical) > 0) meta_info$categorical[1] else all_vars[1]
        } else {
          default_selection <- current_selection
        }

        updateSelectInput(session, "metadata_var",
                         choices = choices,
                         selected = default_selection)

        # Enable button
        shinyjs::enable("update_analysis")
        shinyjs::enable("metadata_var")
      } else {
        shinyjs::disable("metadata_var")
        shinyjs::disable("update_analysis")
        showNotification("No valid metadata variables found", type = "warning")
      }
    })

    # Helper function to update active tab styling
    update_active_tab <- function(active_id) {
      tabs <- c("tab_overview", "tab_comparison", "tab_taxonomy", "tab_quality", "tab_data")
      for (tab in tabs) {
        shinyjs::removeClass(id = tab, class = "active")
      }
      shinyjs::addClass(id = active_id, class = "active")
    }

    # Header tab click handlers to control navset_hidden
    observeEvent(input$tab_overview, ignoreInit = TRUE, {
      nav_select("main_tabs", selected = "Overview", session = session)
      update_active_tab("tab_overview")
    })
    observeEvent(input$tab_comparison, ignoreInit = TRUE, {
      nav_select("main_tabs", selected = "Group Comparison", session = session)
      update_active_tab("tab_comparison")
    })
    observeEvent(input$tab_taxonomy, ignoreInit = TRUE, {
      nav_select("main_tabs", selected = "Taxonomy Analysis", session = session)
      update_active_tab("tab_taxonomy")
    })
    # COMMENTED OUT: Quality Control tab event handler disabled
    # observeEvent(input$tab_quality, ignoreInit = TRUE, {
    #   nav_select("main_tabs", selected = "Quality Control", session = session)
    #   update_active_tab("tab_quality")
    # })
    observeEvent(input$tab_data, ignoreInit = TRUE, {
      nav_select("main_tabs", selected = "Data Table", session = session)
      update_active_tab("tab_data")
    })

    # UI control - When data is not available
    observe({
      if (is.null(integrated_data())) {
        shinyjs::show("analysis_message")
        shinyjs::disable("metadata_var")
        shinyjs::disable("update_analysis")
      } else {
        # Navigate to Overview tab when data is loaded
        nav_select("main_tabs", selected = "Overview", session = session)
      }
    })
    
    # Reactive to calculate sample counts per taxon for dynamic slider
    taxa_sample_counts <- reactive({
      req(integrated_data(), input$metadata_var, input$tax_level)

      tryCatch({
        data <- integrated_data()

        # Use diversity_data
        if("diversity_data" %in% names(data) && !is.null(data$diversity_data)) {
          diversity_data <- data$diversity_data
        } else {
          return(data.frame(Taxon = character(), n = integer()))
        }

        # Get basic filtered data (before min_samples filter)
        basic_filtered <- diversity_data %>%
          filter(
            coverage >= input$min_coverage,
            breadth >= input$min_breadth,
            !is.na(.data[[input$tax_level]]),
            !is.na(.data[[input$metadata_var]])
          )

        # Calculate sample counts per taxon
        taxa_counts <- basic_filtered %>%
          group_by(.data[[input$tax_level]]) %>%
          summarise(n_samples = n(), .groups = "drop") %>%
          rename(Taxon = .data[[input$tax_level]])

        return(taxa_counts)

      }, error = function(e) {
        return(data.frame(Taxon = character(), n = integer()))
      })
    })

    # Render the minimum sample slider dynamically
    output$min_samples_slider_ui <- renderUI({
      req(taxa_sample_counts())

      if(nrow(taxa_sample_counts()) == 0) {
        return(sliderInput(session$ns("min_samples_per_taxon"),
                          "Minimum Samples per Taxon:",
                          min = 3, max = 20, value = 5, step = 1))
      }

      max_count <- max(taxa_sample_counts()$n_samples, 3, na.rm = TRUE)
      default_value <- min(5, max_count)

      sliderInput(session$ns("min_samples_per_taxon"),
                  "Minimum Samples per Taxon:",
                  min = 1,
                  max = max_count,
                  value = default_value,
                  step = 1)
    })

    # Render individual taxa range sliders for plots dynamically
    output$snv_taxa_range_ui <- renderUI({
      req(taxa_sample_counts())

      if(nrow(taxa_sample_counts()) == 0) {
        return(sliderInput(session$ns("snv_taxa_range"),
                          "Select Taxa Range (SNVs per kbp):",
                          min = 1, max = 1, value = c(1, 1), step = 1))
      }

      max_taxa <- nrow(taxa_sample_counts())

      sliderInput(session$ns("snv_taxa_range"),
                  "Select Taxa Range (SNVs per kbp):",
                  min = 1,
                  max = max_taxa,
                  value = c(1, min(50, max_taxa)),
                  step = 1)
    })

    output$nucl_taxa_range_ui <- renderUI({
      req(taxa_sample_counts())

      if(nrow(taxa_sample_counts()) == 0) {
        return(sliderInput(session$ns("nucl_taxa_range"),
                          "Select Taxa Range (Nucleotide Diversity):",
                          min = 1, max = 1, value = c(1, 1), step = 1))
      }

      max_taxa <- nrow(taxa_sample_counts())

      sliderInput(session$ns("nucl_taxa_range"),
                  "Select Taxa Range (Nucleotide Diversity):",
                  min = 1,
                  max = max_taxa,
                  value = c(1, min(50, max_taxa)),
                  step = 1)
    })

    # Process analysis data
    processed_data <- eventReactive(input$update_analysis, {
      req(integrated_data(), input$metadata_var)
      
      tryCatch({
        message("Starting nucleotide diversity analysis...")
        data <- integrated_data()
        # Removed problematic nrow() call on list structure

        # Use diversity_data
        if("diversity_data" %in% names(data) && !is.null(data$diversity_data)) {
          diversity_data <- data$diversity_data
          message("Using diversity_data, rows: ", nrow(diversity_data))
        } else {
          showNotification("No diversity data found in integrated data", type = "error")
          return(NULL)
        }

        # Additional quality filtering - Include both metrics
        filtered <- diversity_data %>%
          filter(
            coverage >= input$min_coverage,
            breadth >= input$min_breadth,
            !is.na(.data[[input$tax_level]]),
            !is.na(.data[[input$metadata_var]]),
            !is.na(nucl_diversity),
            !is.na(snvs_per_kbp),
            nucl_diversity > 0,
            snvs_per_kbp > 0
          )

        # Filter by minimum sample count per taxon
        taxon_counts <- filtered %>%
          group_by(.data[[input$tax_level]]) %>%
          summarise(n_samples = n(), .groups = "drop") %>%
          filter(n_samples >= input$min_samples_per_taxon)

        filtered <- filtered %>%
          filter(.data[[input$tax_level]] %in% taxon_counts[[input$tax_level]])

        # Prepare both metrics to be available
        filtered$nucl_diversity_transformed <- filtered$nucl_diversity
        filtered$snvs_per_kbp_transformed <- filtered$snvs_per_kbp
        
        showNotification("Nucleotide diversity analysis completed successfully!", type = "message", duration = 3)
        
        return(list(
          filtered_data = filtered,
          diversity_data = diversity_data,
          tax_level = input$tax_level,
          meta_var = input$metadata_var,
          var_type = current_var_type(),
          summary_stats = list(
            total_genomes = nrow(diversity_data),
            filtered_genomes = nrow(filtered),
            unique_taxa = length(unique(filtered[[input$tax_level]])),
            total_samples = length(unique(filtered$sample_id))
          )
        ))
        
      }, error = function(e) {
        showNotification(paste("Analysis error:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # =============================================================================
    # Statistical Test Selection Logic
    # =============================================================================

    # Significance Notation Function
    get_significance_notation <- function(p_value) {
      if(is.na(p_value)) return("")
      
      if(p_value < 0.001) {
        return("***")  # p < 0.001
      } else if(p_value < 0.01) {
        return("**")   # p < 0.01
      } else if(p_value < 0.05) {
        return("*")    # p < 0.05
      } else if(p_value < 0.1) {
        return(".")     # p < 0.1 (trend)
      } else {
        return("ns")   # p >= 0.1 (not significant)
      }
    }
    
    # Significance notation range information
    significance_legend <- function() {
      HTML(paste0(
        "<strong>Significance Notation (BH-corrected):</strong><br>",
        "• <strong>***</strong> p < 0.001<br>",
        "• <strong>**</strong> p < 0.01<br>",
        "• <strong>*</strong> p < 0.05<br>",
        "• <strong>.</strong> p < 0.1 (trend)<br>",
        "• <strong>ns</strong> p ≥ 0.1 (not significant)<br>",
        "<em>Note: All p-values are corrected using Benjamini-Hochberg (BH) method</em>"
      ))
    }
    #
    # 1. BINARY VARIABLES:
    #    - Mann-Whitney U test: Non-parametric test for comparing 2 groups
    #
    # 2. CATEGORICAL/NUMERICAL VARIABLES:
    #    - Pairwise Wilcoxon test: Non-parametric test for all group pairs
    #    - Reason: Robust non-parametric approach that requires no distribution assumptions
    #
    # 3. EFFECT SIZE CALCULATION:
    #    - Cohen's d (2 groups): Standardized mean difference
    #      * 0.2 = small effect, 0.5 = medium effect, 0.8 = large effect
    #      * Mean difference between two groups divided by standard deviation
    #      * Biological meaning: Magnitude of diversity difference between groups
    #    - Eta-squared (>2 groups): Proportion of variance explained (0-1 range)
    #      * What percentage of diversity variation the metadata variable explains
    #      * 0.01 = small, 0.06 = medium, 0.14 = large effect
    #
    # 4. MULTIPLE TESTING CORRECTION:
    #    - Uses BH correction (Benjamini-Hochberg FDR correction)
    #    - BH correction also applied in pairwise comparisons
    # =============================================================================

    # Statistical analysis execution - For Nucleotide Diversity only
    nucl_stats_results <- reactive({
      req(processed_data(), input$metadata_var, input$tax_level)

      tryCatch({
        data <- processed_data()$filtered_data
        tax_col <- input$tax_level
        meta_col <- input$metadata_var

        cat("DEBUG nucl_stats_results: meta_col =", meta_col, "\n")

        # === Detect if variable is numerical or categorical ===
        meta_info <- get_diversity_metadata_info(data)
        is_numerical <- meta_col %in% meta_info$numerical

        cat("DEBUG nucl_stats_results: is_numerical =", is_numerical, "\n")

        unique_taxa <- unique(data[[tax_col]])
        results <- list()

        for(taxon in unique_taxa) {
          taxon_data <- data[data[[tax_col]] == taxon, ]
          if(nrow(taxon_data) < 3) next

          taxon_data <- taxon_data[!is.na(taxon_data[[meta_col]]), ]
          if(nrow(taxon_data) < 3) next

          if(is_numerical) {
            # === Numerical variable: Spearman Correlation ===
            meta_values <- as.numeric(as.character(taxon_data[[meta_col]]))
            nucl_values <- taxon_data$nucl_diversity_transformed
            valid_indices <- !is.na(meta_values) & !is.na(nucl_values)

            if(sum(valid_indices) >= 3) {
              test_result <- tryCatch({
                cor.test(meta_values[valid_indices], nucl_values[valid_indices],
                        method = "spearman", exact = FALSE)
              }, error = function(e) list(p.value = NA, estimate = NA))

              if(!is.null(test_result) && !is.na(test_result$p.value)) {
                results[[paste0(taxon, "_nucl_diversity")]] <- data.frame(
                  taxon = taxon,
                  metric = "Nucleotide Diversity",
                  test = "Spearman Correlation",
                  test_type = "correlation",
                  p_value = test_result$p.value,
                  p_adjusted = NA,  # Will be set later (NO BH for correlation)
                  rho = as.numeric(test_result$estimate),
                  n_samples = sum(valid_indices),
                  variable_min = min(meta_values[valid_indices]),
                  variable_max = max(meta_values[valid_indices]),
                  # Dummy columns for compatibility
                  mean_group1 = NA, mean_group2 = NA,
                  median_group1 = NA, median_group2 = NA,
                  n_group1 = NA, n_group2 = NA,
                  group1_name = NA, group2_name = NA,
                  group_info = NA,
                  stringsAsFactors = FALSE
                )
              }
            }
          } else {
            # === Categorical variable: Group comparisons ===
            groups <- unique(taxon_data[[meta_col]])
            groups <- groups[!is.na(groups)]
            if(length(groups) < 2) next

            if(length(groups) == 2) {
              # Two groups: Mann-Whitney U
              group1_data <- taxon_data[taxon_data[[meta_col]] == groups[1], "nucl_diversity_transformed"]
              group2_data <- taxon_data[taxon_data[[meta_col]] == groups[2], "nucl_diversity_transformed"]

              if(length(group1_data) >= 2 && length(group2_data) >= 2) {
                nucl_result <- tryCatch({
                  wilcox.test(group1_data, group2_data)
                }, error = function(e) list(p.value = NA, statistic = NA))

                if(!is.null(nucl_result) && !is.na(nucl_result$p.value)) {
                  results[[paste0(taxon, "_nucl_diversity")]] <- data.frame(
                    taxon = taxon,
                    metric = "Nucleotide Diversity",
                    test = "Mann-Whitney U",
                    test_type = "two_group",
                    p_value = nucl_result$p.value,
                    p_adjusted = NA,  # Will be set later (NO BH for 2-group)
                    mean_group1 = mean(group1_data, na.rm = TRUE),
                    mean_group2 = mean(group2_data, na.rm = TRUE),
                    median_group1 = median(group1_data, na.rm = TRUE),
                    median_group2 = median(group2_data, na.rm = TRUE),
                    n_group1 = length(group1_data),
                    n_group2 = length(group2_data),
                    group1_name = as.character(groups[1]),
                    group2_name = as.character(groups[2]),
                    group_info = NA,
                    # Dummy columns for compatibility
                    rho = NA, n_samples = nrow(taxon_data),
                    variable_min = NA, variable_max = NA,
                    stringsAsFactors = FALSE
                  )
                }
              }
            } else {
              # Multiple groups: Kruskal-Wallis
              nucl_result <- tryCatch({
                kruskal.test(taxon_data$nucl_diversity_transformed, taxon_data[[meta_col]])
              }, error = function(e) list(p.value = NA, statistic = NA))

              if(!is.null(nucl_result) && !is.na(nucl_result$p.value)) {
                results[[paste0(taxon, "_nucl_diversity")]] <- data.frame(
                  taxon = taxon,
                  metric = "Nucleotide Diversity",
                  test = "Kruskal-Wallis",
                  test_type = "multi_group",
                  p_value = nucl_result$p.value,
                  p_adjusted = NA,  # Will be set later (BH if multiple taxa)
                  group_info = paste(groups, collapse = ", "),
                  n_samples = nrow(taxon_data),
                  # Dummy columns for compatibility
                  mean_group1 = NA, mean_group2 = NA,
                  median_group1 = NA, median_group2 = NA,
                  n_group1 = NA, n_group2 = NA,
                  group1_name = NA, group2_name = NA,
                  rho = NA, variable_min = NA, variable_max = NA,
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }

        if(length(results) > 0) {
          result_df <- do.call(rbind, results)

          # === Apply BH correction ONLY for multi-group tests with multiple taxa ===
          multi_group_rows <- which(result_df$test_type == "multi_group")

          if(length(multi_group_rows) > 1) {
            # Multiple KW tests → Apply BH correction
            result_df$p_adjusted[multi_group_rows] <- p.adjust(result_df$p_value[multi_group_rows], method = "BH")
            cat("DEBUG: Applied BH correction to", length(multi_group_rows), "Kruskal-Wallis tests\n")
          } else if(length(multi_group_rows) == 1) {
            # Single KW test → NO BH correction
            result_df$p_adjusted[multi_group_rows] <- result_df$p_value[multi_group_rows]
          }

          # For correlation and two-group tests: p_adjusted = p_value
          other_rows <- which(result_df$test_type %in% c("correlation", "two_group"))
          if(length(other_rows) > 0) {
            result_df$p_adjusted[other_rows] <- result_df$p_value[other_rows]
          }

          result_df$significant <- result_df$p_adjusted < input$alpha_level
          result_df$significance_notation <- sapply(result_df$p_adjusted, get_significance_notation)

          return(result_df)
        } else {
          return(data.frame())
        }

      }, error = function(e) {
        showNotification(paste("Nucleotide Diversity stats error:", e$message), type = "error")
        cat("ERROR in nucl_stats_results:", e$message, "\n")
        return(data.frame())
      })
    })

    # Statistical analysis execution - For SNVs per kbp only
    snv_stats_results <- reactive({
      req(processed_data(), input$metadata_var, input$tax_level)

      tryCatch({
        data <- processed_data()$filtered_data
        tax_col <- input$tax_level
        meta_col <- input$metadata_var

        cat("DEBUG snv_stats_results: meta_col =", meta_col, "\n")

        # === Detect if variable is numerical or categorical ===
        meta_info <- get_diversity_metadata_info(data)
        is_numerical <- meta_col %in% meta_info$numerical

        cat("DEBUG snv_stats_results: is_numerical =", is_numerical, "\n")

        unique_taxa <- unique(data[[tax_col]])
        results <- list()

        for(taxon in unique_taxa) {
          taxon_data <- data[data[[tax_col]] == taxon, ]
          if(nrow(taxon_data) < 3) next

          taxon_data <- taxon_data[!is.na(taxon_data[[meta_col]]), ]
          if(nrow(taxon_data) < 3) next

          if(is_numerical) {
            # === Numerical variable: Spearman Correlation ===
            meta_values <- as.numeric(as.character(taxon_data[[meta_col]]))
            snv_values <- taxon_data$snvs_per_kbp_transformed
            valid_indices <- !is.na(meta_values) & !is.na(snv_values)

            if(sum(valid_indices) >= 3) {
              test_result <- tryCatch({
                cor.test(meta_values[valid_indices], snv_values[valid_indices],
                        method = "spearman", exact = FALSE)
              }, error = function(e) list(p.value = NA, estimate = NA))

              if(!is.null(test_result) && !is.na(test_result$p.value)) {
                results[[paste0(taxon, "_snvs_per_kbp")]] <- data.frame(
                  taxon = taxon,
                  metric = "SNVs per kbp",
                  test = "Spearman Correlation",
                  test_type = "correlation",
                  p_value = test_result$p.value,
                  p_adjusted = NA,  # Will be set later (NO BH for correlation)
                  rho = as.numeric(test_result$estimate),
                  n_samples = sum(valid_indices),
                  variable_min = min(meta_values[valid_indices]),
                  variable_max = max(meta_values[valid_indices]),
                  # Dummy columns for compatibility
                  mean_group1 = NA, mean_group2 = NA,
                  median_group1 = NA, median_group2 = NA,
                  n_group1 = NA, n_group2 = NA,
                  group1_name = NA, group2_name = NA,
                  group_info = NA,
                  stringsAsFactors = FALSE
                )
              }
            }
          } else {
            # === Categorical variable: Group comparisons ===
            groups <- unique(taxon_data[[meta_col]])
            groups <- groups[!is.na(groups)]
            if(length(groups) < 2) next

            if(length(groups) == 2) {
              # Two groups: Mann-Whitney U
              group1_data <- taxon_data[taxon_data[[meta_col]] == groups[1], "snvs_per_kbp_transformed"]
              group2_data <- taxon_data[taxon_data[[meta_col]] == groups[2], "snvs_per_kbp_transformed"]

              if(length(group1_data) >= 2 && length(group2_data) >= 2) {
                snv_result <- tryCatch({
                  wilcox.test(group1_data, group2_data)
                }, error = function(e) list(p.value = NA, statistic = NA))

                if(!is.null(snv_result) && !is.na(snv_result$p.value)) {
                  results[[paste0(taxon, "_snvs_per_kbp")]] <- data.frame(
                    taxon = taxon,
                    metric = "SNVs per kbp",
                    test = "Mann-Whitney U",
                    test_type = "two_group",
                    p_value = snv_result$p.value,
                    p_adjusted = NA,  # Will be set later (NO BH for 2-group)
                    mean_group1 = mean(group1_data, na.rm = TRUE),
                    mean_group2 = mean(group2_data, na.rm = TRUE),
                    median_group1 = median(group1_data, na.rm = TRUE),
                    median_group2 = median(group2_data, na.rm = TRUE),
                    n_group1 = length(group1_data),
                    n_group2 = length(group2_data),
                    group1_name = as.character(groups[1]),
                    group2_name = as.character(groups[2]),
                    group_info = NA,
                    # Dummy columns for compatibility
                    rho = NA, n_samples = nrow(taxon_data),
                    variable_min = NA, variable_max = NA,
                    stringsAsFactors = FALSE
                  )
                }
              }
            } else {
              # Multiple groups: Kruskal-Wallis
              snv_result <- tryCatch({
                kruskal.test(taxon_data$snvs_per_kbp_transformed, taxon_data[[meta_col]])
              }, error = function(e) list(p.value = NA, statistic = NA))

              if(!is.null(snv_result) && !is.na(snv_result$p.value)) {
                results[[paste0(taxon, "_snvs_per_kbp")]] <- data.frame(
                  taxon = taxon,
                  metric = "SNVs per kbp",
                  test = "Kruskal-Wallis",
                  test_type = "multi_group",
                  p_value = snv_result$p.value,
                  p_adjusted = NA,  # Will be set later (BH if multiple taxa)
                  group_info = paste(groups, collapse = ", "),
                  n_samples = nrow(taxon_data),
                  # Dummy columns for compatibility
                  mean_group1 = NA, mean_group2 = NA,
                  median_group1 = NA, median_group2 = NA,
                  n_group1 = NA, n_group2 = NA,
                  group1_name = NA, group2_name = NA,
                  rho = NA, variable_min = NA, variable_max = NA,
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }

        if(length(results) > 0) {
          result_df <- do.call(rbind, results)

          # === Apply BH correction ONLY for multi-group tests with multiple taxa ===
          multi_group_rows <- which(result_df$test_type == "multi_group")

          if(length(multi_group_rows) > 1) {
            # Multiple KW tests → Apply BH correction
            result_df$p_adjusted[multi_group_rows] <- p.adjust(result_df$p_value[multi_group_rows], method = "BH")
            cat("DEBUG: Applied BH correction to", length(multi_group_rows), "Kruskal-Wallis tests (SNVs)\n")
          } else if(length(multi_group_rows) == 1) {
            # Single KW test → NO BH correction
            result_df$p_adjusted[multi_group_rows] <- result_df$p_value[multi_group_rows]
          }

          # For correlation and two-group tests: p_adjusted = p_value
          other_rows <- which(result_df$test_type %in% c("correlation", "two_group"))
          if(length(other_rows) > 0) {
            result_df$p_adjusted[other_rows] <- result_df$p_value[other_rows]
          }

          result_df$significant <- result_df$p_adjusted < input$alpha_level
          result_df$significance_notation <- sapply(result_df$p_adjusted, get_significance_notation)

          return(result_df)
        } else {
          return(data.frame())
        }

      }, error = function(e) {
        showNotification(paste("SNVs per kbp stats error:", e$message), type = "error")
        cat("ERROR in snv_stats_results:", e$message, "\n")
        return(data.frame())
      })
    })

    # For backward compatibility - combined results (deprecated)
    stats_results <- reactive({
      # Combine both separate results for backward compatibility
      nucl_results <- nucl_stats_results()
      snv_results <- snv_stats_results()

      if(nrow(nucl_results) > 0 && nrow(snv_results) > 0) {
        return(rbind(nucl_results, snv_results))
      } else if(nrow(nucl_results) > 0) {
        return(nucl_results)
      } else if(nrow(snv_results) > 0) {
        return(snv_results)
      } else {
        return(data.frame())
      }
    })

        
    # Overview metrics with text outputs
    output$overview_total_genomes_text <- renderText({
      data <- processed_data()
      if(is.null(data)) return("0")
      return(as.character(data$summary_stats$total_genomes))
    })

    output$overview_total_samples_text <- renderText({
      data <- processed_data()
      if(is.null(data)) return("0")
      return(as.character(data$summary_stats$total_samples))
    })

    output$overview_filtered_genomes_text <- renderText({
      data <- processed_data()
      if(is.null(data)) return("0")
      return(as.character(data$summary_stats$filtered_genomes))
    })

    output$overview_unique_taxa_text <- renderText({
      data <- processed_data()
      if(is.null(data)) return("0")
      return(as.character(data$summary_stats$unique_taxa))
    })

    # Overview plots
    output$overview_diversity_distribution_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      if(is.null(filtered_data) || nrow(filtered_data) == 0) {
        return(plot_ly() %>%
                 add_text(x = 0.5, y = 0.5, text = "No diversity data available",
                          textfont = list(size = 16, color = "#6c757d")) %>%
                 layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      # Create two distribution plots for both metrics
      p1 <- ggplot(filtered_data, aes(x = nucl_diversity_transformed)) +
        geom_histogram(bins = 50, fill = "#2FA4E7", alpha = 0.7, color = "white") +
        labs(title = "Distribution of Nucleotide Diversity",
             x = "Nucleotide Diversity",
             y = "Count") +
        theme_minimal()

      p2 <- ggplot(filtered_data, aes(x = snvs_per_kbp_transformed)) +
        geom_histogram(bins = 50, fill = "#E762F2", alpha = 0.7, color = "white") +
        labs(title = "Distribution of SNVs per kbp",
             x = "SNVs per kbp",
             y = "Count") +
        theme_minimal()

      # Convert ggplot objects to plotly and combine
      combined_plot <- plotly::subplot(ggplotly(p1), ggplotly(p2), nrows = 2) %>%
        layout(title = "Diversity Metrics Distribution")

      combined_plot
    })

    output$overview_quality_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      if(is.null(filtered_data) || nrow(filtered_data) == 0) {
        return(plot_ly() %>%
                 add_text(x = 0.5, y = 0.5, text = "No quality data available",
                          textfont = list(size = 16, color = "#6c757d")) %>%
                 layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      p <- ggplot(filtered_data, aes(x = coverage, y = breadth, color = .data[[data$meta_var]])) +
        geom_point(alpha = 0.6) +
        geom_hline(yintercept = input$min_breadth, linetype = "dashed", color = "red") +
        geom_vline(xintercept = input$min_coverage, linetype = "dashed", color = "red") +
        labs(title = "Coverage vs Breadth",
             x = "Coverage",
             y = "Breadth",
             color = data$meta_var) +
        theme_minimal()

      ggplotly(p)
    })
    
    # Metadata comparison plot - supports both categorical and numerical variables
    output$metadata_comparison_plot <- renderPlotly({
      data <- processed_data()
      req(data)

      tryCatch({
        filtered_data <- data$filtered_data
        meta_var <- data$meta_var
        var_type <- data$var_type

        # Get statistical results
        stats <- stats_results()

        if(var_type == "categorical") {
          # Categorical variable: create boxplots for both metrics without faceting
          # First plot keeps the legend
          p1 <- ggplot(filtered_data, aes_string(x = meta_var, y = "nucl_diversity_transformed",
                                                fill = meta_var)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
            geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
            labs(title = "Nucleotide Diversity by Groups",
                 x = meta_var,
                 y = "Nucleotide Diversity") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "right")

          # Second plot hides the legend (same as first plot)
          p2 <- ggplot(filtered_data, aes_string(x = meta_var, y = "snvs_per_kbp_transformed",
                                                fill = meta_var)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
            geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
            labs(title = "SNVs per kbp by Groups",
                 x = meta_var,
                 y = "SNVs per kbp") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")  # Hide legend for second plot

          # Convert to plotly and hide legend from second subplot
          p1_plotly <- ggplotly(p1)
          p2_plotly <- ggplotly(p2) %>% layout(showlegend = FALSE)

          # Combine plots vertically with single legend
          combined_plot <- plotly::subplot(p1_plotly, p2_plotly, nrows = 2, shareX = TRUE) %>%
            layout(title = paste("Diversity Metrics by", meta_var))

        } else if(var_type == "numerical") {
          # Numerical variable: create scatter plots
          p1 <- ggplot(filtered_data, aes_string(x = meta_var, y = "nucl_diversity_transformed")) +
            geom_point(alpha = 0.6, color = "#2FA4E7") +
            geom_smooth(method = "lm", se = TRUE, color = "red") +
            labs(title = "Nucleotide Diversity vs Numerical Variable",
                 x = meta_var,
                 y = "Nucleotide Diversity") +
            theme_minimal()

          p2 <- ggplot(filtered_data, aes_string(x = meta_var, y = "snvs_per_kbp_transformed")) +
            geom_point(alpha = 0.6, color = "#E762F2") +
            geom_smooth(method = "lm", se = TRUE, color = "red") +
            labs(title = "SNVs per kbp vs Numerical Variable",
                 x = meta_var,
                 y = "SNVs per kbp") +
            theme_minimal()

          # Combine plots vertically
          combined_plot <- plotly::subplot(ggplotly(p1), ggplotly(p2), nrows = 2, shareX = TRUE) %>%
            layout(title = paste("Diversity Metrics vs", meta_var))
        }

        # Statistical results are shown in the table below, not on the plot
        combined_plot

      }, error = function(e) {
        showNotification(paste("Comparison plot error:", e$message), type = "error")
        return(plotly_empty())
      })
    })

    # Statistical Results Table for Group Comparison
    output$stats_results_table <- DT::renderDataTable({
      stats <- stats_results()

      if(is.null(stats) || nrow(stats) == 0) {
        return(DT::datatable(
          data.frame(Message = "No statistical results available. Run analysis first."),
          options = list(dom = 't', paging = FALSE),
          rownames = FALSE
        ))
      }

      # Format the table for display
      display_df <- stats %>%
        dplyr::select(
          Taxon = taxon,
          Metric = metric,
          Test = test,
          `P-value` = p_value,
          `Adj. P-value` = p_adjusted
        ) %>%
        dplyr::mutate(
          `P-value` = signif(`P-value`, 3),
          `Adj. P-value` = ifelse(is.na(`Adj. P-value`), "-", signif(`Adj. P-value`, 3)),
          Significance = dplyr::case_when(
            `P-value` < 0.001 ~ "***",
            `P-value` < 0.01 ~ "**",
            `P-value` < 0.05 ~ "*",
            TRUE ~ "ns"
          )
        )

      # Add group comparison info if available
      if("mean_group1" %in% names(stats) && "mean_group2" %in% names(stats)) {
        group_info <- stats %>%
          dplyr::mutate(
            `Group Means` = ifelse(
              !is.na(mean_group1) & !is.na(mean_group2),
              paste0(signif(mean_group1, 3), " vs ", signif(mean_group2, 3)),
              "-"
            )
          ) %>%
          dplyr::pull(`Group Means`)

        display_df$`Group Means` <- group_info
      }

      # Add correlation info if available
      if("rho" %in% names(stats)) {
        rho_info <- stats %>%
          dplyr::mutate(
            Correlation = ifelse(!is.na(rho), paste0("ρ = ", signif(rho, 3)), "-")
          ) %>%
          dplyr::pull(Correlation)

        display_df$Correlation <- rho_info
      }

      DT::datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          order = list(list(3, 'asc'))  # Sort by p-value
        ),
        rownames = FALSE,
        filter = 'top',
        extensions = 'Buttons'
      ) %>%
        DT::formatStyle(
          'Significance',
          backgroundColor = DT::styleEqual(
            c('***', '**', '*', 'ns'),
            c('#d4edda', '#fff3cd', '#ffeeba', '#f8f9fa')
          )
        )
    })

    # COMMENTED OUT: Quality Control plots disabled
    # # Coverage vs Breadth plot
    # output$coverage_breadth_plot <- renderPlotly({
    #   data <- processed_data()
    #   req(data)
    #
    #   tryCatch({
    #     filtered_data <- data$filtered_data
    #     p <- ggplot(filtered_data, aes(x = coverage, y = breadth, color = .data[[data$meta_var]])) +
    #       geom_point(alpha = 0.6) +
    #       geom_hline(yintercept = input$min_breadth, linetype = "dashed", color = "red") +
    #       geom_vline(xintercept = input$min_coverage, linetype = "dashed", color = "red") +
    #       labs(title = "Coverage vs Breadth",
    #            x = "Coverage",
    #            y = "Breadth",
    #            color = data$meta_var)
    #
    #     ggplotly(p)
    #
    #   }, error = function(e) {
    #     showNotification(paste("Coverage plot error:", e$message), type = "error")
    #     return(plotly_empty())
    #   })
    # })
    #
    # # Genome length distribution
    # output$genome_length_plot <- renderPlotly({
    #   data <- processed_data()
    #   req(data)
    #
    #   tryCatch({
    #     filtered_data <- data$filtered_data
    #     p <- ggplot(filtered_data, aes(x = genome_size/1e6)) +
    #       geom_histogram(bins = 30, fill = "#2FA4E7", alpha = 0.7, color = "white") +
    #       labs(title = "Genome Length Distribution",
    #            x = "Genome Length (Mbp)",
    #            y = "Count")
    #
    #     ggplotly(p)
    #
    #   }, error = function(e) {
    #     showNotification(paste("Genome length plot error:", e$message), type = "error")
    #     return(plotly_empty())
    #   })
    # })
    #
    # # Quality metrics plot
    # output$quality_metrics_plot <- renderPlotly({
    #   data <- processed_data()
    #   req(data)
    #
    #   tryCatch({
    #     filtered_data <- data$filtered_data
    #     quality_summary <- filtered_data %>%
    #       group_by(.data[[data$meta_var]]) %>%
    #       summarise(
    #         mean_coverage = mean(coverage, na.rm = TRUE),
    #         mean_breadth = mean(breadth, na.rm = TRUE),
    #         mean_genome_length = mean(genome_size, na.rm = TRUE),
    #         .groups = "drop"
    #       ) %>%
    #       pivot_longer(cols = starts_with("mean_"), names_to = "metric", values_to = "value") %>%
    #       mutate(metric = gsub("mean_", "", metric))
    #
    #     p <- ggplot(quality_summary, aes(x = .data[[data$meta_var]], y = value, fill = .data[[data$meta_var]])) +
    #       geom_col() +
    #       facet_wrap(~metric, scales = "free_y") +
    #       labs(title = "Quality Metrics by Group",
    #            x = data$meta_var,
    #            y = "Mean Value") +
    #       theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #
    #     ggplotly(p)
    #
    #   }, error = function(e) {
    #     showNotification(paste("Quality metrics plot error:", e$message), type = "error")
    #     return(plotly_empty())
    #   })
    # })
    
    # SNVs per kbp - Statistical Test Results Table
    output$snvs_stats_table <- DT::renderDataTable({
      stats <- stats_results()
      if(nrow(stats) == 0) {
        return(DT::datatable(data.frame(Message = "No statistical results available")))
      }

      tryCatch({
        # Filter for SNVs per kbp results only
        snv_stats <- stats %>% filter(metric == "SNVs per kbp")

        if(nrow(snv_stats) == 0) {
          return(DT::datatable(data.frame(Message = "No SNVs per kbp results available")))
        }

        # Check test_type to determine display columns
        test_types <- unique(snv_stats$test_type)

        if("correlation" %in% test_types) {
          # === Numerical variable: Show correlation columns ===
          display_stats <- snv_stats %>%
            mutate(
              p_value = round(p_value, 6),
              p_adjusted = round(p_adjusted, 6),
              rho = round(rho, 3),
              variable_range = paste0(round(variable_min, 1), " - ", round(variable_max, 1)),
              p_adjusted_with_notation = paste0(round(p_adjusted, 4), " ", significance_notation)
            ) %>%
            select(taxon, test, rho, n_samples, variable_range,
                   p_value, p_adjusted, p_adjusted_with_notation, significant)

        } else if("two_group" %in% test_types) {
          # === 2-group categorical: Show group comparison columns ===
          display_stats <- snv_stats %>%
            mutate(
              p_value = round(p_value, 6),
              p_adjusted = round(p_adjusted, 6),
              mean_group1 = round(mean_group1, 4),
              mean_group2 = round(mean_group2, 4),
              median_group1 = round(median_group1, 4),
              median_group2 = round(median_group2, 4),
              p_adjusted_with_notation = paste0(round(p_adjusted, 4), " ", significance_notation)
            ) %>%
            select(taxon, test, group1_name, n_group1, mean_group1, median_group1,
                   group2_name, n_group2, mean_group2, median_group2,
                   p_value, p_adjusted, p_adjusted_with_notation, significant)

        } else {
          # === Multi-group categorical: Show Kruskal-Wallis columns ===
          display_stats <- snv_stats %>%
            mutate(
              p_value = round(p_value, 6),
              p_adjusted = round(p_adjusted, 6),
              p_adjusted_with_notation = paste0(round(p_adjusted, 4), " ", significance_notation)
            ) %>%
            select(taxon, test, n_samples, group_info,
                   p_value, p_adjusted, p_adjusted_with_notation, significant)
        }

        DT::datatable(
          display_stats,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel'),
            scrollY = "300px"
          ),
          rownames = FALSE,
          extensions = 'Buttons',
          caption = "SNVs per kbp - Statistical Test Results"
        ) %>%
        DT::formatStyle(
          'significant',
          backgroundColor = DT::styleEqual(TRUE, '#d4edda'),
          color = DT::styleEqual(TRUE, '#155724')
        ) %>%
        DT::formatStyle(
          'p_adjusted_with_notation',
          backgroundColor = DT::styleEqual(c("***", "**", "*", "."),
                                        c('#d73027', '#fc8d59', '#fee08b', '#d9ef8b')),
          color = DT::styleEqual(c("***", "**", "*", "."), c('white', 'white', 'black', 'black'))
        )

      }, error = function(e) {
        showNotification(paste("SNVs stats table error:", e$message), type = "error")
        cat("ERROR in snvs_stats_table:", e$message, "\n")
        return(DT::datatable(data.frame(Error = "Failed to generate SNVs table")))
      })
    })

    # Nucleotide Diversity - Statistical Test Results Table
    output$nucl_stats_table <- DT::renderDataTable({
      stats <- stats_results()
      if(nrow(stats) == 0) {
        return(DT::datatable(data.frame(Message = "No statistical results available")))
      }

      tryCatch({
        # Filter for Nucleotide Diversity results only
        nucl_stats <- stats %>% filter(metric == "Nucleotide Diversity")

        if(nrow(nucl_stats) == 0) {
          return(DT::datatable(data.frame(Message = "No Nucleotide Diversity results available")))
        }

        # Check test_type to determine display columns
        test_types <- unique(nucl_stats$test_type)

        if("correlation" %in% test_types) {
          # === Numerical variable: Show correlation columns ===
          display_stats <- nucl_stats %>%
            mutate(
              p_value = round(p_value, 6),
              p_adjusted = round(p_adjusted, 6),
              rho = round(rho, 3),
              variable_range = paste0(round(variable_min, 1), " - ", round(variable_max, 1)),
              p_adjusted_with_notation = paste0(round(p_adjusted, 4), " ", significance_notation)
            ) %>%
            select(taxon, test, rho, n_samples, variable_range,
                   p_value, p_adjusted, p_adjusted_with_notation, significant)

        } else if("two_group" %in% test_types) {
          # === 2-group categorical: Show group comparison columns ===
          display_stats <- nucl_stats %>%
            mutate(
              p_value = round(p_value, 6),
              p_adjusted = round(p_adjusted, 6),
              mean_group1 = round(mean_group1, 4),
              mean_group2 = round(mean_group2, 4),
              median_group1 = round(median_group1, 4),
              median_group2 = round(median_group2, 4),
              p_adjusted_with_notation = paste0(round(p_adjusted, 4), " ", significance_notation)
            ) %>%
            select(taxon, test, group1_name, n_group1, mean_group1, median_group1,
                   group2_name, n_group2, mean_group2, median_group2,
                   p_value, p_adjusted, p_adjusted_with_notation, significant)

        } else {
          # === Multi-group categorical: Show Kruskal-Wallis columns ===
          display_stats <- nucl_stats %>%
            mutate(
              p_value = round(p_value, 6),
              p_adjusted = round(p_adjusted, 6),
              p_adjusted_with_notation = paste0(round(p_adjusted, 4), " ", significance_notation)
            ) %>%
            select(taxon, test, n_samples, group_info,
                   p_value, p_adjusted, p_adjusted_with_notation, significant)
        }

        DT::datatable(
          display_stats,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel'),
            scrollY = "300px"
          ),
          rownames = FALSE,
          extensions = 'Buttons',
          caption = "Nucleotide Diversity - Statistical Test Results"
        ) %>%
        DT::formatStyle(
          'significant',
          backgroundColor = DT::styleEqual(TRUE, '#d4edda'),
          color = DT::styleEqual(TRUE, '#155724')
        ) %>%
        DT::formatStyle(
          'p_adjusted_with_notation',
          backgroundColor = DT::styleEqual(c("***", "**", "*", "."),
                                        c('#d73027', '#fc8d59', '#fee08b', '#d9ef8b')),
          color = DT::styleEqual(c("***", "**", "*", "."), c('white', 'white', 'black', 'black'))
        )

      }, error = function(e) {
        showNotification(paste("Nucleotide Diversity stats table error:", e$message), type = "error")
        cat("ERROR in nucl_stats_table:", e$message, "\n")
        return(DT::datatable(data.frame(Error = "Failed to generate Nucleotide Diversity table")))
      })
    })
    
    # Output significance notation range information
    output$significance_legend <- renderText({
      significance_legend()
    })
    
      
    
    # Data table
    output$diversity_data_table <- DT::renderDataTable({
      data <- processed_data()
      req(data)

      tryCatch({
        filtered_data <- data$filtered_data
        # Select main columns only - include both transformed metrics
        display_cols <- c("sample_id", "clean_genome", data$tax_level,
                         data$meta_var, "nucl_diversity", "snvs_per_kbp",
                         "nucl_diversity_transformed", "snvs_per_kbp_transformed",
                         "SNV_count", "coverage", "breadth", "genome_size")

        available_cols <- intersect(display_cols, colnames(filtered_data))
        display_data <- filtered_data[, available_cols]

        # Round numeric columns
        numeric_cols <- sapply(display_data, is.numeric)
        display_data[numeric_cols] <- lapply(display_data[numeric_cols], function(x) round(x, 4))

        DT::datatable(
          display_data,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
          ),
          rownames = FALSE,
          extensions = 'Buttons'
        )

      }, error = function(e) {
        showNotification(paste("Data table error:", e$message), type = "error")
        return(DT::datatable(data.frame(Error = "Failed to generate data table")))
      })
    }, options = list(pageLength = 15, scrollX = TRUE))
    
    # Download handler for Selected Point Details (CSV)
    output$download_selected_details <- downloadHandler(
      filename = function() {
        paste0("nucleotide_diversity_details_", Sys.Date(), ".csv")
      },
      content = function(file) {
        dd <- clicked_data_rv$data
        if (is.null(dd)) {
          write.csv(data.frame(), file, row.names = FALSE)
        } else {
          utils::write.csv(dd, file, row.names = FALSE, fileEncoding = "UTF-8")
        }
      }
    )
    
    # =============================================================================
    # TAXONOMY ANALYSIS REACTIVES AND OUTPUTS
    # =============================================================================

    # Taxonomy data processing for boxplots with metadata grouping
    taxonomy_data_raw <- eventReactive(input$update_analysis, {
      req(processed_data())

      tryCatch({
        data <- processed_data()$filtered_data
        tax_col <- input$tax_level
        meta_col <- input$metadata_var

        # Create summaries for both metrics
        taxa_summary_nucl <- data %>%
          group_by(.data[[tax_col]]) %>%
          summarise(
            n_samples = n(),
            n_metadata_groups = n_distinct(.data[[meta_col]]),
            median_diversity = median(nucl_diversity_transformed, na.rm = TRUE),
            mean_diversity = mean(nucl_diversity_transformed, na.rm = TRUE),
            q25_diversity = quantile(nucl_diversity_transformed, 0.25, na.rm = TRUE),
            q75_diversity = quantile(nucl_diversity_transformed, 0.75, na.rm = TRUE),
            min_diversity = min(nucl_diversity_transformed, na.rm = TRUE),
            max_diversity = max(nucl_diversity_transformed, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          filter(n_samples >= input$min_samples_per_taxon,
                 n_metadata_groups >= 2) %>%
          arrange(desc(median_diversity))

        taxa_summary_snv <- data %>%
          group_by(.data[[tax_col]]) %>%
          summarise(
            n_samples = n(),
            n_metadata_groups = n_distinct(.data[[meta_col]]),
            median_diversity = median(snvs_per_kbp_transformed, na.rm = TRUE),
            mean_diversity = mean(snvs_per_kbp_transformed, na.rm = TRUE),
            q25_diversity = quantile(snvs_per_kbp_transformed, 0.25, na.rm = TRUE),
            q75_diversity = quantile(snvs_per_kbp_transformed, 0.75, na.rm = TRUE),
            min_diversity = min(snvs_per_kbp_transformed, na.rm = TRUE),
            max_diversity = max(snvs_per_kbp_transformed, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          filter(n_samples >= input$min_samples_per_taxon,
                 n_metadata_groups >= 2) %>%
          arrange(desc(median_diversity))

        # Filter data to only include qualified taxa (intersection of both metrics)
        qualified_taxa <- intersect(taxa_summary_nucl[[tax_col]], taxa_summary_snv[[tax_col]])
        plot_data <- data[data[[tax_col]] %in% qualified_taxa, ]

        
        # Create synchronized taxa ordering based on nucleotide diversity median values
        # Safety check: ensure the sorting column exists
        if ("median_diversity" %in% colnames(taxa_summary_nucl)) {
          taxa_order <- taxa_summary_nucl %>%
            arrange(desc(median_diversity)) %>%
            pull(.data[[tax_col]])
        } else {
          # Fallback: order alphabetically if sorting column not found
          taxa_order <- taxa_summary_nucl %>%
            arrange(.data[[tax_col]]) %>%
            pull(.data[[tax_col]])
          cat("Warning: median_diversity column not found, using alphabetical ordering\n")
        }

        # Get metadata group info
        metadata_groups <- unique(plot_data[[meta_col]])
        metadata_type <- current_var_type()

        # Return processed data for plotting
        list(
          plot_data = plot_data,
          summary_data_nucl = taxa_summary_nucl,
          summary_data_snv = taxa_summary_snv,
          tax_level = tax_col,
          meta_col = meta_col,
          meta_groups = metadata_groups,
          meta_type = metadata_type,
          metric_names = c("Nucleotide Diversity", "SNVs per kbp"),
          taxa_order = taxa_order  # Add synchronized ordering
        )

      }, error = function(e) {
        showNotification(paste("Taxonomy data processing error:", e$message), type = "error")
        return(NULL)
      })
    })

    # Filtered taxonomy data reactive that applies "Show Significant Only" filtering
    taxonomy_data <- reactive({
      req(taxonomy_data_raw())

      # Explicitly depend on the checkbox to trigger reactivity when it changes
      input$show_significant_only

      tryCatch({
        data <- taxonomy_data_raw()

        # Apply "Show Significant Only" filtering without creating circular dependency
        if (input$show_significant_only) {
          sig_taxa <- significant_taxa()
          cat("Debug: Show Significant Only =", input$show_significant_only, "\n")
          cat("Debug: Found", length(sig_taxa), "significant taxa\n")

          if (length(sig_taxa) > 0 && !is.null(data)) {
            # Filter the baseline data to include only significant taxa
            # For general taxonomy filtering: require taxa to be in both summary datasets for consistency
            tax_col <- data$tax_level
            qualified_taxa <- intersect(
              intersect(data$summary_data_nucl[[tax_col]], data$summary_data_snv[[tax_col]]),
              sig_taxa
            )

            cat("Debug: Qualified taxa after filtering:", length(qualified_taxa), "\n")
            cat("Debug: Original plot data rows:", nrow(data$plot_data), "\n")

            # Filter plot data
            plot_data <- data$plot_data[data$plot_data[[tax_col]] %in% qualified_taxa, ]
            cat("Debug: Filtered plot data rows:", nrow(plot_data), "\n")

            # Filter summaries
            taxa_summary_nucl <- data$summary_data_nucl[data$summary_data_nucl[[tax_col]] %in% qualified_taxa, ]
            taxa_summary_snv <- data$summary_data_snv[data$summary_data_snv[[tax_col]] %in% qualified_taxa, ]

            # Update taxa ordering based on filtered data
            taxa_order <- taxa_summary_nucl %>%
              arrange(desc(median_diversity)) %>%
              pull(.data[[tax_col]])

            # Return filtered data
            return(list(
              plot_data = plot_data,
              summary_data_nucl = taxa_summary_nucl,
              summary_data_snv = taxa_summary_snv,
              tax_level = data$tax_level,
              meta_col = data$meta_col,
              meta_groups = data$meta_groups,
              meta_type = data$meta_type,
              metric_names = data$metric_names,
              taxa_order = taxa_order
            ))
          } else {
            # Checkbox is checked but no significant taxa found - return empty dataset
            cat("Debug: Show Significant Only is ON but no significant taxa found. Returning empty dataset.\n")
            tax_col <- data$tax_level

            # Return empty data structure with same format but zero rows
            return(list(
              plot_data = data$plot_data[0, ],
              summary_data_nucl = data$summary_data_nucl[0, ],
              summary_data_snv = data$summary_data_snv[0, ],
              tax_level = data$tax_level,
              meta_col = data$meta_col,
              meta_groups = data$meta_groups,
              meta_type = data$meta_type,
              metric_names = data$metric_names,
              taxa_order = character(0)  # Empty taxa order
            ))
          }
        }

        # Return baseline data when not filtering
        return(data)

      }, error = function(e) {
        showNotification(paste("Taxonomy filtering error:", e$message), type = "error")
        return(taxonomy_data_raw())  # Fallback to baseline data
      })
    })

    # Taxon-specific statistical tests comparing metadata groups - SIMPLIFIED (RAW P-VALUES FOR WILCOXON)
    taxonomy_metadata_stats <- reactive({
      req(taxonomy_data_raw())

      tryCatch({
        tax_data <- taxonomy_data_raw()
        data <- tax_data$plot_data
        tax_col <- tax_data$tax_level
        meta_col <- tax_data$meta_col
        meta_type <- tax_data$meta_type

        cat("DEBUG taxonomy_metadata_stats: meta_col =", meta_col, "\n")
        cat("DEBUG taxonomy_metadata_stats: meta_type =", meta_type, "\n")

        # Get unique taxa
        unique_taxa <- unique(data[[tax_col]])
        results <- list()

        # Check if metadata is numerical or categorical
        is_numerical <- (meta_type == "numerical")
        cat("DEBUG taxonomy_metadata_stats: is_numerical =", is_numerical, "\n")

        for(taxon in unique_taxa) {
          taxon_data <- data[data[[tax_col]] == taxon, ]

          if(nrow(taxon_data) < 3) next

          # Remove NA values from metadata column
          taxon_data <- taxon_data[!is.na(taxon_data[[meta_col]]), ]
          if(nrow(taxon_data) < 3) next

          if(is_numerical) {
            # === For numerical variables: Use Spearman Correlation ===
            # Convert metadata to numeric if it's not already
            meta_values <- as.numeric(as.character(taxon_data[[meta_col]]))

            # Test for Nucleotide Diversity using Spearman correlation
            nucl_values <- taxon_data$nucl_diversity_transformed
            valid_indices <- !is.na(meta_values) & !is.na(nucl_values)

            if(sum(valid_indices) >= 3) {
              test_result <- tryCatch({
                cor.test(meta_values[valid_indices], nucl_values[valid_indices],
                        method = "spearman", exact = FALSE)
              }, error = function(e) {
                return(list(p.value = NA, estimate = NA))
              })

              if(!is.null(test_result) && !is.na(test_result$p.value)) {
                results[[paste0(taxon, "_nucl_diversity")]] <- list(
                  taxon = taxon,
                  metric = "Nucleotide Diversity",
                  test_name = "Spearman Correlation",
                  test_type = "correlation",  # For table display
                  statistic = test_result$estimate,  # rho value
                  p_value = test_result$p.value,
                  p_adjusted = test_result$p.value,  # NO BH correction for correlations
                  significant = test_result$p.value < input$alpha_level,
                  n_samples = sum(valid_indices),
                  variable_min = min(meta_values[valid_indices]),
                  variable_max = max(meta_values[valid_indices]),
                  rho = as.numeric(test_result$estimate)
                )
              }
            }

            # Test for SNVs per kbp using Spearman correlation
            snv_values <- taxon_data$snvs_per_kbp_transformed
            valid_indices <- !is.na(meta_values) & !is.na(snv_values)

            if(sum(valid_indices) >= 3) {
              test_result <- tryCatch({
                cor.test(meta_values[valid_indices], snv_values[valid_indices],
                        method = "spearman", exact = FALSE)
              }, error = function(e) {
                return(list(p.value = NA, estimate = NA))
              })

              if(!is.null(test_result) && !is.na(test_result$p.value)) {
                results[[paste0(taxon, "_snvs_per_kbp")]] <- list(
                  taxon = taxon,
                  metric = "SNVs per kbp",
                  test_name = "Spearman Correlation",
                  test_type = "correlation",  # For table display
                  statistic = test_result$estimate,  # rho value
                  p_value = test_result$p.value,
                  p_adjusted = test_result$p.value,  # NO BH correction for correlations
                  significant = test_result$p.value < input$alpha_level,
                  n_samples = sum(valid_indices),
                  variable_min = min(meta_values[valid_indices]),
                  variable_max = max(meta_values[valid_indices]),
                  rho = as.numeric(test_result$estimate)
                )
              }
            }
          } else {
            # For categorical variables, use group comparison tests
            # Get metadata groups for this taxon
            meta_groups <- unique(taxon_data[[meta_col]])
            meta_groups <- meta_groups[!is.na(meta_groups)]

            if(length(meta_groups) < 2) next  # Need at least 2 groups for comparison

          # Test for Nucleotide Diversity
          if(length(meta_groups) == 2) {
            # === Two groups: Mann-Whitney U test (NO BH correction) ===
            group1_data <- taxon_data[taxon_data[[meta_col]] == meta_groups[1], "nucl_diversity_transformed"]
            group2_data <- taxon_data[taxon_data[[meta_col]] == meta_groups[2], "nucl_diversity_transformed"]

            if(length(group1_data) >= 2 && length(group2_data) >= 2) {
              test_result <- tryCatch({
                wilcox.test(group1_data, group2_data)
              }, error = function(e) {
                return(list(p.value = NA, statistic = NA))
              })

              if(!is.null(test_result) && !is.na(test_result$p.value)) {
                # For Wilcoxon: significance based on raw p-value only (single test)
                results[[paste0(taxon, "_nucl_diversity")]] <- list(
                  taxon = taxon,
                  metric = "Nucleotide Diversity",
                  test_name = "Mann-Whitney U",
                  test_type = "two_group",  # For table display
                  statistic = test_result$statistic,
                  p_value = test_result$p.value,
                  p_adjusted = test_result$p.value,  # NO BH correction for single test
                  significant = test_result$p.value < input$alpha_level,
                  n_groups = 2,
                  n_samples = nrow(taxon_data),
                  group1_name = as.character(meta_groups[1]),
                  group2_name = as.character(meta_groups[2]),
                  n_group1 = length(group1_data),
                  n_group2 = length(group2_data),
                  mean_group1 = mean(group1_data, na.rm = TRUE),
                  mean_group2 = mean(group2_data, na.rm = TRUE),
                  median_group1 = median(group1_data, na.rm = TRUE),
                  median_group2 = median(group2_data, na.rm = TRUE)
                )
              }
            }
          } else {
            # === Multiple groups (3+): Kruskal-Wallis test ===
            # BH correction applied later if multiple taxa tested
            test_result <- tryCatch({
              kruskal.test(taxon_data$nucl_diversity_transformed, taxon_data[[meta_col]])
            }, error = function(e) {
              return(list(p.value = NA, statistic = NA))
            })

            if(!is.null(test_result) && !is.na(test_result$p.value)) {
              results[[paste0(taxon, "_nucl_diversity")]] <- list(
                taxon = taxon,
                metric = "Nucleotide Diversity",
                test_name = "Kruskal-Wallis",
                test_type = "multi_group",  # For table display and BH correction
                statistic = test_result$statistic,
                p_value = test_result$p.value,
                p_adjusted = NA,  # Will be set after BH correction (if multiple taxa)
                significant = FALSE,  # Will be set after determining p_adjusted
                n_groups = length(meta_groups),
                n_samples = nrow(taxon_data),
                group_info = paste(meta_groups, collapse = ", "),
                group_means = sapply(meta_groups, function(g) {
                  mean(taxon_data[taxon_data[[meta_col]] == g, "nucl_diversity_transformed"], na.rm = TRUE)
                }),
                group_medians = sapply(meta_groups, function(g) {
                  median(taxon_data[taxon_data[[meta_col]] == g, "nucl_diversity_transformed"], na.rm = TRUE)
                }),
                group_sizes = sapply(meta_groups, function(g) {
                  sum(taxon_data[[meta_col]] == g, na.rm = TRUE)
                })
              )
            }
          }

          # Test for SNVs per kbp
          if(length(meta_groups) == 2) {
            # === Two groups: Mann-Whitney U test (NO BH correction) ===
            group1_data <- taxon_data[taxon_data[[meta_col]] == meta_groups[1], "snvs_per_kbp_transformed"]
            group2_data <- taxon_data[taxon_data[[meta_col]] == meta_groups[2], "snvs_per_kbp_transformed"]

            if(length(group1_data) >= 2 && length(group2_data) >= 2) {
              test_result <- tryCatch({
                wilcox.test(group1_data, group2_data)
              }, error = function(e) {
                return(list(p.value = NA, statistic = NA))
              })

              if(!is.null(test_result) && !is.na(test_result$p.value)) {
                # For Wilcoxon: significance based on raw p-value only (single test)
                results[[paste0(taxon, "_snvs_per_kbp")]] <- list(
                  taxon = taxon,
                  metric = "SNVs per kbp",
                  test_name = "Mann-Whitney U",
                  test_type = "two_group",  # For table display
                  statistic = test_result$statistic,
                  p_value = test_result$p.value,
                  p_adjusted = test_result$p.value,  # NO BH correction for single test
                  significant = test_result$p.value < input$alpha_level,
                  n_groups = 2,
                  n_samples = nrow(taxon_data),
                  group1_name = as.character(meta_groups[1]),
                  group2_name = as.character(meta_groups[2]),
                  n_group1 = length(group1_data),
                  n_group2 = length(group2_data),
                  mean_group1 = mean(group1_data, na.rm = TRUE),
                  mean_group2 = mean(group2_data, na.rm = TRUE),
                  median_group1 = median(group1_data, na.rm = TRUE),
                  median_group2 = median(group2_data, na.rm = TRUE)
                )
              }
            }
          } else {
            # === Multiple groups (3+): Kruskal-Wallis test ===
            # BH correction applied later if multiple taxa tested
            test_result <- tryCatch({
              kruskal.test(taxon_data$snvs_per_kbp_transformed, taxon_data[[meta_col]])
            }, error = function(e) {
              return(list(p.value = NA, statistic = NA))
            })

            if(!is.null(test_result) && !is.na(test_result$p.value)) {
              results[[paste0(taxon, "_snvs_per_kbp")]] <- list(
                taxon = taxon,
                metric = "SNVs per kbp",
                test_name = "Kruskal-Wallis",
                test_type = "multi_group",  # For table display and BH correction
                statistic = test_result$statistic,
                p_value = test_result$p.value,
                p_adjusted = NA,  # Will be set after BH correction (if multiple taxa)
                significant = FALSE,  # Will be set after determining p_adjusted
                n_groups = length(meta_groups),
                n_samples = nrow(taxon_data),
                group_info = paste(meta_groups, collapse = ", "),
                group_means = sapply(meta_groups, function(g) {
                  mean(taxon_data[taxon_data[[meta_col]] == g, "snvs_per_kbp_transformed"], na.rm = TRUE)
                }),
                group_medians = sapply(meta_groups, function(g) {
                  median(taxon_data[taxon_data[[meta_col]] == g, "snvs_per_kbp_transformed"], na.rm = TRUE)
                }),
                group_sizes = sapply(meta_groups, function(g) {
                  sum(taxon_data[[meta_col]] == g, na.rm = TRUE)
                })
              )
            }
          }
          }  # End of else block (categorical variables)
        }  # End of for loop

        # === Apply BH correction ONLY for multiple Kruskal-Wallis tests ===
        # Correlation tests (numerical): NO BH correction - already set p_adjusted = p_value
        # Two-group tests (Wilcoxon): NO BH correction - already set p_adjusted = p_value
        # Multi-group tests (Kruskal-Wallis): BH correction ONLY if testing multiple taxa

        # Extract multi-group tests
        multi_group_tests <- results[sapply(results, function(x) {
          !is.null(x$test_type) && x$test_type == "multi_group"
        })]

        if(length(multi_group_tests) > 1) {
          # Multiple KW tests across taxa → Apply BH correction
          kw_p_values <- sapply(multi_group_tests, function(x) x$p_value)
          kw_p_adjusted <- p.adjust(kw_p_values, method = "BH")

          kw_names <- names(multi_group_tests)
          for(i in seq_along(kw_names)) {
            result_name <- kw_names[i]
            results[[result_name]]$p_adjusted <- kw_p_adjusted[i]
            results[[result_name]]$significant <- kw_p_adjusted[i] < input$alpha_level
          }
          cat("DEBUG: Applied BH correction to", length(multi_group_tests), "Kruskal-Wallis tests\n")
        } else if(length(multi_group_tests) == 1) {
          # Single KW test → NO BH correction needed
          result_name <- names(multi_group_tests)[1]
          results[[result_name]]$p_adjusted <- results[[result_name]]$p_value
          results[[result_name]]$significant <- results[[result_name]]$p_value < input$alpha_level
          cat("DEBUG: Single Kruskal-Wallis test - using raw p-value (no BH correction)\n")
        }

        cat("DEBUG: Final results summary by test_type:\n")
        cat("  - Correlation tests:", sum(sapply(results, function(x) !is.null(x$test_type) && x$test_type == "correlation")), "\n")
        cat("  - Two-group tests:", sum(sapply(results, function(x) !is.null(x$test_type) && x$test_type == "two_group")), "\n")
        cat("  - Multi-group tests:", sum(sapply(results, function(x) !is.null(x$test_type) && x$test_type == "multi_group")), "\n")

        return(results)

      }, error = function(e) {
        showNotification(paste("Taxonomy metadata statistics error:", e$message), type = "error")
        return(list())
      })
    })

    # Reactive to get significant taxa for SNV metric (used by SNV plot)
    significant_taxa_snv <- reactive({
      req(taxonomy_metadata_stats())

      tryCatch({
        stats_results <- taxonomy_metadata_stats()
        if(length(stats_results) > 0) {
          # Find results where significant is TRUE for SNV metric only
          cat("Debug: Checking", length(stats_results), "statistical results for SNV significance\n")
          significant_results <- stats_results[sapply(stats_results, function(x) {
            is_sig <- x$significant && x$metric == "SNVs per kbp"
            if (is_sig) {
              cat("Debug: Found significant SNV taxon:", x$taxon, "p_adj:", x$p_adjusted, "\n")
            }
            return(is_sig)
          })]
          if(length(significant_results) > 0) {
            sig_taxa <- unique(sapply(significant_results, function(x) x$taxon))
            cat("Debug: Found", length(sig_taxa), "significant SNV taxa:", paste(sig_taxa, collapse = ", "), "\n")
            return(sig_taxa)
          } else {
            cat("Debug: No significant SNV taxa found in statistical results\n")
            return(character(0))
          }
        } else {
          cat("Debug: No statistical results available\n")
          return(character(0))
        }
      }, error = function(e) {
        cat("Error in significant_taxa_snv:", e$message, "\n")
        return(character(0))
      })
    })

    # Reactive to get significant taxa for Nucleotide Diversity metric (used by Nucleotide Diversity plot)
    significant_taxa_nucl <- reactive({
      req(taxonomy_metadata_stats())

      tryCatch({
        stats_results <- taxonomy_metadata_stats()
        if(length(stats_results) > 0) {
          # Find results where significant is TRUE for Nucleotide Diversity metric only
          cat("Debug: Checking", length(stats_results), "statistical results for nucleotide diversity significance\n")
          significant_results <- stats_results[sapply(stats_results, function(x) {
            is_sig <- x$significant && x$metric == "Nucleotide Diversity"
            if (is_sig) {
              cat("Debug: Found significant nucleotide diversity taxon:", x$taxon, "p_adj:", x$p_adjusted, "\n")
            }
            return(is_sig)
          })]
          if(length(significant_results) > 0) {
            sig_taxa <- unique(sapply(significant_results, function(x) x$taxon))
            cat("Debug: Found", length(sig_taxa), "significant nucleotide diversity taxa:", paste(sig_taxa, collapse = ", "), "\n")
            return(sig_taxa)
          } else {
            cat("Debug: No significant nucleotide diversity taxa found in statistical results\n")
            return(character(0))
          }
        } else {
          cat("Debug: No statistical results available\n")
          return(character(0))
        }
      }, error = function(e) {
        cat("Error in significant_taxa_nucl:", e$message, "\n")
        return(character(0))
      })
    })

    # Reactive to get all significant taxa (union of both metrics, used for summary table)
    significant_taxa <- reactive({
      req(taxonomy_metadata_stats())

      tryCatch({
        stats_results <- taxonomy_metadata_stats()
        if(length(stats_results) > 0) {
          # Find results where significant is TRUE for any metric
          significant_results <- stats_results[sapply(stats_results, function(x) x$significant)]
          if(length(significant_results) > 0) {
            sig_taxa <- unique(sapply(significant_results, function(x) x$taxon))
            cat("Debug: Found", length(sig_taxa), "significant taxa (any metric):", paste(sig_taxa, collapse = ", "), "\n")
            return(sig_taxa)
          }
        }
        cat("Debug: No significant taxa found\n")
        return(character(0))
      }, error = function(e) {
        cat("Debug: Error in significant_taxa:", e$message, "\n")
        return(character(0))
      })
    })

    # Filtered statistics for "Show Significant Only" functionality (all significant taxa)
    filtered_taxonomy_metadata_stats <- reactive({
      req(taxonomy_metadata_stats())

      # Explicitly depend on the checkbox to trigger reactivity when it changes
      input$show_significant_only

      tryCatch({
        all_stats <- taxonomy_metadata_stats()

        if (input$show_significant_only) {
          # Filter to show only significant results (any metric)
          sig_taxa <- significant_taxa()
          if (length(sig_taxa) > 0) {
            filtered_stats <- all_stats[sapply(all_stats, function(x) x$taxon %in% sig_taxa)]
            cat("Debug: Filtered stats from", length(all_stats), "to", length(filtered_stats), "significant results\n")
            return(filtered_stats)
          } else {
            # No significant taxa found - return empty list
            cat("Debug: Show Significant Only is ON but no significant taxa found. Returning empty stats.\n")
            return(list())
          }
        } else {
          # Return all stats when not filtering
          return(all_stats)
        }

      }, error = function(e) {
        showNotification(paste("Filtered statistics error:", e$message), type = "error")
        return(taxonomy_metadata_stats())  # Fallback to all stats
      })
    })

    # Independent filtered data for SNV plot when "Show Significant Only" is checked
    taxonomy_plot_data_snv <- reactive({
      req(taxonomy_data_raw())

      # Explicitly depend on the checkbox to trigger reactivity when it changes
      input$show_significant_only

      tryCatch({
        data <- taxonomy_data_raw()

        # If not filtering, return baseline data
        if (!input$show_significant_only) {
          return(data)
        }

        # Get SNV-specific significant taxa
        sig_taxa <- significant_taxa_snv()
        cat("Debug: SNV plot - Show Significant Only =", input$show_significant_only, "\n")
        cat("Debug: SNV plot - Found", length(sig_taxa), "significant SNV taxa\n")

        if (length(sig_taxa) > 0) {
          cat("Debug: SNV plot - Significant taxa:", paste(sig_taxa, collapse = ", "), "\n")
        }

        # If no significant SNV taxa, return empty structured data
        if (length(sig_taxa) == 0 || is.null(data)) {
          cat("Debug: SNV plot - No significant taxa found. Returning empty dataset.\n")
          tax_col <- data$tax_level
          return(list(
            plot_data = data$plot_data[0, ],
            summary_data_nucl = data$summary_data_nucl[0, ],
            summary_data_snv = data$summary_data_snv[0, ],
            tax_level = data$tax_level,
            meta_col = data$meta_col,
            meta_groups = data$meta_groups,
            meta_type = data$meta_type,
            metric_names = data$metric_names,
            taxa_order = character(0)
          ))
        }

        # Filter data to include only SNV-significant taxa
        tax_col <- data$tax_level
        qualified_taxa <- intersect(data$summary_data_snv[[tax_col]], sig_taxa)

        cat("Debug: SNV plot - Qualified taxa after filtering:", length(qualified_taxa), "\n")

        # Filter plot data and summaries
        plot_data <- data$plot_data[data$plot_data[[tax_col]] %in% qualified_taxa, ]
        taxa_summary_nucl <- data$summary_data_nucl[data$summary_data_nucl[[tax_col]] %in% qualified_taxa, ]
        taxa_summary_snv <- data$summary_data_snv[data$summary_data_snv[[tax_col]] %in% qualified_taxa, ]

        # Update taxa ordering based on SNV means
        # Safety check: ensure the sorting column exists
        if ("mean_diversity" %in% colnames(taxa_summary_snv)) {
          taxa_order <- taxa_summary_snv %>%
            arrange(desc(mean_diversity)) %>%
            pull(.data[[tax_col]])
        } else {
          # Fallback: order alphabetically if sorting column not found
          taxa_order <- taxa_summary_snv %>%
            arrange(.data[[tax_col]]) %>%
            pull(.data[[tax_col]])
          cat("Warning: mean_diversity column not found, using alphabetical ordering\n")
        }

        return(list(
          plot_data = plot_data,
          summary_data_nucl = taxa_summary_nucl,
          summary_data_snv = taxa_summary_snv,
          tax_level = data$tax_level,
          meta_col = data$meta_col,
          meta_groups = data$meta_groups,
          meta_type = data$meta_type,
          metric_names = data$metric_names,
          taxa_order = taxa_order
        ))

      }, error = function(e) {
        showNotification(paste("SNV taxonomy filtering error:", e$message), type = "error")
        return(taxonomy_data_raw())
      })
    })

    # Independent filtered data for Nucleotide Diversity plot when "Show Significant Only" is checked
    taxonomy_plot_data_nucl <- reactive({
      req(taxonomy_data_raw())

      # Explicitly depend on the checkbox to trigger reactivity when it changes
      input$show_significant_only

      tryCatch({
        data <- taxonomy_data_raw()

        # If not filtering, return baseline data
        if (!input$show_significant_only) {
          return(data)
        }

        # Get Nucleotide Diversity-specific significant taxa
        sig_taxa <- significant_taxa_nucl()
        cat("Debug: Nucleotide Diversity plot - Show Significant Only =", input$show_significant_only, "\n")
        cat("Debug: Nucleotide Diversity plot - Found", length(sig_taxa), "significant nucleotide diversity taxa\n")

        if (length(sig_taxa) > 0) {
          cat("Debug: Nucleotide Diversity plot - Significant taxa:", paste(sig_taxa, collapse = ", "), "\n")
        }

        # If no significant nucleotide diversity taxa, return empty structured data
        if (length(sig_taxa) == 0 || is.null(data)) {
          cat("Debug: Nucleotide Diversity plot - No significant taxa found. Returning empty dataset.\n")
          tax_col <- data$tax_level
          return(list(
            plot_data = data$plot_data[0, ],
            summary_data_nucl = data$summary_data_nucl[0, ],
            summary_data_snv = data$summary_data_snv[0, ],
            tax_level = data$tax_level,
            meta_col = data$meta_col,
            meta_groups = data$meta_groups,
            meta_type = data$meta_type,
            metric_names = data$metric_names,
            taxa_order = character(0)
          ))
        }

        # Filter data to include only nucleotide diversity-significant taxa
        tax_col <- data$tax_level
        qualified_taxa <- intersect(data$summary_data_nucl[[tax_col]], sig_taxa)

        cat("Debug: Nucleotide Diversity plot - Qualified taxa after filtering:", length(qualified_taxa), "\n")

        # Filter plot data and summaries
        plot_data <- data$plot_data[data$plot_data[[tax_col]] %in% qualified_taxa, ]
        taxa_summary_nucl <- data$summary_data_nucl[data$summary_data_nucl[[tax_col]] %in% qualified_taxa, ]
        taxa_summary_snv <- data$summary_data_snv[data$summary_data_snv[[tax_col]] %in% qualified_taxa, ]

        # Update taxa ordering based on median diversity
        # Safety check: ensure the sorting column exists
        if ("median_diversity" %in% colnames(taxa_summary_nucl)) {
          taxa_order <- taxa_summary_nucl %>%
            arrange(desc(median_diversity)) %>%
            pull(.data[[tax_col]])
        } else {
          # Fallback: order alphabetically if sorting column not found
          taxa_order <- taxa_summary_nucl %>%
            arrange(.data[[tax_col]]) %>%
            pull(.data[[tax_col]])
          cat("Warning: median_diversity column not found, using alphabetical ordering\n")
        }

        return(list(
          plot_data = plot_data,
          summary_data_nucl = taxa_summary_nucl,
          summary_data_snv = taxa_summary_snv,
          tax_level = data$tax_level,
          meta_col = data$meta_col,
          meta_groups = data$meta_groups,
          meta_type = data$meta_type,
          metric_names = data$metric_names,
          taxa_order = taxa_order
        ))

      }, error = function(e) {
        showNotification(paste("Nucleotide Diversity taxonomy filtering error:", e$message), type = "error")
        return(taxonomy_data_raw())
      })
    })

    # Filtered data for SNV plot with plot-specific significant taxa filtering
    taxonomy_data_snv <- reactive({
      req(taxonomy_data_raw())

      # Explicitly depend on the checkbox to trigger reactivity when it changes
      input$show_significant_only

      tryCatch({
        data <- taxonomy_data_raw()

        # Apply "Show Significant Only" filtering for SNV metric
        if (input$show_significant_only) {
          sig_taxa <- significant_taxa_snv()
          cat("Debug: SNV plot - Show Significant Only =", input$show_significant_only, "\n")
          cat("Debug: SNV plot - Found", length(sig_taxa), "significant SNV taxa\n")
          if (length(sig_taxa) > 0) {
            cat("Debug: SNV plot - Significant taxa:", paste(sig_taxa, collapse = ", "), "\n")
          }

          if (length(sig_taxa) > 0 && !is.null(data)) {
            # Filter the baseline data to include only SNV-significant taxa
            tax_col <- data$tax_level
            # For SNV plot: only require taxa to be in SNV summary data AND be significant for SNV
            qualified_taxa <- intersect(data$summary_data_snv[[tax_col]], sig_taxa)

            # Also filter the nucleotide summary data to match the same taxa for consistency
            if (length(qualified_taxa) > 0) {
              # Ensure taxa are also in nucleotide summary data for consistent plotting
              qualified_taxa <- intersect(qualified_taxa, data$summary_data_nucl[[tax_col]])
            }

            cat("Debug: SNV plot - Qualified taxa after filtering:", length(qualified_taxa), "\n")
            cat("Debug: SNV plot - Original plot data rows:", nrow(data$plot_data), "\n")

            # Filter plot data
            plot_data <- data$plot_data[data$plot_data[[tax_col]] %in% qualified_taxa, ]
            cat("Debug: SNV plot - Filtered plot data rows:", nrow(plot_data), "\n")

            # Filter summaries
            taxa_summary_nucl <- data$summary_data_nucl[data$summary_data_nucl[[tax_col]] %in% qualified_taxa, ]
            taxa_summary_snv <- data$summary_data_snv[data$summary_data_snv[[tax_col]] %in% qualified_taxa, ]

            # Update taxa ordering based on filtered data
            # Safety check: ensure the sorting column exists
            if ("mean_diversity" %in% colnames(taxa_summary_snv)) {
              taxa_order <- taxa_summary_snv %>%
                arrange(desc(mean_diversity)) %>%
                pull(.data[[tax_col]])
            } else {
              # Fallback: order alphabetically if sorting column not found
              taxa_order <- taxa_summary_snv %>%
                arrange(.data[[tax_col]]) %>%
                pull(.data[[tax_col]])
              cat("Warning: mean_diversity column not found, using alphabetical ordering\n")
            }

            # Return filtered data
            return(list(
              plot_data = plot_data,
              summary_data_nucl = taxa_summary_nucl,
              summary_data_snv = taxa_summary_snv,
              tax_level = data$tax_level,
              meta_col = data$meta_col,
              meta_groups = data$meta_groups,
              meta_type = data$meta_type,
              metric_names = data$metric_names,
              taxa_order = taxa_order
            ))
          } else {
            # Checkbox is checked but no significant SNV taxa found - return empty dataset
            cat("Debug: SNV plot - Show Significant Only is ON but no significant SNV taxa found. Returning empty dataset.\n")
            tax_col <- data$tax_level

            # Return empty data structure with same format but zero rows
            return(list(
              plot_data = data$plot_data[0, ],
              summary_data_nucl = data$summary_data_nucl[0, ],
              summary_data_snv = data$summary_data_snv[0, ],
              tax_level = data$tax_level,
              meta_col = data$meta_col,
              meta_groups = data$meta_groups,
              meta_type = data$meta_type,
              metric_names = data$metric_names,
              taxa_order = character(0)  # Empty taxa order
            ))
          }
        }

        # Return baseline data when not filtering
        return(data)

      }, error = function(e) {
        showNotification(paste("SNV taxonomy filtering error:", e$message), type = "error")
        return(taxonomy_data_raw())  # Fallback to baseline data
      })
    })

    # Filtered data for Nucleotide Diversity plot with plot-specific significant taxa filtering
    taxonomy_data_nucl <- reactive({
      req(taxonomy_data_raw())

      # Explicitly depend on the checkbox to trigger reactivity when it changes
      input$show_significant_only

      tryCatch({
        data <- taxonomy_data_raw()

        # Apply "Show Significant Only" filtering for Nucleotide Diversity metric
        if (input$show_significant_only) {
          sig_taxa <- significant_taxa_nucl()
          cat("Debug: Nucleotide Diversity plot - Show Significant Only =", input$show_significant_only, "\n")
          cat("Debug: Nucleotide Diversity plot - Found", length(sig_taxa), "significant nucleotide diversity taxa\n")
          if (length(sig_taxa) > 0) {
            cat("Debug: Nucleotide Diversity plot - Significant taxa:", paste(sig_taxa, collapse = ", "), "\n")
          }

          if (length(sig_taxa) > 0 && !is.null(data)) {
            # Filter the baseline data to include only nucleotide diversity-significant taxa
            tax_col <- data$tax_level
            # For Nucleotide Diversity plot: only require taxa to be in nucleotide summary data AND be significant for nucleotide diversity
            qualified_taxa <- intersect(data$summary_data_nucl[[tax_col]], sig_taxa)

            # Also filter the SNV summary data to match the same taxa for consistency
            if (length(qualified_taxa) > 0) {
              # Ensure taxa are also in SNV summary data for consistent plotting
              qualified_taxa <- intersect(qualified_taxa, data$summary_data_snv[[tax_col]])
            }

            cat("Debug: Nucleotide Diversity plot - Qualified taxa after filtering:", length(qualified_taxa), "\n")
            cat("Debug: Nucleotide Diversity plot - Original plot data rows:", nrow(data$plot_data), "\n")

            # Filter plot data
            plot_data <- data$plot_data[data$plot_data[[tax_col]] %in% qualified_taxa, ]
            cat("Debug: Nucleotide Diversity plot - Filtered plot data rows:", nrow(plot_data), "\n")

            # Filter summaries
            taxa_summary_nucl <- data$summary_data_nucl[data$summary_data_nucl[[tax_col]] %in% qualified_taxa, ]
            taxa_summary_snv <- data$summary_data_snv[data$summary_data_snv[[tax_col]] %in% qualified_taxa, ]

            # Update taxa ordering based on filtered data
            taxa_order <- taxa_summary_nucl %>%
              arrange(desc(median_diversity)) %>%
              pull(.data[[tax_col]])

            # Return filtered data
            return(list(
              plot_data = plot_data,
              summary_data_nucl = taxa_summary_nucl,
              summary_data_snv = taxa_summary_snv,
              tax_level = data$tax_level,
              meta_col = data$meta_col,
              meta_groups = data$meta_groups,
              meta_type = data$meta_type,
              metric_names = data$metric_names,
              taxa_order = taxa_order
            ))
          } else {
            # Checkbox is checked but no significant nucleotide diversity taxa found - return empty dataset
            cat("Debug: Nucleotide Diversity plot - Show Significant Only is ON but no significant nucleotide diversity taxa found. Returning empty dataset.\n")
            tax_col <- data$tax_level

            # Return empty data structure with same format but zero rows
            return(list(
              plot_data = data$plot_data[0, ],
              summary_data_nucl = data$summary_data_nucl[0, ],
              summary_data_snv = data$summary_data_snv[0, ],
              tax_level = data$tax_level,
              meta_col = data$meta_col,
              meta_groups = data$meta_groups,
              meta_type = data$meta_type,
              metric_names = data$metric_names,
              taxa_order = character(0)  # Empty taxa order
            ))
          }
        }

        # Return baseline data when not filtering
        return(data)

      }, error = function(e) {
        showNotification(paste("Nucleotide Diversity taxonomy filtering error:", e$message), type = "error")
        return(taxonomy_data_raw())  # Fallback to baseline data
      })
    })


    # Taxonomy boxplot output - SNVs per kbp ONLY (individual plot)
    output$taxonomy_boxplot <- renderPlotly({
      req(taxonomy_plot_data_snv())

      tryCatch({
        # Get SNV data only
        tax_data_snv <- taxonomy_plot_data_snv()
        data_snv <- tax_data_snv$plot_data

        tax_col <- tax_data_snv$tax_level
        meta_col <- tax_data_snv$meta_col
        meta_type <- tax_data_snv$meta_type

        if (is.null(meta_col) || is.null(meta_type)) {
          raw_tax_data <- taxonomy_data_raw()
          if (is.null(meta_col)) meta_col <- raw_tax_data$meta_col
          if (is.null(meta_type)) meta_type <- raw_tax_data$meta_type
        }

        # Get statistical results
        stats_results <- taxonomy_metadata_stats()
        snv_results <- stats_results[sapply(stats_results, function(x) x$metric == "SNVs per kbp")]

        # Check if metadata is numerical
        is_numerical <- (meta_type == "numerical")

        if (is_numerical) {
          # For numerical variables: Create scatter plot
          data_snv$meta_numeric <- as.numeric(as.character(data_snv[[meta_col]]))
          n_sig_snv <- sum(sapply(snv_results, function(x) x$significant), na.rm = TRUE)

          p <- ggplot(data_snv, aes(x = meta_numeric, y = snvs_per_kbp_transformed,
                                    color = .data[[tax_col]], fill = .data[[tax_col]])) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
            labs(
              title = "SNVs per kbp by Taxa and Metadata",
              subtitle = paste0(n_sig_snv, " significant taxa (α=", input$alpha_level, ")",
                               if (input$show_significant_only) " | Filter: ACTIVE" else ""),
              x = meta_col,
              y = "SNVs per kbp"
            ) +
            theme_minimal() +
            theme(
              axis.text = element_text(size = 9),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 9)
            ) +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")

          return(ggplotly(p, height = input$plot_height, width = input$plot_width) %>%
                  layout(margin = list(l = input$plot_margin * 100,
                                      r = input$plot_margin * 100,
                                      t = input$plot_margin * 100,
                                      b = input$plot_margin * 100)))
        }

        # For categorical variables: Create boxplot
        if (!is.null(tax_data_snv$taxa_order)) {
          data_snv[[tax_col]] <- factor(data_snv[[tax_col]], levels = tax_data_snv$taxa_order)
        }

        # Apply SNV-specific taxa range filtering
        if (!is.null(input$snv_taxa_range) && length(input$snv_taxa_range) == 2 && !is.null(data_snv[[tax_col]])) {
          taxa_levels_snv <- levels(data_snv[[tax_col]])
          range_start <- input$snv_taxa_range[1]
          range_end <- min(input$snv_taxa_range[2], length(taxa_levels_snv))
          if (range_start <= range_end && range_start <= length(taxa_levels_snv)) {
            selected_taxa <- taxa_levels_snv[range_start:range_end]
            data_snv <- data_snv[data_snv[[tax_col]] %in% selected_taxa, ]
            data_snv[[tax_col]] <- factor(data_snv[[tax_col]], levels = selected_taxa)
          }
        }

        # Fix factor levels
        if (!is.null(meta_col) && meta_col %in% names(data_snv)) {
          data_snv[[meta_col]] <- droplevels(factor(data_snv[[meta_col]]))
        }
        dodge <- position_dodge2(width = 0.75, preserve = "single", padding = 0.1)

        # Get displayed taxa after filtering
        plot_taxa_snv <- unique(data_snv[[tax_col]])

        # Create annotation dataframe for SNVs per kbp
        annotation_df_snv <- data.frame()
        for (taxon_name in as.character(plot_taxa_snv)) {
          result_key <- paste0(taxon_name, "_snvs_per_kbp")
          if (result_key %in% names(snv_results)) {
            result <- snv_results[[result_key]]
            p_val <- result$p_adjusted
            if (!is.na(p_val) && p_val < input$alpha_level) {
              label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else ""
              if (label != "") {
                taxon_data <- data_snv[data_snv[[tax_col]] == taxon_name, ]
                y_pos <- max(taxon_data$snvs_per_kbp_transformed, na.rm = TRUE) * 1.1
                annotation_df_snv <- rbind(annotation_df_snv, data.frame(
                  taxon = taxon_name, y_position = y_pos, label = label, stringsAsFactors = FALSE
                ))
              }
            }
          }
        }

        # Create SNVs per kbp boxplot
        p <- ggplot(data_snv, aes_string(x = tax_col, y = "snvs_per_kbp_transformed", fill = meta_col)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.8, position = dodge) +
          geom_jitter(alpha = 0.3, size = 0.4, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))

        if (nrow(annotation_df_snv) > 0) {
          p <- p + geom_text(data = annotation_df_snv,
                            aes(x = taxon, y = y_position, label = label),
                            inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 0)
        }

        n_displayed_snv <- length(plot_taxa_snv)
        n_sig_snv <- nrow(annotation_df_snv)

        p <- p + labs(
            title = "SNVs per kbp by Taxa and Metadata Groups",
            subtitle = paste0(n_sig_snv, "/", n_displayed_snv, " significant (α=", input$alpha_level, ")",
                             if (input$show_significant_only) " | Filter: ACTIVE" else ""),
            x = "Taxa",
            y = "SNVs per kbp"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9)
          )

        ggplotly(p, height = input$plot_height, width = input$plot_width) %>%
          layout(margin = list(l = input$plot_margin * 100,
                              r = input$plot_margin * 100,
                              t = input$plot_margin * 100,
                              b = input$plot_margin * 100),
                  boxmode = "group")

      }, error = function(e) {
        showNotification(paste("SNV taxonomy boxplot error:", e$message), type = "error")
        return(plotly_empty())
      })
    })

    # Median Nucleotide Diversity plot output (UNIFIED STYLE - Boxplot)
    output$median_nucl_plot <- renderPlotly({
      req(taxonomy_plot_data_nucl())

      tryCatch({
        tax_data <- taxonomy_plot_data_nucl()
        data <- tax_data$plot_data
        tax_col <- tax_data$tax_level
        meta_col <- tax_data$meta_col

        if(nrow(data) == 0) {
          return(plot_ly() %>%
                   add_text(x = 0.5, y = 0.5, text = "No data available for median nucleotide diversity plot",
                            textfont = list(size = 16, color = "#6c757d")) %>%
                   layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }

        # Use the same taxa ordering as taxonomy_data() for synchronization
        if (!is.null(tax_data$taxa_order)) {
          data[[tax_col]] <- factor(data[[tax_col]], levels = tax_data$taxa_order)
        }

        # Apply Nucleotide Diversity-specific taxa range filtering
        if (!is.null(input$nucl_taxa_range) && length(input$nucl_taxa_range) == 2) {
          taxa_levels <- levels(data[[tax_col]])
          range_start <- input$nucl_taxa_range[1]
          range_end <- min(input$nucl_taxa_range[2], length(taxa_levels))

          if (range_start <= range_end && range_start <= length(taxa_levels)) {
            selected_taxa <- taxa_levels[range_start:range_end]
            data <- data[data[[tax_col]] %in% selected_taxa, ]
            data[[tax_col]] <- factor(data[[tax_col]], levels = selected_taxa)
          }
        }

        # Get statistical results from taxonomy_metadata_stats for manual annotations
        stats_results <- taxonomy_metadata_stats()
        nucl_results <- stats_results[sapply(stats_results, function(x) x$metric == "Nucleotide Diversity")]

        # Get displayed taxa after all filtering
        plot_taxa <- unique(data[[tax_col]])

        # Check if metadata is numerical or categorical
        meta_type <- tax_data$meta_type
        is_numerical <- (meta_type == "numerical")

        if (is_numerical) {
          # For numerical variables: Create combined scatter plot with colored lines per taxon
          # Convert metadata to numeric
          data$meta_numeric <- as.numeric(as.character(data[[meta_col]]))

          # Create correlation annotation dataframe for legend/table
          annotation_df <- data.frame()
          for (taxon_name in as.character(plot_taxa)) {
            result_key <- paste0(taxon_name, "_nucl_diversity")
            if (result_key %in% names(nucl_results)) {
              result <- nucl_results[[result_key]]
              rho <- result$statistic  # Correlation coefficient
              p_val <- result$p_adjusted

              if (!is.na(p_val)) {
                # Create annotation text with rho and significance
                sig_label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"

                annotation_df <- rbind(annotation_df, data.frame(
                  taxon = taxon_name,
                  rho = rho,
                  p_adjusted = p_val,
                  sig_label = sig_label,
                  significant = (p_val < input$alpha_level),
                  label = paste0(taxon_name, ": ρ=", round(rho, 3), " ", sig_label),
                  stringsAsFactors = FALSE
                ))
              }
            }
          }

          # Create combined scatter plot with all taxa, each in different color
          p <- ggplot(data, aes(x = meta_numeric, y = nucl_diversity_transformed,
                               color = .data[[tax_col]], fill = .data[[tax_col]])) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2)

          # Add text annotations showing correlation for each taxon
          if (nrow(annotation_df) > 0) {
            # Position annotations vertically along the right side
            annotation_df$y_pos <- seq(
              max(data$nucl_diversity_transformed, na.rm = TRUE) * 0.95,
              max(data$nucl_diversity_transformed, na.rm = TRUE) * 0.75,
              length.out = nrow(annotation_df)
            )
            annotation_df$x_pos <- max(data$meta_numeric, na.rm = TRUE) * 0.98

            p <- p + geom_text(data = annotation_df,
                              aes(x = x_pos, y = y_pos, label = label),
                              inherit.aes = FALSE, size = 3, fontface = "bold",
                              hjust = 1, vjust = 1, color = "black")
          }

        } else {
          # For categorical variables: Create boxplot (existing code)
          data[[meta_col]] <- droplevels(factor(data[[meta_col]]))
          dodge <- position_dodge2(width = 0.75, preserve = "single", padding = 0.1)

          # Create annotation dataframe based on taxonomy_metadata_stats results
          annotation_df <- data.frame()
          for (taxon_name in as.character(plot_taxa)) {
            result_key <- paste0(taxon_name, "_nucl_diversity")
            if (result_key %in% names(nucl_results)) {
              result <- nucl_results[[result_key]]
              p_val <- result$p_adjusted

              if (!is.na(p_val) && p_val < input$alpha_level) {
                label <- if (p_val < 0.001) {
                  "***"
                } else if (p_val < 0.01) {
                  "**"
                } else if (p_val < 0.05) {
                  "*"
                } else {
                  ""
                }

                if (label != "") {
                  taxon_data <- data[data[[tax_col]] == taxon_name, ]
                  y_pos <- max(taxon_data$nucl_diversity_transformed, na.rm = TRUE) * 1.1

                  annotation_df <- rbind(annotation_df, data.frame(
                    taxon = taxon_name,
                    y_position = y_pos,
                    label = label,
                    stringsAsFactors = FALSE
                  ))
                }
              }
            }
          }

          # Create boxplot
          p <- ggplot(data, aes_string(x = tax_col, y = "nucl_diversity_transformed", fill = meta_col)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.8, position = dodge) +
            geom_jitter(alpha = 0.3, size = 0.4, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))

          # Add manual significance annotations if any exist
          if (nrow(annotation_df) > 0) {
            p <- p + geom_text(data = annotation_df,
                              aes(x = taxon, y = y_position, label = label),
                              inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 0)
          }
        }

        # Create subtitle to show filtering status
        n_displayed_taxa <- length(plot_taxa)
        n_sig_shown <- nrow(annotation_df)

        if (is_numerical) {
          # Subtitle for numerical variables (correlation)
          n_sig_correlations <- sum(annotation_df$significant, na.rm = TRUE)

          subtitle_text <- if (input$show_significant_only) {
            all_sig_taxa <- significant_taxa_nucl()
            total_sig_taxa <- length(all_sig_taxa)
            paste0("Showing ", n_displayed_taxa, " significant taxa (",
                  n_sig_correlations, " with p < ", input$alpha_level, ")")
          } else {
            paste0("Showing ", n_displayed_taxa, " taxa (",
                  n_sig_correlations, " significant correlations, α=", input$alpha_level, ")")
          }

          # Add labs and theme for scatter plot
          p <- p + labs(
            title = paste("Nucleotide Diversity vs", meta_col),
            subtitle = subtitle_text,
            x = meta_col,
            y = "Nucleotide Diversity",
            color = "Taxon",
            fill = "Taxon"
          ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 9),
              plot.subtitle = element_text(size = 9, color = "#666666", face = "italic")
            ) +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")
        } else {
          # Subtitle for categorical variables (group comparison)
          subtitle_text <- if (input$show_significant_only) {
            all_sig_taxa <- significant_taxa_nucl()
            total_sig_taxa <- length(all_sig_taxa)
            paste0("Showing ", n_displayed_taxa, " significant taxa with ", n_sig_shown, " asterisks (out of ",
                  total_sig_taxa, " significant for nucleotide diversity, α=", input$alpha_level, ")")
          } else {
            paste0("Showing all ", n_displayed_taxa, " taxa with ", n_sig_shown, " significant (α=", input$alpha_level, ")")
          }

          # Add labs and theme for boxplot
          p <- p + labs(
            title = "Nucleotide Diversity comparison by Taxon",
            subtitle = subtitle_text,
            x = "Taxa",
            y = "Nucleotide Diversity"
          ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              legend.position = "right",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              plot.subtitle = element_text(size = 9, color = "#666666", face = "italic")
            )
        }

        ggplotly(p, height = 500) %>%
          layout(margin = list(l = 100, r = 100, t = 80, b = 150),
                  boxmode = "group")

      }, error = function(e) {
        showNotification(paste("Median nucl plot error:", e$message), type = "error")
        return(plotly_empty())
      })
    })

    # Taxonomy statistics table showing metadata group comparisons within taxa - SIMPLIFIED (NO POST-HOC)
    output$taxonomy_stats_table <- DT::renderDataTable({
      req(filtered_taxonomy_metadata_stats())

      tryCatch({
        # Use filtered stats to respect "Show Significant Only" checkbox
        stats_results <- filtered_taxonomy_metadata_stats()

        if(length(stats_results) == 0) {
          return(DT::datatable(data.frame(Message = "No statistical results available")))
        }

        # Convert results list to data frame, handling different test types
        stats_df <- do.call(rbind, lapply(names(stats_results), function(result_key) {
          result <- stats_results[[result_key]]

          # Get test type (default to "unknown" if not set)
          test_type <- if(!is.null(result$test_type)) result$test_type else "unknown"

          # Different columns based on test_type
          if(test_type == "correlation") {
            # === Numerical variable: Spearman Correlation ===
            data.frame(
              Metric = result$metric,
              Taxon = result$taxon,
              Test = result$test_name,
              Test_Type = "Correlation",
              Rho = round(result$rho, 3),  # Correlation coefficient
              N_Samples = result$n_samples,
              Variable_Range = paste0(round(result$variable_min, 1), " - ", round(result$variable_max, 1)),
              P_Value = format.pval(result$p_value, digits = 3),
              P_Adjusted = format.pval(result$p_adjusted, digits = 3),
              BH_Applied = "No",  # No BH correction for correlations
              Significant = result$significant,
              stringsAsFactors = FALSE
            )
          } else if(test_type == "two_group") {
            # === 2-group categorical: Mann-Whitney U ===
            data.frame(
              Metric = result$metric,
              Taxon = result$taxon,
              Test = result$test_name,
              Test_Type = "Two-Group",
              N_Groups = 2,
              Group1 = result$group1_name,
              N1 = result$n_group1,
              Mean1 = round(result$mean_group1, 3),
              Median1 = round(result$median_group1, 3),
              Group2 = result$group2_name,
              N2 = result$n_group2,
              Mean2 = round(result$mean_group2, 3),
              Median2 = round(result$median_group2, 3),
              P_Value = format.pval(result$p_value, digits = 3),
              P_Adjusted = format.pval(result$p_adjusted, digits = 3),
              BH_Applied = "No",  # No BH correction for single test
              Significant = result$significant,
              stringsAsFactors = FALSE
            )
          } else if(test_type == "multi_group") {
            # === 3+ group categorical: Kruskal-Wallis ===
            # Check if BH was applied (p_adjusted != p_value)
            bh_applied <- ifelse(abs(result$p_adjusted - result$p_value) > 1e-10, "Yes", "No")

            data.frame(
              Metric = result$metric,
              Taxon = result$taxon,
              Test = result$test_name,
              Test_Type = "Multi-Group",
              N_Groups = result$n_groups,
              N_Samples = result$n_samples,
              Groups = result$group_info,
              P_Value = format.pval(result$p_value, digits = 3),
              P_Adjusted = format.pval(result$p_adjusted, digits = 3),
              BH_Applied = bh_applied,  # Shows if BH was applied
              Significant = result$significant,
              stringsAsFactors = FALSE
            )
          } else {
            # Fallback for unknown test types
            data.frame(
              Metric = result$metric,
              Taxon = result$taxon,
              Test = result$test_name,
              Test_Type = "Unknown",
              P_Value = format.pval(result$p_value, digits = 3),
              P_Adjusted = format.pval(result$p_adjusted, digits = 3),
              Significant = result$significant,
              stringsAsFactors = FALSE
            )
          }
        }))

        # Add significance notation using NUMERIC p-values before formatting
        # Extract numeric p-values from stats_results for significance calculation
        numeric_p_adjusted <- sapply(stats_results, function(result) result$p_adjusted)
        stats_df$Significance_Notation <- sapply(numeric_p_adjusted, get_significance_notation)

        # Reorder columns for better display based on test type
        test_types <- unique(stats_df$Test_Type)

        if(length(test_types) == 1) {
          # Single test type - show only relevant columns
          if(test_types[1] == "Correlation") {
            # All correlation tests
            display_df <- stats_df %>%
              select(Metric, Taxon, Test, Rho, N_Samples, Variable_Range,
                     P_Value, P_Adjusted, BH_Applied, Significance_Notation, Significant) %>%
              arrange(Metric, Taxon)
          } else if(test_types[1] == "Two-Group") {
            # All two-group tests
            display_df <- stats_df %>%
              select(Metric, Taxon, Test, Group1, N1, Mean1, Median1,
                     Group2, N2, Mean2, Median2,
                     P_Value, P_Adjusted, BH_Applied, Significance_Notation, Significant) %>%
              arrange(Metric, Taxon)
          } else if(test_types[1] == "Multi-Group") {
            # All multi-group tests
            display_df <- stats_df %>%
              select(Metric, Taxon, Test, N_Groups, N_Samples, Groups,
                     P_Value, P_Adjusted, BH_Applied, Significance_Notation, Significant) %>%
              arrange(Metric, Taxon)
          } else {
            # Unknown type - show all
            display_df <- stats_df %>%
              select(Metric, Taxon, Test, Test_Type, everything()) %>%
              arrange(Metric, Taxon)
          }
        } else {
          # Mixed test types - include Test_Type column and key columns
          display_df <- stats_df %>%
            select(Metric, Taxon, Test, Test_Type, P_Value, P_Adjusted, BH_Applied,
                   Significance_Notation, Significant, everything()) %>%
            arrange(Metric, Taxon)
        }

        # Apply "Show Significant Only" filtering to the stats table
        if (input$show_significant_only) {
          display_df <- display_df %>% filter(Significant == TRUE)
        }

        dt_table <- DT::datatable(
          display_df,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            scrollY = "300px",
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel'),
            autoWidth = TRUE
          ),
          rownames = FALSE,
          extensions = c('Buttons')
        ) %>%
        DT::formatStyle(
          'Significant',
          backgroundColor = DT::styleEqual(TRUE, '#d4edda'),
          color = DT::styleEqual(TRUE, '#155724')
        ) %>%
        DT::formatStyle(
          'Significance_Notation',
          backgroundColor = DT::styleEqual(c("***", "**", "*", "."),
                                        c('#d73027', '#fc8d59', '#fee08b', '#d9ef8b')),
          color = DT::styleEqual(c("***", "**", "*", "."), c('white', 'white', 'black', 'black'))
        )

        # Add BH_Applied styling if the column exists
        if("BH_Applied" %in% names(display_df)) {
          dt_table <- dt_table %>%
            DT::formatStyle(
              'BH_Applied',
              backgroundColor = DT::styleEqual(c("Yes", "No"), c('#fff3cd', '#e7f3ff')),
              color = DT::styleEqual(c("Yes", "No"), c('#856404', '#004085')),
              fontWeight = DT::styleEqual("Yes", 'bold')
            )
        }

        return(dt_table)

      }, error = function(e) {
        showNotification(paste("Taxonomy stats table error:", e$message), type = "error")
        return(DT::datatable(data.frame(Error = "Failed to generate statistics table")))
      })
    })

    # SNVs vs Nucleotide Diversity Correlation Plot
    output$snv_nucl_correlation_plot <- renderPlotly({
      req(taxonomy_data_raw())

      tryCatch({
        tax_data <- taxonomy_data_raw()
        data <- tax_data$plot_data
        tax_col <- tax_data$tax_level

        if(is.null(data) || nrow(data) == 0) {
          return(plotly_empty("No data available for correlation analysis"))
        }

        # Get unique taxa
        unique_taxa <- unique(data[[tax_col]])

        # Limit to top taxa if too many
        if(length(unique_taxa) > 20) {
          # Use top 20 taxa by sample count
          top_taxa <- data %>%
            group_by(.data[[tax_col]]) %>%
            summarise(n = n(), .groups = "drop") %>%
            arrange(desc(n)) %>%
            head(20) %>%
            pull(.data[[tax_col]])

          data <- data %>% filter(.data[[tax_col]] %in% top_taxa)
          unique_taxa <- top_taxa
        }

        # Create scatter plot with regression lines per taxon
        p <- ggplot(data, aes(x = snvs_per_kbp_transformed,
                             y = nucl_diversity_transformed,
                             color = .data[[tax_col]],
                             fill = .data[[tax_col]])) +
          geom_point(alpha = 0.5, size = 2) +
          geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
          scale_color_brewer(palette = "Set1", name = tax_col) +
          scale_fill_brewer(palette = "Set1", name = tax_col) +
          labs(
            title = "SNVs per kbp vs Nucleotide Diversity by Taxon",
            x = "SNVs per kbp",
            y = "Nucleotide Diversity"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            legend.position = "right",
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10)
          )

        # Convert to plotly
        ggplotly(p, height = 600) %>%
          layout(
            legend = list(
              orientation = "v",
              y = 0.5,
              yanchor = "middle"
            )
          )

      }, error = function(e) {
        cat("ERROR in snv_nucl_correlation_plot:", e$message, "\n")
        showNotification(paste("Correlation plot error:", e$message), type = "error")
        return(plotly_empty("Error generating correlation plot"))
      })
    })

    # SNVs vs Nucleotide Diversity Correlation Statistics Table
    output$snv_nucl_correlation_table <- DT::renderDataTable({
      req(taxonomy_data_raw())

      tryCatch({
        tax_data <- taxonomy_data_raw()
        data <- tax_data$plot_data
        tax_col <- tax_data$tax_level

        if(is.null(data) || nrow(data) == 0) {
          return(DT::datatable(data.frame(Message = "No data available for correlation analysis")))
        }

        # Calculate correlation for each taxon
        unique_taxa <- unique(data[[tax_col]])
        correlation_results <- list()

        for(taxon in unique_taxa) {
          taxon_data <- data[data[[tax_col]] == taxon, ]

          if(nrow(taxon_data) < 3) next  # Need at least 3 points for correlation

          # Remove NA values
          valid_data <- taxon_data[!is.na(taxon_data$snvs_per_kbp_transformed) &
                                   !is.na(taxon_data$nucl_diversity_transformed), ]

          if(nrow(valid_data) < 3) next

          # Calculate Spearman correlation
          cor_test <- tryCatch({
            cor.test(valid_data$snvs_per_kbp_transformed,
                    valid_data$nucl_diversity_transformed,
                    method = "spearman", exact = FALSE)
          }, error = function(e) NULL)

          if(!is.null(cor_test) && !is.na(cor_test$p.value)) {
            correlation_results[[taxon]] <- data.frame(
              Taxon = taxon,
              N_Samples = nrow(valid_data),
              Rho = as.numeric(cor_test$estimate),
              P_Value = cor_test$p.value,
              Significant = cor_test$p.value < 0.05,
              stringsAsFactors = FALSE
            )
          }
        }

        if(length(correlation_results) == 0) {
          return(DT::datatable(data.frame(Message = "No correlation results available")))
        }

        # Combine results
        correlation_df <- do.call(rbind, correlation_results)

        # Add significance notation
        correlation_df$Significance <- sapply(correlation_df$P_Value, get_significance_notation)

        # Format for display
        display_df <- correlation_df %>%
          mutate(
            Rho = round(Rho, 4),
            P_Value = format.pval(P_Value, digits = 3),
            Direction = ifelse(Rho > 0, "Positive", "Negative")
          ) %>%
          arrange(desc(abs(Rho))) %>%  # Sort by absolute correlation strength
          select(Taxon, N_Samples, Rho, Direction, P_Value, Significance, Significant)

        # Create DataTable
        DT::datatable(
          display_df,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            scrollY = "350px",
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel'),
            autoWidth = TRUE
          ),
          rownames = FALSE,
          extensions = c('Buttons'),
          caption = "Spearman Correlation: SNVs per kbp vs Nucleotide Diversity"
        ) %>%
        DT::formatStyle(
          'Significant',
          backgroundColor = DT::styleEqual(TRUE, '#d4edda'),
          color = DT::styleEqual(TRUE, '#155724')
        ) %>%
        DT::formatStyle(
          'Significance',
          backgroundColor = DT::styleEqual(c("***", "**", "*", "."),
                                        c('#d73027', '#fc8d59', '#fee08b', '#d9ef8b')),
          color = DT::styleEqual(c("***", "**", "*", "."),
                                c('white', 'white', 'black', 'black'))
        ) %>%
        DT::formatStyle(
          'Direction',
          backgroundColor = DT::styleEqual(c("Positive", "Negative"),
                                        c('#d1ecf1', '#f8d7da')),
          color = DT::styleEqual(c("Positive", "Negative"),
                              c('#0c5460', '#721c24'))
        )

      }, error = function(e) {
        cat("ERROR in snv_nucl_correlation_table:", e$message, "\n")
        showNotification(paste("Correlation table error:", e$message), type = "error")
        return(DT::datatable(data.frame(Error = "Failed to generate correlation table")))
      })
    })

    # Taxonomy summary table - updated for both metrics
    output$taxonomy_summary_table <- DT::renderDataTable({
      req(taxonomy_data())

      tryCatch({
        tax_data <- taxonomy_data()
        summary_data_nucl <- tax_data$summary_data_nucl
        summary_data_snv <- tax_data$summary_data_snv
        tax_col <- tax_data$tax_level

        if(is.null(summary_data_nucl) || is.null(summary_data_snv) ||
           nrow(summary_data_nucl) == 0 || nrow(summary_data_snv) == 0) {
          return(DT::datatable(data.frame(Message = "No summary data available")))
        }

        # Join summary data for both metrics
        summary_combined <- summary_data_nucl %>%
          select(all_of(tax_col), n_samples, median_diversity, mean_diversity,
                 q25_diversity, q75_diversity, min_diversity, max_diversity) %>%
          rename(
            N_Samples = n_samples,
            Median_Nucl = median_diversity,
            Mean_Nucl = mean_diversity,
            Q25_Nucl = q25_diversity,
            Q75_Nucl = q75_diversity,
            Min_Nucl = min_diversity,
            Max_Nucl = max_diversity
          ) %>%
          left_join(
            summary_data_snv %>%
              select(all_of(tax_col), median_diversity, mean_diversity,
                     q25_diversity, q75_diversity, min_diversity, max_diversity) %>%
              rename(
                Median_SNV = median_diversity,
                Mean_SNV = mean_diversity,
                Q25_SNV = q25_diversity,
                Q75_SNV = q75_diversity,
                Min_SNV = min_diversity,
                Max_SNV = max_diversity
              ),
            by = tax_col
          ) %>%
          mutate(
            IQR_Nucl = round(Q75_Nucl - Q25_Nucl, 4),
            IQR_SNV = round(Q75_SNV - Q25_SNV, 4),
            Median_Nucl = round(Median_Nucl, 4),
            Mean_Nucl = round(Mean_Nucl, 4),
            Q25_Nucl = round(Q25_Nucl, 4),
            Q75_Nucl = round(Q75_Nucl, 4),
            Min_Nucl = round(Min_Nucl, 4),
            Max_Nucl = round(Max_Nucl, 4),
            Median_SNV = round(Median_SNV, 4),
            Mean_SNV = round(Mean_SNV, 4),
            Q25_SNV = round(Q25_SNV, 4),
            Q75_SNV = round(Q75_SNV, 4),
            Min_SNV = round(Min_SNV, 4),
            Max_SNV = round(Max_SNV, 4)
          ) %>%
          arrange(desc(Median_Nucl))

        # Reorder columns for better display
        display_df <- summary_combined %>%
          select(all_of(tax_col), N_Samples,
                 Median_Nucl, Mean_Nucl, Q25_Nucl, Q75_Nucl, IQR_Nucl, Min_Nucl, Max_Nucl,
                 Median_SNV, Mean_SNV, Q25_SNV, Q75_SNV, IQR_SNV, Min_SNV, Max_SNV)

        colnames(display_df) <- c(
          "Taxon", "N_Samples",
          "Median_Nucl", "Mean_Nucl", "Q25_Nucl", "Q75_Nucl", "IQR_Nucl", "Min_Nucl", "Max_Nucl",
          "Median_SNV", "Mean_SNV", "Q25_SNV", "Q75_SNV", "IQR_SNV", "Min_SNV", "Max_SNV"
        )

        DT::datatable(
          display_df,
          options = list(
            pageLength = 20,
            scrollX = TRUE,
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
          ),
          rownames = FALSE,
          extensions = 'Buttons'
        ) %>%
        DT::formatStyle(
          'N_Samples',
          backgroundColor = DT::styleInterval(c(5, 10, 20), c('#fee5d9', '#fcae91', '#fb6a4a', '#de2d26'))
        )

      }, error = function(e) {
        showNotification(paste("Taxonomy summary table error:", e$message), type = "error")
        return(DT::datatable(data.frame(Error = "Failed to generate summary table")))
      })
    })

    # ========== Download Modal Handlers ==========

    # Overview Diversity Distribution ggplot
    overview_diversity_ggplot <- reactive({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      req(nrow(filtered_data) > 0)

      p <- ggplot(filtered_data, aes(x = nucl_diversity_transformed)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "#2FA4E7", alpha = 0.7) +
        geom_density(color = "darkblue", linewidth = 1) +
        labs(title = "Nucleotide Diversity Distribution", x = "Nucleotide Diversity", y = "Density") +
        theme_minimal()
      apply_plot_style(p, input)
    })

    # Overview Quality ggplot
    overview_quality_ggplot <- reactive({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      req(nrow(filtered_data) > 0)

      p <- ggplot(filtered_data, aes(x = coverage, y = breadth)) +
        geom_point(alpha = 0.5, color = "#28a745") +
        labs(title = "Coverage vs Breadth", x = "Coverage", y = "Breadth") +
        theme_minimal()
      apply_plot_style(p, input)
    })

    # Metadata Comparison ggplot
    metadata_comparison_ggplot <- reactive({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      meta_var <- data$meta_var
      var_type <- data$var_type
      req(nrow(filtered_data) > 0)

      if (var_type == "categorical") {
        p <- ggplot(filtered_data, aes_string(x = meta_var, y = "nucl_diversity_transformed", fill = meta_var)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.8) +
          geom_jitter(width = 0.2, alpha = 0.3, size = input$point_size) +
          labs(title = paste("Nucleotide Diversity by", meta_var), x = meta_var, y = "Nucleotide Diversity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      } else {
        p <- ggplot(filtered_data, aes_string(x = meta_var, y = "nucl_diversity_transformed")) +
          geom_point(alpha = 0.6, color = "#2FA4E7", size = input$point_size) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(title = paste("Nucleotide Diversity vs", meta_var), x = meta_var, y = "Nucleotide Diversity") +
          theme_minimal()
      }
      apply_plot_style(p, input)
    })

    # Taxonomy Boxplot ggplot - SNVs per kbp ONLY (matches screen plot for download)
    taxonomy_boxplot_ggplot <- reactive({
      req(taxonomy_plot_data_snv())

      tryCatch({
        # Get SNV data only (same as screen plot)
        tax_data_snv <- taxonomy_plot_data_snv()
        data_snv <- tax_data_snv$plot_data

        tax_col <- tax_data_snv$tax_level
        meta_col <- tax_data_snv$meta_col
        meta_type <- tax_data_snv$meta_type

        if (is.null(meta_col) || is.null(meta_type)) {
          raw_tax_data <- taxonomy_data_raw()
          if (is.null(meta_col)) meta_col <- raw_tax_data$meta_col
          if (is.null(meta_type)) meta_type <- raw_tax_data$meta_type
        }

        # Get statistical results
        stats_results <- taxonomy_metadata_stats()
        snv_results <- stats_results[sapply(stats_results, function(x) x$metric == "SNVs per kbp")]

        # Check if metadata is numerical
        is_numerical <- (meta_type == "numerical")

        if (is_numerical) {
          # For numerical variables: Create scatter plot
          data_snv$meta_numeric <- as.numeric(as.character(data_snv[[meta_col]]))
          n_sig_snv <- sum(sapply(snv_results, function(x) x$significant), na.rm = TRUE)

          p <- ggplot(data_snv, aes(x = meta_numeric, y = snvs_per_kbp_transformed,
                                    color = .data[[tax_col]], fill = .data[[tax_col]])) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
            labs(
              title = "SNVs per kbp by Taxa and Metadata",
              subtitle = paste0(n_sig_snv, " significant taxa (α=", input$alpha_level, ")",
                               if (input$show_significant_only) " | Filter: ACTIVE" else ""),
              x = meta_col,
              y = "SNVs per kbp"
            ) +
            theme_minimal() +
            theme(
              axis.text = element_text(size = 9),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 9)
            ) +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")

          return(apply_plot_style(p, input))
        }

        # For categorical variables: Create boxplot
        if (!is.null(tax_data_snv$taxa_order)) {
          data_snv[[tax_col]] <- factor(data_snv[[tax_col]], levels = tax_data_snv$taxa_order)
        }

        # Apply SNV-specific taxa range filtering
        if (!is.null(input$snv_taxa_range) && length(input$snv_taxa_range) == 2 && !is.null(data_snv[[tax_col]])) {
          taxa_levels_snv <- levels(data_snv[[tax_col]])
          range_start <- input$snv_taxa_range[1]
          range_end <- min(input$snv_taxa_range[2], length(taxa_levels_snv))
          if (range_start <= range_end && range_start <= length(taxa_levels_snv)) {
            selected_taxa <- taxa_levels_snv[range_start:range_end]
            data_snv <- data_snv[data_snv[[tax_col]] %in% selected_taxa, ]
            data_snv[[tax_col]] <- factor(data_snv[[tax_col]], levels = selected_taxa)
          }
        }

        # Fix factor levels
        if (!is.null(meta_col) && meta_col %in% names(data_snv)) {
          data_snv[[meta_col]] <- droplevels(factor(data_snv[[meta_col]]))
        }
        dodge <- position_dodge2(width = 0.75, preserve = "single", padding = 0.1)

        # Get displayed taxa after filtering
        plot_taxa_snv <- unique(data_snv[[tax_col]])

        # Create annotation dataframe for SNVs per kbp
        annotation_df_snv <- data.frame()
        for (taxon_name in as.character(plot_taxa_snv)) {
          result_key <- paste0(taxon_name, "_snvs_per_kbp")
          if (result_key %in% names(snv_results)) {
            result <- snv_results[[result_key]]
            p_val <- result$p_adjusted
            if (!is.na(p_val) && p_val < input$alpha_level) {
              label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else ""
              if (label != "") {
                taxon_data <- data_snv[data_snv[[tax_col]] == taxon_name, ]
                y_pos <- max(taxon_data$snvs_per_kbp_transformed, na.rm = TRUE) * 1.1
                annotation_df_snv <- rbind(annotation_df_snv, data.frame(
                  taxon = taxon_name, y_position = y_pos, label = label, stringsAsFactors = FALSE
                ))
              }
            }
          }
        }

        # Create SNVs per kbp boxplot
        p <- ggplot(data_snv, aes_string(x = tax_col, y = "snvs_per_kbp_transformed", fill = meta_col)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.8, position = dodge) +
          geom_jitter(alpha = 0.3, size = 0.4, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))

        if (nrow(annotation_df_snv) > 0) {
          p <- p + geom_text(data = annotation_df_snv,
                            aes(x = taxon, y = y_position, label = label),
                            inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 0)
        }

        n_displayed_snv <- length(plot_taxa_snv)
        n_sig_snv <- nrow(annotation_df_snv)

        p <- p + labs(
            title = "SNVs per kbp by Taxa and Metadata Groups",
            subtitle = paste0(n_sig_snv, "/", n_displayed_snv, " significant (α=", input$alpha_level, ")",
                             if (input$show_significant_only) " | Filter: ACTIVE" else ""),
            x = "Taxa",
            y = "SNVs per kbp"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9)
          )

        apply_plot_style(p, input)

      }, error = function(e) {
        cat("Error in taxonomy_boxplot_ggplot:", e$message, "\n")
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
          theme_void()
      })
    })

    # Median Nucleotide Diversity ggplot - matches screen plot for download
    median_nucl_ggplot <- reactive({
      req(taxonomy_plot_data_nucl())

      tryCatch({
        tax_data <- taxonomy_plot_data_nucl()
        data <- tax_data$plot_data
        tax_col <- tax_data$tax_level
        meta_col <- tax_data$meta_col
        meta_type <- tax_data$meta_type

        if (is.null(meta_col) || is.null(meta_type)) {
          raw_tax_data <- taxonomy_data_raw()
          if (is.null(meta_col)) meta_col <- raw_tax_data$meta_col
          if (is.null(meta_type)) meta_type <- raw_tax_data$meta_type
        }

        # Use taxa ordering
        if (!is.null(tax_data$taxa_order)) {
          data[[tax_col]] <- factor(data[[tax_col]], levels = tax_data$taxa_order)
        }

        # Apply Nucleotide Diversity-specific taxa range filtering
        if (!is.null(input$nucl_taxa_range) && length(input$nucl_taxa_range) == 2 && !is.null(data[[tax_col]])) {
          taxa_levels <- levels(data[[tax_col]])
          range_start <- input$nucl_taxa_range[1]
          range_end <- min(input$nucl_taxa_range[2], length(taxa_levels))
          if (range_start <= range_end && range_start <= length(taxa_levels)) {
            selected_taxa <- taxa_levels[range_start:range_end]
            data <- data[data[[tax_col]] %in% selected_taxa, ]
            data[[tax_col]] <- factor(data[[tax_col]], levels = selected_taxa)
          }
        }

        # Get statistical results
        stats_results <- taxonomy_metadata_stats()
        nucl_results <- stats_results[sapply(stats_results, function(x) x$metric == "Nucleotide Diversity")]

        plot_taxa <- unique(data[[tax_col]])
        is_numerical <- (meta_type == "numerical")

        if (is_numerical) {
          # For numerical variables: Create scatter plot
          data$meta_numeric <- as.numeric(as.character(data[[meta_col]]))
          n_sig_nucl <- sum(sapply(nucl_results, function(x) x$significant), na.rm = TRUE)

          p <- ggplot(data, aes(x = meta_numeric, y = nucl_diversity_transformed,
                               color = .data[[tax_col]], fill = .data[[tax_col]])) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
            labs(
              title = paste("Nucleotide Diversity vs", meta_col),
              subtitle = paste0(n_sig_nucl, " significant taxa (α=", input$alpha_level, ")",
                               if (input$show_significant_only) " | Filter: ACTIVE" else ""),
              x = meta_col,
              y = "Nucleotide Diversity"
            ) +
            theme_minimal() +
            theme(
              axis.text = element_text(size = 9),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 9)
            ) +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")

          return(apply_plot_style(p, input))
        }

        # For categorical variables: Create boxplot
        if (!is.null(meta_col) && meta_col %in% names(data)) {
          data[[meta_col]] <- droplevels(factor(data[[meta_col]]))
        }
        dodge <- position_dodge2(width = 0.75, preserve = "single", padding = 0.1)

        # Create annotation dataframe
        annotation_df <- data.frame()
        for (taxon_name in as.character(plot_taxa)) {
          result_key <- paste0(taxon_name, "_nucl_diversity")
          if (result_key %in% names(nucl_results)) {
            result <- nucl_results[[result_key]]
            p_val <- result$p_adjusted
            if (!is.na(p_val) && p_val < input$alpha_level) {
              label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else ""
              if (label != "") {
                taxon_data <- data[data[[tax_col]] == taxon_name, ]
                y_pos <- max(taxon_data$nucl_diversity_transformed, na.rm = TRUE) * 1.1
                annotation_df <- rbind(annotation_df, data.frame(
                  taxon = taxon_name, y_position = y_pos, label = label, stringsAsFactors = FALSE
                ))
              }
            }
          }
        }

        # Create boxplot
        p <- ggplot(data, aes_string(x = tax_col, y = "nucl_diversity_transformed", fill = meta_col)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.8, position = dodge) +
          geom_jitter(alpha = 0.3, size = 0.4, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))

        if (nrow(annotation_df) > 0) {
          p <- p + geom_text(data = annotation_df,
                            aes(x = taxon, y = y_position, label = label),
                            inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 0)
        }

        n_displayed <- length(plot_taxa)
        n_sig <- nrow(annotation_df)

        p <- p + labs(
            title = "Nucleotide Diversity by Taxa and Metadata Groups",
            subtitle = paste0(n_sig, "/", n_displayed, " significant (α=", input$alpha_level, ")",
                             if (input$show_significant_only) " | Filter: ACTIVE" else ""),
            x = "Taxa",
            y = "Nucleotide Diversity"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9)
          )

        apply_plot_style(p, input)

      }, error = function(e) {
        cat("Error in median_nucl_ggplot:", e$message, "\n")
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
          theme_void()
      })
    })

    # SNV-Nucleotide Correlation ggplot - matches screen plot for download
    snv_nucl_correlation_ggplot <- reactive({
      req(taxonomy_data_raw())

      tryCatch({
        tax_data <- taxonomy_data_raw()
        data <- tax_data$plot_data
        tax_col <- tax_data$tax_level

        if (is.null(data) || nrow(data) == 0) {
          return(ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "No data available") +
            theme_void())
        }

        # Get unique taxa
        unique_taxa <- unique(data[[tax_col]])

        # Limit to top 20 taxa if too many (same as screen)
        if (length(unique_taxa) > 20) {
          top_taxa <- data %>%
            group_by(.data[[tax_col]]) %>%
            summarise(n = n(), .groups = "drop") %>%
            arrange(desc(n)) %>%
            head(20) %>%
            pull(.data[[tax_col]])

          data <- data %>% filter(.data[[tax_col]] %in% top_taxa)
        }

        # Create scatter plot with regression lines per taxon (same as screen)
        p <- ggplot(data, aes(x = snvs_per_kbp_transformed,
                             y = nucl_diversity_transformed,
                             color = .data[[tax_col]],
                             fill = .data[[tax_col]])) +
          geom_point(alpha = 0.5, size = 2) +
          geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
          scale_color_brewer(palette = "Set1", name = tax_col) +
          scale_fill_brewer(palette = "Set1", name = tax_col) +
          labs(
            title = "SNVs per kbp vs Nucleotide Diversity by Taxon",
            x = "SNVs per kbp",
            y = "Nucleotide Diversity"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            legend.position = "right",
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10)
          )

        apply_plot_style(p, input)

      }, error = function(e) {
        cat("Error in snv_nucl_correlation_ggplot:", e$message, "\n")
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
          theme_void()
      })
    })

    # Coverage Breadth ggplot
    coverage_breadth_ggplot <- reactive({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      req(nrow(filtered_data) > 0)

      p <- ggplot(filtered_data, aes(x = coverage, y = breadth)) +
        geom_point(alpha = 0.5, color = "#28a745", size = input$point_size) +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(title = "Coverage vs Breadth", x = "Coverage", y = "Breadth") +
        theme_minimal()
      apply_plot_style(p, input)
    })

    # Genome Length ggplot
    genome_length_ggplot <- reactive({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      req(nrow(filtered_data) > 0)

      p <- ggplot(filtered_data, aes(x = genome_length / 1e6)) +
        geom_histogram(bins = 30, fill = "#6f42c1", alpha = 0.7) +
        labs(title = "Genome Length Distribution", x = "Genome Length (Mb)", y = "Count") +
        theme_minimal()
      apply_plot_style(p, input)
    })

    # Quality Metrics ggplot
    quality_metrics_ggplot <- reactive({
      data <- processed_data()
      req(data)
      filtered_data <- data$filtered_data
      req(nrow(filtered_data) > 0)

      p <- ggplot(filtered_data, aes(x = coverage, y = nucl_diversity_transformed)) +
        geom_point(alpha = 0.5, color = "#fd7e14", size = input$point_size) +
        geom_smooth(method = "loess", se = TRUE, color = "blue") +
        labs(title = "Coverage vs Nucleotide Diversity", x = "Coverage", y = "Nucleotide Diversity") +
        theme_minimal()
      apply_plot_style(p, input)
    })

    # Setup download modal handlers
    setup_download_modal_handler(input, output, session, "overview_diversity", overview_diversity_ggplot, "nucdiv_overview_diversity")
    setup_download_modal_handler(input, output, session, "overview_quality", overview_quality_ggplot, "nucdiv_overview_quality")
    setup_download_modal_handler(input, output, session, "metadata_comparison", metadata_comparison_ggplot, "nucdiv_metadata_comparison")
    setup_download_modal_handler(input, output, session, "taxonomy_boxplot", taxonomy_boxplot_ggplot, "nucdiv_taxonomy_boxplot")
    setup_download_modal_handler(input, output, session, "median_nucl", median_nucl_ggplot, "nucdiv_median_nucl")
    setup_download_modal_handler(input, output, session, "snv_nucl_correlation", snv_nucl_correlation_ggplot, "nucdiv_snv_correlation")
    setup_download_modal_handler(input, output, session, "coverage_breadth", coverage_breadth_ggplot, "nucdiv_coverage_breadth")
    setup_download_modal_handler(input, output, session, "genome_length", genome_length_ggplot, "nucdiv_genome_length")
    setup_download_modal_handler(input, output, session, "quality_metrics", quality_metrics_ggplot, "nucdiv_quality_metrics")

  })
}
