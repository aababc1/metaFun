# =============================================================================
# Network Analysis Shiny App - SIF Container Version
# Based on ggClusterNet2 for group-wise network comparison
# =============================================================================

library(shiny)
library(shinyWidgets)
library(shinybusy)
library(shinyjs)
library(bslib)
library(plotly)
library(DT)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(igraph)
library(Cairo)
library(tidyr)
library(RColorBrewer)
library(bsicons)
library(future)
library(promises)
library(callr)  # For killable background processes

# Optional: ggClusterNet for advanced features (not required for basic analysis)
# library(ggClusterNet)
# Optional: microeco for additional data transformations (not required)
# library(microeco)

# Future configuration for async processing
future::plan(multisession)

# Source modules - Use current working directory (set by container)
app_dir <- getwd()
source(file.path(app_dir, "modules/data_loading_module.R"))
source(file.path(app_dir, "modules/filtering_module.R"))
source(file.path(app_dir, "modules/network_analysis_module.R"))
source(file.path(app_dir, "modules/network_visualization_module.R"))
source(file.path(app_dir, "modules/influential_nodes_module.R"))
source(file.path(app_dir, "modules/network_robustness_module.R"))

# Helper functions
source(file.path(app_dir, "helper/data_conversion.R"))
source(file.path(app_dir, "helper/network_utils.R"))
source(file.path(app_dir, "helper/plot_download_helper.R"))

# =============================================================================
# UI Definition
# =============================================================================

ui <- page_navbar(
  header = tagList(
    shinyjs::useShinyjs(),
    shinybusy::use_busy_spinner(spin = "fading-circle", position = "top-right"),
    tags$head(
      tags$meta(name="viewport", content="width=device-width, initial-scale=1.0"),
      tags$style(HTML("
        .sidebar {
          min-width: 250px;
          width: 250px;
          padding: 0 !important;
          margin: 0 !important;
          flex: 0 0 250px !important;
          -webkit-flex: 0 0 250px !important;
          flex-shrink: 0;
          position: relative;
          background: transparent !important;
        }

        .sidebar-content {
          max-height: 100%;
          overflow-y: auto;
          position: sticky;
          top: 0;
          margin: 0;
          padding: 0;
          width: 100%;
          background: #f8f9fa;
        }

        .nav-underline {
          margin-bottom: 0 !important;
          padding-bottom: 0 !important;
        }

        .page-navbar {
          margin-top: 0 !important;
        }

        .nav-panel {
          padding-top: 0 !important;
        }

        .card {
          margin-top: 0 !important;
        }

        .card-body {
          flex: 1;
          display: flex;
          padding: 0 !important;
          flex-direction: column;
          overflow: hidden;
        }

        .main-card {
          display: flex;
          flex-direction: column;
          height: calc(100vh - 60px);
          overflow-y: auto;
          padding: 0;
          margin: 0;
        }

        .inner-card {
          flex: 0 0 auto;
          margin-bottom: 20px;
          min-height: fit-content;
        }

        .inner-card .card-body {
          overflow-y: auto;
          height: auto !important;
          display: block;
          min-height: fit-content;
        }

        .main-card > .card-body {
          flex: 1;
          overflow-y: auto;
        }

        .main-card .inner-card {
          overflow: visible;
        }

        .blink {
          animation: blinker 1.8s linear infinite;
        }

        @keyframes blinker {
          50% { opacity: 0; }
        }

        #load_data {
          width: 100%;
          font-size: 16px;
          padding: 10px;
          margin: 10px 0;
          z-index: 1000;
          position: relative;
        }

        #master_run_analysis {
          background-color: #2FA4E7;
          color: white;
          border-color: #2FA4E7;
        }

        .network-analysis {
          color: #2FA4E7 !important;
        }

        [data-bs-theme='light'] .nav-link.active .network-analysis {
          color: black !important;
        }

        [data-bs-theme='dark'] .nav-link.active .network-analysis {
          color: white !important;
        }

        .accordion-button {
          font-weight: bold;
        }

        /* Analysis modal with stop button */
        .analysis-modal-content {
          text-align: center;
          padding: 20px;
        }
        .analysis-spinner {
          display: inline-block;
          width: 50px;
          height: 50px;
          border: 4px solid #f3f3f3;
          border-top: 4px solid #2FA4E7;
          border-radius: 50%;
          animation: spin 1s linear infinite;
          margin-bottom: 15px;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        .analysis-text {
          font-size: 16px;
          color: #666;
          margin-bottom: 20px;
        }
        .stop-analysis-btn {
          background-color: #dc3545;
          color: white;
          border: none;
          padding: 10px 30px;
          border-radius: 5px;
          cursor: pointer;
          font-size: 14px;
          transition: background-color 0.3s;
        }
        .stop-analysis-btn:hover {
          background-color: #c82333;
        }
        .stop-analysis-btn:disabled {
          background-color: #6c757d;
          cursor: not-allowed;
        }

        /* Sample data warning - blinking */
        .sample-data-warning {
          background-color: #fff3cd;
          border: 2px solid #ffc107;
          border-radius: 8px;
          padding: 10px 15px;
          margin: 10px 0;
          animation: warning-blink 1.5s ease-in-out infinite;
        }

        @keyframes warning-blink {
          0%, 100% {
            background-color: #fff3cd;
            border-color: #ffc107;
          }
          50% {
            background-color: #ffe69c;
            border-color: #ff9800;
          }
        }

        .sample-data-warning .warning-icon {
          color: #856404;
          font-size: 18px;
          margin-right: 8px;
        }

        .sample-data-warning .warning-text {
          color: #856404;
          font-weight: bold;
          font-size: 12px;
        }
      "))
    )
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "cerulean",
    primary = "#2FA4E7",
    navbar_bg = "#2FA4E7",
    navbar_fg = "white"
  ),
  title = tags$span(
    tags$i(class = "fas fa-project-diagram", style = "margin-right: 10px;"),
    "metaFun: INTERACTIVE_NETWORK"
  ),
  # GitHub link
  nav_item(
    tags$a(
      href = "https://github.com/aababc1/metaFun",
      target = "_blank",
      class = "nav-link",
      tags$i(class = "fab fa-github", style = "font-size:20px; color:white;"),
      " GitHub"
    )
  ),
  # Documentation link
  nav_item(
    tags$a(
      href = "https://metafun-doc.readthedocs.io/en/latest/index.html",
      target = "_blank",
      class = "nav-link",
      tags$i(class = "fas fa-book", style = "font-size:20px; color:white;"),
      " Docs"
    )
  ),
  # MGSSB Database link
  nav_item(
    tags$a(
      href = "http://www.microbiome.re.kr/home_design/Database.html",
      target = "_blank",
      class = "nav-link",
      tags$i(class = "fas fa-database", style = "font-size:20px; color:white;"),
      " MGSSB"
    )
  ),
  # Dark mode toggle
  nav_item(
    class = "ms-auto",
    input_dark_mode(id = "dark_mode", mode = "light", label = "Dark Mode")
  ),

  # Master Controller sidebar
  bslib::layout_sidebar(
    sidebar = sidebar(
      id = "master-sidebar",
      title = "Master Controller",
      open = "always",  # Keep sidebar always visible
      gap = "20px",  # Add spacing between sidebar and main content

      accordion(
        open = c("Data Loading & Configuration", "Analysis Parameters"),  # Open both panels by default
        # Data Loading Panel
        accordion_panel(
          title = "Data Loading & Configuration",
          icon = bsicons::bs_icon("upload"),

          h5("Phyloseq RDS File"),
          textInput("phyloseq_rds_file", NULL,
                   value = "data/sample_phyloseq.RDS"),

          # Sample data warning (shown when using sample data)
          uiOutput("sample_data_warning"),

          h5("Group Column"),
          selectInput("group_column", NULL, choices = NULL),

          h5("Taxonomic Aggregation"),
          selectInput("tax_rank", "Aggregate at Taxonomic Rank:",
                     choices = c("Species" = "Species",
                               "Genus" = "Genus",
                               "Family" = "Family",
                               "Order" = "Order",
                               "Class" = "Class",
                               "Phylum" = "Phylum",
                               "Kingdom" = "Kingdom"),
                     selected = "Species") |>
            tooltip("Select taxonomic rank for aggregation. Aggregates abundances by the selected rank before network analysis."),

          # Info box explaining aggregation
          div(style = "background-color: #e7f3ff; padding: 10px; border-radius: 4px; margin: 5px 0; font-size: 12px; border-left: 4px solid #2196F3;",
            icon("info-circle", style = "color: #1976D2;"),
            HTML(" Aggregates abundances by taxonomic rank before correlation analysis")
          ),

          hr(),

          div(style = "background-color: #e7f3ff; padding: 15px; border-radius: 8px; margin: 10px 0;",
            actionButton("load_data",
                        "Load All Data",
                        icon = icon("upload"),
                        class = "btn-primary btn-lg",
                        width = "100%",
                        style = "font-weight: bold;") |>
              tooltip("Load taxonomic data and metadata"),

            br(), br(),

            actionButton("reset_data",
                        "Reset All Data",
                        icon = icon("redo"),
                        class = "btn-warning",
                        width = "100%") |>
              tooltip("Reset all loaded data")
          )
        ),

        # Analysis Parameters Panel
        accordion_panel(
          title = "Analysis Parameters",
          icon = bsicons::bs_icon("sliders"),

          h5("Feature Filtering"),
          sliderInput("prevalence_threshold", "Prevalence (%):",
                     min = 1, max = 50, value = 10, step = 1),
          sliderInput("ra_threshold", "Min Relative Abundance (%):",
                     min = 0.01, max = 1, value = 0.1, step = 0.01),

          selectInput("correlation_method", "Method:",
                     choices = c("FastSpar (Compositional)" = "fastspar",
                               "FlashWeave (Conditional Independence)" = "flashweave"),
                     selected = "fastspar"),

          # FlashWeave info - explain that it uses its own thresholds
          conditionalPanel(
            condition = "input.correlation_method == 'flashweave'",
            div(
              style = "background-color: #d1ecf1; padding: 10px; border-radius: 4px; margin: 10px 0; font-size: 12px; border-left: 4px solid #17a2b8;",
              icon("info-circle", style = "color: #0c5460;"),
              strong(" FlashWeave Thresholds", style = "color: #0c5460;"),
              br(),
              "FlashWeave uses its own statistical testing (alpha parameter below).",
              br(),
              "Correlation threshold is not applicable - edges are determined by conditional independence testing."
            )
          ),

          # FastSpar-specific parameters
          conditionalPanel(
            condition = "input.correlation_method == 'fastspar'",
            h5("FastSpar Parameters", style = "margin-top: 15px; color: #2FA4E7;"),

            sliderInput("fastspar_iterations", "Iterations:",
                       min = 10, max = 100, value = 20, step = 10) |>
              tooltip("Number of iterations for median calculation (default: 20)"),
            sliderInput("fastspar_bootstraps", "Bootstraps:",
                       min = 50, max = 2000, value = 100, step = 50) |>
              tooltip("Number of bootstrap replicates for p-values (default: 100, max: 2000)"),
            numericInput("fastspar_threads", "CPU Threads:",
                        value = min(30, parallel::detectCores(logical = FALSE)),
                        min = 1, max = min(64, parallel::detectCores()), step = 1) |>
              tooltip(paste0("Number of CPU threads (detected: ", parallel::detectCores(), " cores, max: 64)"))
          ),

          # FlashWeave-specific parameters
          conditionalPanel(
            condition = "input.correlation_method == 'flashweave'",
            h5("FlashWeave Parameters", style = "margin-top: 15px; color: #17a2b8;"),

            # Mode selection
            selectInput("flashweave_mode", "Algorithm Mode:",
                       choices = c("Fast (F): Univariate, all associations" = "FALSE",
                                 "Sensitive (S): Multivariate, direct only" = "TRUE"),
                       selected = "FALSE"),

            # Max_k slider
            sliderInput("flashweave_max_k", "max_k (confounding variables):",
                       min = 0, max = 5, value = 0, step = 1),

            # Alpha
            selectInput("flashweave_alpha", "alpha (threshold):",
                       choices = c("0.001" = "0.001",
                                 "0.01" = "0.01",
                                 "0.05" = "0.05",
                                 "0.1" = "0.1"),
                       selected = "0.05"),

            # Concise info box
            div(style = "background-color: #e7f3ff; padding: 10px; border-radius: 4px; margin: 5px 0; font-size: 11px; border-left: 4px solid #17a2b8;",
              HTML("<strong>Modes:</strong> Fast (F): Univariate; rapid, all associations. Sensitive (S): Multivariate; ~50% slower, direct only<br>
                   <strong>max_k:</strong> Max confounders (0=univariate/Fast, 3=Sensitive default); higher = sparser<br>
                   <strong>Alpha:</strong> Edge threshold (default: 0.05); lower = sparser, higher confidence<br>
                   <strong>HE:</strong> Multi-habitat data. <strong>FDR:</strong> Applied by default")
            ),

            sliderInput("flashweave_threads", "Julia Threads:",
                       min = 1, max = 16, value = 4, step = 1) |>
              tooltip("Number of Julia threads to use (default: 4)"),

            checkboxInput("flashweave_heterogeneous", "Heterogeneous Mode",
                         value = FALSE) |>
              tooltip("Enable for mixed data types. Usually FALSE for microbiome data"),

            checkboxInput("flashweave_fdr", "FDR Correction",
                         value = FALSE) |>
              tooltip("Apply False Discovery Rate correction. Reduces edges. Usually FALSE for initial analysis")
          )
        )
      ),

      hr(),

      # Run button always visible
      div(style = "background-color: #d4edda; padding: 15px; border-radius: 8px; margin: 10px 0;",
        # Instructional message (will disappear after first click)
        div(
          id = "run_instruction",
          style = "color: red; font-weight: bold; margin-bottom: 10px; text-align: center;",
          "Adjust parameters and metadata of your interest to run network analysis."
        ),

        actionButton(
          inputId = "master_run_analysis",
          label = "Run Network Analysis",
          icon = icon("play"),
          class = "btn-primary btn-lg blink",
          width = "100%",
          style = "font-weight: bold;"
        ) |>
        tooltip("Run network analysis with current settings")
      ),

      # Interactive Threshold Panel - appears after FastSpar analysis completes
      uiOutput("interactive_threshold_panel")
    ),

    # Main navigation panels
    div(
      style = "padding-left: 15px;",  # Add left padding to main content
      navset_underline(
        nav_panel(
        tags$span("Data Overview", class = "network-analysis"),
        data_loading_UI("data_overview")
      ),
      nav_panel(
        tags$span("Network Analysis", class = "network-analysis"),
        network_analysis_UI("network_analysis")
      ),
      nav_panel(
        tags$span("Influential Nodes", class = "network-analysis"),
        influential_nodes_UI("influential_nodes")
      ),
      nav_panel(
        tags$span("Network Robustness", class = "network-analysis"),
        network_robustness_UI("network_robustness")
      ),
      nav_panel(
        tags$span("Network Visualization", class = "network-analysis"),
        network_visualization_UI("network_viz")
      )
      )  # Close navset_underline
    )  # Close div with padding
  )  # Close layout_sidebar
)  # Close page_navbar

# =============================================================================
# Server Definition
# =============================================================================

server <- function(input, output, session) {

  # Reactive values for data management
  rv <- reactiveValues(
    raw_data = NULL,
    metadata = NULL,
    phyloseq_obj = NULL,
    filtered_ps = NULL,
    network_results = NULL,
    raw_network_results = NULL,   # Store original results before interactive filtering
    analysis_method = NULL,       # Track which method was used (fastspar, flashweave)
    current_r_threshold = NULL,   # Current interactive r threshold
    current_p_threshold = NULL,   # Current interactive p threshold
    data_loaded = FALSE,
    analysis_completed = FALSE,
    analysis_running = FALSE,
    stop_analysis_requested = FALSE,
    analysis_process = NULL,      # callr process handle for killing
    analysis_temp_dir = NULL,     # temp directory for cleanup
    reset_stability_test = NULL   # Trigger to reset stability test results
  )

  # Initial UI state configuration
  shinyjs::disable("master_run_analysis")
  shinyjs::disable("reset_data")

  # Sample data warning - shows blinking warning when sample data is used
  output$sample_data_warning <- renderUI({
    # Check if the file path contains "sample_phyloseq.RDS"
    if (!is.null(input$phyloseq_rds_file) && grepl("sample_phyloseq\\.RDS$", input$phyloseq_rds_file)) {
      div(
        class = "sample-data-warning",
        icon("exclamation-triangle", class = "warning-icon"),
        span(class = "warning-text", "SAMPLE DATA - For testing only. Please use your own phyloseq RDS file for real analysis.")
      )
    } else {
      NULL
    }
  })

  # Load metadata and update group column choices
  observeEvent(input$metadata_file, {
    if (file.exists(input$metadata_file)) {
      tryCatch({
        # Detect delimiter
        sep_char <- ifelse(grepl("\\.csv$", input$metadata_file), ",", "\t")

        meta <- read.table(input$metadata_file, header = TRUE, sep = sep_char,
                          stringsAsFactors = FALSE, check.names = FALSE)

        # Find categorical columns (exclude potential ID columns)
        id_cols <- c("SampleID", "Sample_ID", "sample_id", "accession_used_in_analysis",
                    "sample", "Sample", "ID", "bioproject_accession")

        cat_cols <- names(meta)[sapply(names(meta), function(col_name) {
          x <- meta[[col_name]]
          # Exclude ID columns and include categorical columns
          !col_name %in% id_cols &&
            (is.character(x) || is.factor(x) || (is.numeric(x) && length(unique(x)) < 10))
        })]

        message("Available group columns: ", paste(cat_cols, collapse = ", "))

        # Auto-detect common group column names (case-insensitive)
        # Priority order: group-related names first
        common_group_names <- c(
          "group", "Group", "GROUP",
          "condition", "Condition", "CONDITION",
          "treatment", "Treatment", "TREATMENT",
          "status", "Status", "STATUS",
          "disease", "Disease", "DISEASE",
          "disease_group", "Disease_Group",
          "study_condition", "Study_Condition",
          "host_disease", "Host_Disease",
          "case_control", "Case_Control",
          "phenotype", "Phenotype"
        )

        # Find first matching column
        default_col <- NULL
        for (col_name in common_group_names) {
          if (col_name %in% cat_cols) {
            default_col <- col_name
            break
          }
        }

        # If no common name found, use first categorical column
        if (is.null(default_col) && length(cat_cols) > 0) {
          default_col <- cat_cols[1]
        }

        updateSelectInput(session, "group_column",
                         choices = cat_cols,
                         selected = default_col)

        message("Selected group column: ", default_col)
      }, error = function(e) {
        showNotification(paste("Error reading metadata:", e$message), type = "error")
      })
    }
  })

  # Data loading event
  observeEvent(input$load_data, {

    tryCatch({
      show_modal_spinner(spin = "fading-circle", text = "Loading data...")

      # Load data from phyloseq RDS (simplified - only RDS supported now)
      {
        # Load pre-built phyloseq object from RDS
        req(input$phyloseq_rds_file)
        ps <- readRDS(input$phyloseq_rds_file)

        # Extract metadata from phyloseq sample_data
        metadata <- as(sample_data(ps), "data.frame")

        # Extract raw data from phyloseq for reference
        raw_data <- list(
          abundance = as.matrix(otu_table(ps)),
          taxonomy = as.matrix(tax_table(ps))
        )

        # Update group column dropdown with columns from phyloseq
        id_cols <- c("SampleID", "Sample_ID", "sample_id", "accession_used_in_analysis",
                    "sample", "Sample", "ID", "bioproject_accession")

        cat_cols <- names(metadata)[sapply(names(metadata), function(col_name) {
          x <- metadata[[col_name]]
          !col_name %in% id_cols &&
            (is.character(x) || is.factor(x) || (is.numeric(x) && length(unique(x)) < 10))
        })]

        # Auto-detect common group column names (case-insensitive)
        common_group_names <- c(
          "group", "Group", "GROUP",
          "condition", "Condition", "CONDITION",
          "treatment", "Treatment", "TREATMENT",
          "status", "Status", "STATUS",
          "disease", "Disease", "DISEASE",
          "disease_group", "Disease_Group",
          "study_condition", "Study_Condition",
          "host_disease", "Host_Disease",
          "case_control", "Case_Control",
          "phenotype", "Phenotype"
        )

        # Find first matching column
        default_col <- NULL
        for (col_name in common_group_names) {
          if (col_name %in% cat_cols) {
            default_col <- col_name
            break
          }
        }

        # If no common name found, use first categorical column
        if (is.null(default_col) && length(cat_cols) > 0) {
          default_col <- cat_cols[1]
        }

        updateSelectInput(session, "group_column",
                         choices = cat_cols,
                         selected = default_col)
      }

      # Store in reactive values
      rv$raw_data <- raw_data
      rv$metadata <- metadata
      rv$phyloseq_obj <- ps
      rv$data_loaded <- TRUE
      rv$analysis_completed <- FALSE

      # Enable controls
      shinyjs::enable("master_run_analysis")
      shinyjs::enable("reset_data")
      shinyjs::removeClass(selector = "#load_data", class = "blink")

      remove_modal_spinner()

      # Debug: Show current tax_rank selection
      message("=== DATA LOADED ===")
      message("Current tax_rank selection: ", input$tax_rank)
      message("===================")

      # Show notification based on whether group column is set
      if (is.null(input$group_column) || input$group_column == "") {
        showNotification(
          "Data loaded! Please select a Group Column from the dropdown above to enable network analysis.",
          type = "warning",
          duration = 10
        )
      } else {
        showNotification("Data loaded successfully!", type = "message")
      }

    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })

  # Reset data event
  observeEvent(input$reset_data, {
    rv$raw_data <- NULL
    rv$metadata <- NULL
    rv$phyloseq_obj <- NULL
    rv$filtered_ps <- NULL
    rv$network_results <- NULL
    rv$raw_network_results <- NULL
    rv$analysis_method <- NULL
    rv$current_r_threshold <- NULL
    rv$current_p_threshold <- NULL
    rv$data_loaded <- FALSE
    rv$analysis_completed <- FALSE

    shinyjs::disable("master_run_analysis")
    shinyjs::disable("reset_data")
    shinyjs::addClass(selector = "#load_data", class = "blink")

    showNotification("All data has been reset.", type = "message")
  })

  # Helper function to show analysis modal with stop button
  show_analysis_modal <- function() {
    showModal(modalDialog(
      div(
        class = "analysis-modal-content",
        div(class = "analysis-spinner"),
        div(class = "analysis-text", "Running network analysis..."),
        actionButton(
          "stop_analysis_btn",
          "Stop Analysis",
          class = "stop-analysis-btn",
          icon = icon("stop")
        )
      ),
      title = NULL,
      footer = NULL,
      size = "s",
      easyClose = FALSE
    ))
  }

  # Helper function to kill child processes (fastspar, julia, etc.)
  kill_child_processes <- function(parent_pid) {
    tryCatch({
      # Find all child processes of the R background process
      child_pids <- system(sprintf("pgrep -P %s", parent_pid), intern = TRUE)

      # Also find fastspar and julia processes that might be running
      fastspar_pids <- system("pgrep -f fastspar", intern = TRUE)
      julia_pids <- system("pgrep -f 'julia.*flashweave'", intern = TRUE)

      all_pids <- unique(c(child_pids, fastspar_pids, julia_pids))

      for (pid in all_pids) {
        if (nchar(pid) > 0) {
          system(sprintf("kill -9 %s 2>/dev/null", pid), ignore.stdout = TRUE, ignore.stderr = TRUE)
          message("  Killed process: ", pid)
        }
      }
    }, error = function(e) {
      message("  Warning: Could not kill some child processes: ", e$message)
    })
  }

  # Helper function to clean up temp files
  cleanup_temp_files <- function() {
    tryCatch({
      tmp_dir <- tempdir()

      # Clean up FastSpar files
      fastspar_files <- list.files(tmp_dir, pattern = "^fastspar", full.names = TRUE)
      for (f in fastspar_files) {
        if (file.exists(f)) {
          if (dir.exists(f)) {
            unlink(f, recursive = TRUE)
          } else {
            file.remove(f)
          }
        }
      }

      # Clean up FlashWeave files
      flashweave_files <- list.files(tmp_dir, pattern = "^flashweave", full.names = TRUE)
      for (f in flashweave_files) {
        if (file.exists(f)) {
          if (dir.exists(f)) {
            unlink(f, recursive = TRUE)
          } else {
            file.remove(f)
          }
        }
      }

      # Clean up bootstrap directories
      boot_dirs <- list.files(tmp_dir, pattern = "bootstrap", full.names = TRUE)
      for (d in boot_dirs) {
        if (dir.exists(d)) {
          unlink(d, recursive = TRUE)
        }
      }

      message("  Cleaned up temporary files")
    }, error = function(e) {
      message("  Warning: Could not clean some temp files: ", e$message)
    })
  }

  # Stop analysis button observer - KILLS PROCESS IMMEDIATELY
  observeEvent(input$stop_analysis_btn, {
    if (rv$analysis_running && !is.null(rv$analysis_process)) {
      message("\n=== STOPPING ANALYSIS ===")

      # Update modal to show stopping state
      showModal(modalDialog(
        div(
          class = "analysis-modal-content",
          div(class = "analysis-spinner"),
          div(class = "analysis-text", "Stopping analysis..."),
          actionButton(
            "stop_analysis_btn",
            "Stopping...",
            class = "stop-analysis-btn",
            icon = icon("spinner", class = "fa-spin"),
            disabled = TRUE
          )
        ),
        title = NULL,
        footer = NULL,
        size = "s",
        easyClose = FALSE
      ))

      tryCatch({
        # Get the process ID before killing
        if (!is.null(rv$analysis_process) && rv$analysis_process$is_alive()) {
          proc_pid <- rv$analysis_process$get_pid()
          message("  Killing main process PID: ", proc_pid)

          # Kill child processes first (fastspar, julia, etc.)
          kill_child_processes(proc_pid)

          # Kill the main R background process
          rv$analysis_process$kill()
          message("  Main process killed")
        }

        # Clean up temp files
        cleanup_temp_files()

        # Reset state
        rv$analysis_running <- FALSE
        rv$stop_analysis_requested <- FALSE
        rv$analysis_process <- NULL

        removeModal()
        showNotification("Analysis stopped and cleaned up.", type = "warning", duration = 5)
        message("=== ANALYSIS STOPPED SUCCESSFULLY ===\n")

      }, error = function(e) {
        message("  Error during stop: ", e$message)
        rv$analysis_running <- FALSE
        rv$analysis_process <- NULL
        removeModal()
        showNotification(
          paste("Error stopping analysis:", e$message),
          type = "error",
          duration = 10
        )
      })
    }
  })

  # Observer to check background process status
  observe({
    invalidateLater(500)  # Check every 500ms

    if (rv$analysis_running && !is.null(rv$analysis_process)) {
      # Check if process has finished
      if (!rv$analysis_process$is_alive()) {
        tryCatch({
          # Get the result
          result <- rv$analysis_process$get_result()

          if (!is.null(result) && !is.null(result$network_res)) {
            # Validate results
            if (!"comparison" %in% names(result$network_res)) {
              stop("Network results missing 'comparison' dataframe")
            }

            comp_df <- result$network_res$comparison

            # Store results
            rv$filtered_ps <- result$filtered_ps
            rv$network_results <- result$network_res
            rv$raw_network_results <- result$network_res  # Store original for interactive filtering
            rv$analysis_method <- result$method
            rv$current_r_threshold <- result$initial_r_threshold
            rv$current_p_threshold <- result$initial_p_threshold
            rv$analysis_completed <- TRUE
            rv$tax_rank <- result$tax_rank

            # Log success
            message("\n=== NETWORK ANALYSIS COMPLETED ===")
            message("Analysis completed: TRUE")
            message("Results components: ", paste(names(result$network_res), collapse = ", "))
            if (!is.null(comp_df) && nrow(comp_df) > 0) {
              message("Groups analyzed: ", paste(comp_df$group, collapse = ", "))
              message("Total edges: ", paste(comp_df$n_edges, collapse = ", "))
            }
            message("==================================\n")

            removeModal()

            # Provide feedback about edges
            total_edges <- sum(comp_df$n_edges, na.rm = TRUE)
            if (total_edges == 0) {
              showNotification(
                HTML(paste0(
                  "<strong>Network analysis completed</strong><br>",
                  "Warning: No edges detected!<br>",
                  "Try lowering r.threshold or increasing p.threshold"
                )),
                type = "warning",
                duration = 15
              )
            } else {
              showNotification(
                paste0("Network analysis completed! Total edges: ", total_edges),
                type = "message",
                duration = 5
              )
            }
          } else {
            stop("Network analysis returned NULL results")
          }

        }, error = function(e) {
          message("=== NETWORK ANALYSIS ERROR ===")
          message("Error message: ", e$message)

          # Try to get more detailed error from subprocess
          tryCatch({
            if (!is.null(rv$analysis_process)) {
              # Read stderr/stdout for detailed error
              stderr_output <- rv$analysis_process$read_error()
              stdout_output <- rv$analysis_process$read_output()

              if (length(stderr_output) > 0 && nchar(paste(stderr_output, collapse = "")) > 0) {
                message("STDERR from subprocess:")
                message(paste(stderr_output, collapse = "\n"))
              }
              if (length(stdout_output) > 0 && nchar(paste(stdout_output, collapse = "")) > 0) {
                message("STDOUT from subprocess:")
                message(paste(stdout_output, collapse = "\n"))
              }
            }
          }, error = function(e2) {
            message("Could not read subprocess output: ", e2$message)
          })

          message("==============================")

          rv$network_results <- NULL
          rv$analysis_completed <- FALSE

          removeModal()
          showNotification(
            HTML(paste0(
              "<strong>Network Analysis Error:</strong><br>",
              e$message, "<br><br>",
              "<em>Check R console for details</em>"
            )),
            type = "error",
            duration = NULL
          )
        })

        # Reset running state
        rv$analysis_running <- FALSE
        rv$analysis_process <- NULL
      }
    }
  })

  # Master analysis event
  observeEvent(input$master_run_analysis, {
    req(rv$data_loaded, rv$phyloseq_obj)

    # Prevent starting new analysis if one is already running
    if (rv$analysis_running) {
      showNotification("Analysis is already running!", type = "warning")
      return(NULL)
    }

    # Validate group column selection
    if (is.null(input$group_column) || input$group_column == "") {
      showNotification(
        "Please select a Group Column from the Master Controller sidebar!",
        type = "error",
        duration = 10
      )
      return(NULL)
    }

    # Hide instruction message and remove blink after first click
    shinyjs::hide("run_instruction")
    shinyjs::removeClass(selector = "#master_run_analysis", class = "blink")

    # Reset state and set running
    rv$stop_analysis_requested <- FALSE
    rv$analysis_running <- TRUE
    rv$analysis_process <- NULL

    # Show custom modal with stop button
    show_analysis_modal()

    # Log analysis parameters
    message("\n=== STARTING NETWORK ANALYSIS ===")
    message("Method: ", input$correlation_method)
    message("Group column: ", input$group_column)
    message("Taxonomic rank: ", input$tax_rank)
    message("Prevalence: ", input$prevalence_threshold, "%")
    message("Min RA: ", input$ra_threshold, "%")
    message("=================================\n")

    # Capture all inputs for background processing
    ps_obj <- rv$phyloseq_obj
    tax_rank_val <- input$tax_rank
    group_col <- input$group_column
    prev_thresh <- input$prevalence_threshold / 100
    ra_thresh <- input$ra_threshold / 100
    corr_method <- input$correlation_method
    # For FastSpar, use default thresholds since post-hoc filtering is used
    # For FlashWeave, no thresholds are needed (all edges are significant)
    if (corr_method == "fastspar") {
      r_thresh <- 0.2  # Default for post-hoc filtering
      p_thresh <- 0.05  # Default for post-hoc filtering
    } else {
      r_thresh <- 0  # FlashWeave uses all inferred edges
      p_thresh <- 1  # FlashWeave uses all inferred edges
    }
    current_app_dir <- app_dir

    # Prepare method-specific parameters
    method_params <- list()
    if (corr_method == "fastspar") {
      method_params <- list(
        iterations = if (!is.null(input$fastspar_iterations)) input$fastspar_iterations else 50,
        bootstraps = if (!is.null(input$fastspar_bootstraps)) input$fastspar_bootstraps else 100,
        threads = if (!is.null(input$fastspar_threads)) input$fastspar_threads else 4
      )
      message("[APP] FastSpar parameters collected: ", paste(names(method_params), "=", method_params, collapse=", "))
    } else if (corr_method == "flashweave") {
      method_params <- list(
        sensitive = if (!is.null(input$flashweave_mode)) as.logical(input$flashweave_mode) else FALSE,
        heterogeneous = if (!is.null(input$flashweave_heterogeneous)) input$flashweave_heterogeneous else FALSE,
        max_k = if (!is.null(input$flashweave_max_k)) input$flashweave_max_k else 0,
        alpha = if (!is.null(input$flashweave_alpha)) as.numeric(input$flashweave_alpha) else 0.1,
        n_threads = if (!is.null(input$flashweave_threads)) input$flashweave_threads else 4,
        FDR = if (!is.null(input$flashweave_fdr)) input$flashweave_fdr else FALSE
      )
      message("[APP] FlashWeave parameters collected: ", paste(names(method_params), "=", method_params, collapse=", "))
    } else {
      message("[APP] Using method '", corr_method, "' with no custom parameters")
    }

    # Run analysis in a KILLABLE background process using callr
    # IMPORTANT: Pass environment variables for Julia/FastSpar to work in subprocess
    rv$analysis_process <- callr::r_bg(
      func = function(ps_obj, tax_rank_val, group_col, prev_thresh, ra_thresh,
                      r_thresh, p_thresh, corr_method, method_params, app_dir) {
        # Source helper functions in the background process
        source(file.path(app_dir, "helper/data_conversion.R"))
        source(file.path(app_dir, "helper/network_utils.R"))

        # Load required packages
        library(phyloseq)
        library(igraph)
        library(dplyr)

        # Step 1: Aggregate by taxonomic rank (if selected)
        ps_working <- aggregate_taxa(ps_obj, tax_rank_val)

        # Step 2: Apply filtering
        ps_filtered <- apply_filtering(
          ps_working,
          strategy = "prev_ra",
          top_n = NULL,
          prevalence = prev_thresh,
          ra_threshold = ra_thresh,
          group = group_col
        )

        # Run network analysis (group-wise)
        network_res <- run_network_analysis(
          ps_filtered,
          group = group_col,
          r.threshold = r_thresh,
          p.threshold = p_thresh,
          method = corr_method,
          method_params = method_params
        )

        list(
          filtered_ps = ps_filtered,
          network_res = network_res,
          tax_rank = tax_rank_val,
          method = corr_method,
          initial_r_threshold = r_thresh,
          initial_p_threshold = p_thresh
        )
      },
      args = list(
        ps_obj = ps_obj,
        tax_rank_val = tax_rank_val,
        group_col = group_col,
        prev_thresh = prev_thresh,
        ra_thresh = ra_thresh,
        r_thresh = r_thresh,
        p_thresh = p_thresh,
        corr_method = corr_method,
        method_params = method_params,
        app_dir = current_app_dir
      ),
      env = c(
        callr::rcmd_safe_env(),
        PATH = Sys.getenv("PATH"),
        LD_LIBRARY_PATH = Sys.getenv("LD_LIBRARY_PATH"),
        JULIA_DEPOT_PATH = Sys.getenv("JULIA_DEPOT_PATH", "/opt/julia/.julia")
      ),
      stdout = "|",  # Capture stdout
      stderr = "|",  # Capture stderr
      supervise = TRUE  # Allows killing child processes
    )

    message("Background process started with PID: ", rv$analysis_process$get_pid())

    # Return NULL to not block the UI
    NULL
  })

  # Analysis status output
  output$analysis_status <- renderUI({
    if (!rv$data_loaded) {
      div(
        style = "color: red; font-weight: bold;",
        icon("exclamation-triangle"),
        " Load data first"
      )
    } else if (!rv$analysis_completed) {
      div(
        style = "color: orange; font-weight: bold;",
        icon("info-circle"),
        " Ready to analyze"
      )
    } else {
      div(
        style = "color: #2FA4E7; font-weight: bold;",
        icon("check-circle"),
        " Analysis completed"
      )
    }
  })

  # ==========================================================================
  # Interactive Threshold Panel (appears after FastSpar analysis completes)
  # ==========================================================================
  output$interactive_threshold_panel <- renderUI({
    # Only show if analysis completed and method is FastSpar
    req(rv$analysis_completed, rv$analysis_method == "fastspar")

    div(
      style = "background-color: #e3f2fd; padding: 15px; border-radius: 8px; margin-top: 15px; border-left: 4px solid #2196f3;",
      h5(icon("sliders-h"), " Interactive Threshold", style = "color: #1565c0; margin-bottom: 15px;"),

      div(
        style = "background-color: #bbdefb; padding: 10px; border-radius: 4px; margin-bottom: 10px; font-size: 11px;",
        icon("info-circle", style = "color: #1565c0;"),
        " Adjust thresholds without re-running FastSpar analysis. Network will update instantly."
      ),

      sliderInput("interactive_r_threshold", "r (Correlation):",
                  min = 0.05, max = 0.9,
                  value = if (!is.null(isolate(rv$current_r_threshold))) isolate(rv$current_r_threshold) else 0.2,
                  step = 0.01),

      sliderInput("interactive_p_threshold", "p-value:",
                  min = 0.001, max = 0.1,
                  value = if (!is.null(isolate(rv$current_p_threshold))) isolate(rv$current_p_threshold) else 0.05,
                  step = 0.001),

      # Edge count preview
      uiOutput("edge_count_preview"),

      actionButton(
        inputId = "apply_interactive_threshold",
        label = "Apply Thresholds",
        icon = icon("sync"),
        class = "btn-info",
        width = "100%",
        style = "margin-top: 10px;"
      )
    )
  })

  # Node count preview (updates as sliders change)
  output$edge_count_preview <- renderUI({
    req(rv$raw_network_results, input$interactive_r_threshold, input$interactive_p_threshold)

    # Count nodes (non-isolated vertices) with current slider values
    r_thresh <- input$interactive_r_threshold
    p_thresh <- input$interactive_p_threshold

    node_counts <- sapply(names(rv$raw_network_results), function(g) {
      if (is.null(rv$raw_network_results[[g]])) return(0)
      if (g == "comparison") return(NA)  # Skip comparison dataframe

      cor_mat <- rv$raw_network_results[[g]]$correlation_matrix
      p_mat <- rv$raw_network_results[[g]]$p_matrix

      if (is.null(cor_mat) || is.null(p_mat)) return(0)

      # Count edges that pass both thresholds
      pass_filter <- (abs(cor_mat) >= r_thresh) & (p_mat <= p_thresh)
      diag(pass_filter) <- FALSE  # Exclude self-loops

      # Count nodes that have at least one edge (non-isolated)
      nodes_with_edges <- rowSums(pass_filter) > 0 | colSums(pass_filter) > 0
      sum(nodes_with_edges)
    })

    # Remove NA values (comparison)
    node_counts <- node_counts[!is.na(node_counts)]

    div(
      style = "background-color: #fff3e0; padding: 10px; border-radius: 4px; margin-top: 10px; font-size: 12px;",
      strong("Node Preview:"),
      br(),
      lapply(names(node_counts), function(g) {
        div(sprintf("%s: %.0f nodes", g, node_counts[[g]]))
      })
    )
  })

  # Apply interactive thresholds
  observeEvent(input$apply_interactive_threshold, {
    req(rv$raw_network_results, rv$analysis_method == "fastspar")

    r_thresh <- input$interactive_r_threshold
    p_thresh <- input$interactive_p_threshold

    # Show loading modal
    showModal(modalDialog(
      title = "Applying Thresholds",
      div(
        style = "text-align: center; padding: 20px;",
        tags$i(class = "fa fa-spinner fa-spin fa-3x", style = "color: #007bff;"),
        br(), br(),
        h4("Recalculating network with new thresholds..."),
        p(style = "color: #666;",
          sprintf("r >= %.2f, p <= %.3f", r_thresh, p_thresh)),
        p(style = "color: #888; font-size: 12px;",
          "Computing centrality metrics, Zi-Pi values, and IVI scores...")
      ),
      footer = NULL,
      easyClose = FALSE
    ))

    # Trigger stability test reset (observed by robustness module)
    rv$reset_stability_test <- Sys.time()

    message("\n=== APPLYING INTERACTIVE THRESHOLDS ===")
    message("r.threshold: ", r_thresh)
    message("p.threshold: ", p_thresh)

    # Aggressive cleanup of stale connections to prevent "all connections in use" error
    tryCatch({
      all_cons <- showConnections(all = TRUE)
      if (nrow(all_cons) > 3) {  # Keep stdin, stdout, stderr
        for (i in seq_len(nrow(all_cons))) {
          con_id <- as.integer(rownames(all_cons)[i])
          if (con_id >= 3) {  # Don't close stdin(0), stdout(1), stderr(2)
            tryCatch(close(getConnection(con_id)), error = function(e) NULL)
          }
        }
      }
    }, error = function(e) NULL)
    gc(verbose = FALSE)

    # Re-filter network for each group
    updated_results <- rv$raw_network_results

    for (g in names(updated_results)) {
      if (g == "comparison") next
      if (is.null(updated_results[[g]])) next

      tryCatch({
        # Get raw matrices
        cor_matrix <- updated_results[[g]]$correlation_matrix
        p_matrix <- updated_results[[g]]$p_matrix
        ps_group <- updated_results[[g]]$phyloseq

        # Create adjacency matrix with new thresholds
        n_taxa <- nrow(cor_matrix)
        adj_matrix <- matrix(1, nrow = n_taxa, ncol = n_taxa)
        rownames(adj_matrix) <- colnames(adj_matrix) <- rownames(cor_matrix)
        diag(adj_matrix) <- 0

        # Apply thresholds
        adj_matrix[abs(cor_matrix) < r_thresh | p_matrix > p_thresh] <- 0

        # Create weighted adjacency with ABSOLUTE values for igraph
        # (igraph requires positive weights for betweenness)
        weighted_adj_abs <- abs(cor_matrix) * adj_matrix

        # Build new igraph with absolute weights
        g_igraph <- igraph::graph_from_adjacency_matrix(
          weighted_adj_abs,
          mode = "max",  # Use max to handle any asymmetry
          weighted = TRUE,
          diag = FALSE
        )

        # Remove isolated vertices
        isolated <- which(igraph::degree(g_igraph) == 0)
        if (length(isolated) > 0) {
          g_igraph <- igraph::delete_vertices(g_igraph, isolated)
        }

        # Add sign attribute to edges based on original correlation matrix
        if (igraph::ecount(g_igraph) > 0) {
          # Get edge endpoints and determine sign from original correlation
          edge_list <- igraph::as_edgelist(g_igraph, names = TRUE)
          signs <- sapply(1:nrow(edge_list), function(i) {
            from_name <- edge_list[i, 1]
            to_name <- edge_list[i, 2]
            if (from_name %in% rownames(cor_matrix) && to_name %in% colnames(cor_matrix)) {
              cor_val <- cor_matrix[from_name, to_name]
              return(ifelse(cor_val > 0, "positive", "negative"))
            }
            return("positive")  # default
          })
          igraph::E(g_igraph)$sign <- signs

          # Also add signed_weight for visualization (positive/negative correlations)
          signed_weights <- sapply(1:nrow(edge_list), function(i) {
            from_name <- edge_list[i, 1]
            to_name <- edge_list[i, 2]
            if (from_name %in% rownames(cor_matrix) && to_name %in% colnames(cor_matrix)) {
              return(cor_matrix[from_name, to_name])
            }
            return(0)
          })
          igraph::E(g_igraph)$signed_weight <- signed_weights
        }

        # Transfer taxonomy attribute from original graph
        original_g <- updated_results[[g]]$igraph
        if (!is.null(original_g) && !is.null(igraph::V(original_g)$taxonomy)) {
          # Match vertices by name
          new_vertex_names <- igraph::V(g_igraph)$name
          old_vertex_names <- igraph::V(original_g)$name
          old_taxonomy <- igraph::V(original_g)$taxonomy

          matched_taxonomy <- sapply(new_vertex_names, function(v) {
            idx <- match(v, old_vertex_names)
            if (!is.na(idx)) {
              return(old_taxonomy[idx])
            }
            return(NA)
          })
          igraph::V(g_igraph)$taxonomy <- matched_taxonomy
        }

        # Recalculate metrics
        n_nodes <- igraph::vcount(g_igraph)
        n_edges <- igraph::ecount(g_igraph)

        # Calculate modularity (needed for metrics comparison plot)
        mod_value <- NA
        if (n_nodes > 0 && n_edges > 0) {
          tryCatch({
            comm <- igraph::cluster_fast_greedy(g_igraph)
            mod_value <- igraph::modularity(comm)
          }, error = function(e) {})
        }

        # Use consistent column names with original network analysis
        metrics <- data.frame(
          n_nodes = n_nodes,
          n_edges = n_edges,
          density = if (n_nodes > 1) igraph::edge_density(g_igraph) else 0,
          avg_degree = if (n_nodes > 0) mean(igraph::degree(g_igraph)) else 0,
          transitivity = if (n_nodes > 2) igraph::transitivity(g_igraph, type = "global") else 0,
          avg_path_length = NA,
          diameter = NA,
          modularity = mod_value,
          stringsAsFactors = FALSE
        )

        # Calculate basic node metrics without parallel processing
        # (skip get_comprehensive_network_metrics to avoid connection exhaustion)
        node_metrics <- NULL
        zi_pi <- NULL
        ivi <- NULL
        communities <- NULL

        if (n_nodes > 0 && n_edges > 0) {
          # Basic centrality metrics
          degree_vals <- igraph::degree(g_igraph)
          betweenness_vals <- igraph::betweenness(g_igraph, normalized = TRUE)

          # Closeness centrality (only for connected graphs)
          closeness_vals <- tryCatch({
            if (igraph::is_connected(g_igraph)) {
              igraph::closeness(g_igraph, normalized = TRUE)
            } else {
              # For disconnected graphs, calculate per component
              cl <- rep(NA, igraph::vcount(g_igraph))
              comps <- igraph::components(g_igraph)
              for (comp_id in unique(comps$membership)) {
                comp_nodes <- which(comps$membership == comp_id)
                if (length(comp_nodes) > 1) {
                  subg <- igraph::induced_subgraph(g_igraph, comp_nodes)
                  cl[comp_nodes] <- igraph::closeness(subg, normalized = TRUE)
                }
              }
              cl
            }
          }, error = function(e) rep(NA, igraph::vcount(g_igraph)))

          node_metrics <- data.frame(
            node = igraph::V(g_igraph)$name,
            degree = degree_vals,
            betweenness = betweenness_vals,
            closeness = closeness_vals,
            stringsAsFactors = FALSE
          )

          # Community detection for Zi-Pi
          tryCatch({
            communities <- igraph::cluster_fast_greedy(g_igraph)
            membership <- igraph::membership(communities)

            # Calculate Zi (within-module degree z-score)
            Zi_vals <- sapply(seq_along(membership), function(i) {
              module <- membership[i]
              module_nodes <- which(membership == module)
              ki <- sum(g_igraph[i, module_nodes])
              ks <- sapply(module_nodes, function(j) sum(g_igraph[j, module_nodes]))
              ks_mean <- mean(ks)
              ks_sd <- sd(ks)
              if (is.na(ks_sd) || ks_sd == 0) return(0)
              (ki - ks_mean) / ks_sd
            })

            # Calculate Pi (participation coefficient)
            Pi_vals <- sapply(seq_along(membership), function(i) {
              ki <- degree_vals[i]
              if (ki == 0) return(0)
              kis <- sapply(unique(membership), function(m) {
                sum(g_igraph[i, which(membership == m)])
              })
              1 - sum((kis / ki)^2)
            })

            # Classify nodes into roles based on Zi and Pi
            roles <- rep("Peripheral", length(Zi_vals))
            roles[Zi_vals > 2.5] <- "Module hub"
            roles[Pi_vals > 0.62] <- "Connector"
            roles[Zi_vals > 2.5 & Pi_vals > 0.62] <- "Network hub"

            # Use uppercase Zi/Pi to match original calculate_zi_pi function
            zi_pi <- data.frame(
              node = igraph::V(g_igraph)$name,
              module = as.integer(membership),
              Zi = Zi_vals,
              Pi = Pi_vals,
              role = roles,
              stringsAsFactors = FALSE
            )
          }, error = function(e) {
            message("  Zi-Pi calculation skipped: ", e$message)
          })

          # Calculate IVI using influential package
          # Clean up connections before IVI calculation to prevent exhaustion
          tryCatch({
            # Close any open connections to prevent "all connections in use" error
            all_cons <- showConnections(all = TRUE)
            if (nrow(all_cons) > 3) {  # Keep stdin, stdout, stderr
              for (i in 4:nrow(all_cons)) {
                tryCatch(close(getConnection(as.integer(rownames(all_cons)[i]))),
                        error = function(e) NULL)
              }
            }
            gc(verbose = FALSE)
          }, error = function(e) NULL)

          tryCatch({
            if (requireNamespace("influential", quietly = TRUE)) {
              ivi_scores <- influential::ivi(
                graph = g_igraph,
                vertices = igraph::V(g_igraph),
                directed = FALSE,
                mode = "all",
                d = 3,
                scale = "range"
              )

              hubness <- influential::hubness.score(g_igraph, scale = "range")
              spreading <- influential::spreading.score(g_igraph, d = 3, scale = "range")

              ivi <- data.frame(
                node = igraph::V(g_igraph)$name,
                hubness_score = hubness,
                spreading_score = spreading,
                IVI = ivi_scores,
                stringsAsFactors = FALSE
              )
              ivi <- ivi[order(ivi$IVI, decreasing = TRUE), ]
              ivi$rank <- 1:nrow(ivi)

              # Clean up after IVI calculation
              gc(verbose = FALSE)
            }
          }, error = function(e) {
            message("  IVI calculation skipped: ", e$message)
          })
        }

        # Update results
        updated_results[[g]]$igraph <- g_igraph
        updated_results[[g]]$metrics <- metrics
        updated_results[[g]]$n_samples <- updated_results[[g]]$n_samples  # preserve original
        updated_results[[g]]$n_taxa <- n_taxa  # total taxa before filtering
        if (!is.null(node_metrics)) updated_results[[g]]$node_metrics <- node_metrics
        if (!is.null(zi_pi)) updated_results[[g]]$zi_pi <- zi_pi
        if (!is.null(ivi)) updated_results[[g]]$ivi <- ivi
        if (!is.null(communities)) updated_results[[g]]$communities <- communities

        message(sprintf("  %s: %d edges (from %d taxa)", g, n_edges, n_taxa))

      }, error = function(e) {
        message("Error updating group ", g, ": ", e$message)
      })
    }

    # Update comparison dataframe - use consistent column names with original analysis
    groups <- setdiff(names(updated_results), "comparison")
    # Filter out groups that failed to update (missing metrics)
    valid_groups <- groups[sapply(groups, function(g) {
      !is.null(updated_results[[g]]$metrics) && !is.null(updated_results[[g]]$metrics$n_nodes)
    })]

    if (length(valid_groups) > 0) {
      comparison_df <- data.frame(
        group = valid_groups,
        n_nodes = sapply(valid_groups, function(g) updated_results[[g]]$metrics$n_nodes),
        n_edges = sapply(valid_groups, function(g) updated_results[[g]]$metrics$n_edges),
        density = sapply(valid_groups, function(g) updated_results[[g]]$metrics$density),
        avg_degree = sapply(valid_groups, function(g) updated_results[[g]]$metrics$avg_degree),
        transitivity = sapply(valid_groups, function(g) updated_results[[g]]$metrics$transitivity),
        avg_path_length = sapply(valid_groups, function(g) updated_results[[g]]$metrics$avg_path_length),
        diameter = sapply(valid_groups, function(g) updated_results[[g]]$metrics$diameter),
        modularity = sapply(valid_groups, function(g) updated_results[[g]]$metrics$modularity),
        n_samples = sapply(valid_groups, function(g) updated_results[[g]]$n_samples),
        n_taxa = sapply(valid_groups, function(g) updated_results[[g]]$n_taxa),
        stringsAsFactors = FALSE
      )
    } else {
      # Fallback: empty dataframe with correct structure (matches original analysis)
      comparison_df <- data.frame(
        group = character(0), n_nodes = integer(0), n_edges = integer(0),
        density = numeric(0), avg_degree = numeric(0), transitivity = numeric(0),
        avg_path_length = numeric(0), diameter = numeric(0), modularity = numeric(0),
        n_samples = integer(0), n_taxa = integer(0), stringsAsFactors = FALSE
      )
    }
    updated_results$comparison <- comparison_df

    # Store updated results
    rv$network_results <- updated_results
    rv$current_r_threshold <- r_thresh
    rv$current_p_threshold <- p_thresh

    # Final cleanup to prevent connection warnings
    tryCatch({
      all_cons <- showConnections(all = TRUE)
      if (nrow(all_cons) > 3) {
        for (i in seq_len(nrow(all_cons))) {
          con_id <- as.integer(rownames(all_cons)[i])
          if (con_id >= 3) {
            tryCatch(close(getConnection(con_id)), error = function(e) NULL)
          }
        }
      }
    }, error = function(e) NULL)
    gc(verbose = FALSE)

    message("=== INTERACTIVE THRESHOLD UPDATE COMPLETE ===\n")

    # Close loading modal
    removeModal()

    showNotification(
      paste0("Network updated! Total edges: ", sum(comparison_df$n_edges, na.rm = TRUE)),
      type = "message",
      duration = 3
    )
  })

  # Module servers
  # Pass input$tax_rank directly to all modules so they can sync with Master Controller
  data_loading_Server("data_overview",
                     reactive(rv$phyloseq_obj),
                     reactive(rv$data_loaded),
                     reactive(input$tax_rank))

  network_analysis_Server("network_analysis",
                         reactive(rv$filtered_ps),
                         reactive(rv$network_results),
                         reactive(rv$analysis_completed),
                         reactive(input$tax_rank))

  influential_nodes_Server("influential_nodes",
                          reactive(rv$filtered_ps),
                          reactive(rv$network_results),
                          reactive(rv$analysis_completed),
                          reactive(input$tax_rank))

  network_robustness_Server("network_robustness",
                           reactive(rv$filtered_ps),
                           reactive(rv$network_results),
                           reactive(rv$analysis_completed),
                           reactive(rv$reset_stability_test))

  network_visualization_Server("network_viz",
                              reactive(rv$network_results),
                              reactive(input$group_column),
                              reactive(rv$analysis_completed),
                              reactive(input$tax_rank))
}

# Launch app
options(shiny.autoreload = TRUE)
shinyApp(ui = ui, server = server)
