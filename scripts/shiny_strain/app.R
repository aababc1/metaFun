# =============================================================================
# INTERACTIVE_STRAIN - InStrain Analysis Visualization App
# Part of metaFun pipeline (https://github.com/aababc1/metafun)
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
library(ggVennDiagram)
library(UpSetR)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(shinythemes)
library(Cairo)
library(cowplot)
library(grid)
library(bsicons)
library(future)
library(later)
library(promises)
library(shinydashboard)

# Future configuration for async processing
future::plan(multisession)

# Source modules
app_dir <- getwd()  # Current working directory
#source(file.path(app_dir, "modules/data_processing.R")) # will be removed 
source(file.path(app_dir, "modules/shared_taxa_module.R"))
source(file.path(app_dir, "modules/nucleotide_diversity_module.R"))
# source(file.path(app_dir, "modules/prevalence_module.R"))  # Moved to old_code
# source(file.path(app_dir, "modules/diversity_module.R"))   # Moved to old_code (not in use)
source(file.path(app_dir, "modules/pnps_module_quick.R"))
source(file.path(app_dir, "modules/gene_level_pnps_module.R"))
source(file.path(app_dir, "modules/diversity_pnps_correlation_module.R"))

# Always source preprocessing functions (lightweight, only defines functions)
source(file.path(app_dir, "modules/preprocessing_scripts.R"))

# Check if preprocessed data exists
preprocessed_rds <- file.path(app_dir, "data/integrated_microbiome_data.rds")
if (file.exists(preprocessed_rds)) {
  message("Using existing preprocessed data: ", preprocessed_rds)
}
# Helper functions
source(file.path(app_dir, "helper/plot_customization.R"))
source(file.path(app_dir, "helper/error_handlers.R"))

# WWW path configuration
www_path <- file.path(app_dir, "www")
addResourcePath('www', www_path)

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

        .shiny-notification {
          position: fixed;
          top: 100px;
          left: 100px;
        }

        #master_run_analysis {
          background-color: #2FA4E7;
          color: white;
          border-color: #2FA4E7;
        }

        .strain-analysis {
          color: #2FA4E7 !important;
        }
        .diversity-analysis {
          color: #1E88E5 !important;
        }
        .prevalence-analysis {
          color: #1976D2 !important;
        }

        [data-bs-theme='light'] .nav-link.active .strain-analysis,
        [data-bs-theme='light'] .nav-link.active .diversity-analysis,
        [data-bs-theme='light'] .nav-link.active .prevalence-analysis {
          color: black !important;
        }

        [data-bs-theme='dark'] .nav-link.active .strain-analysis,
        [data-bs-theme='dark'] .nav-link.active .diversity-analysis,
        [data-bs-theme='dark'] .nav-link.active .prevalence-analysis {
          color: white !important;
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
  title = "INTERACTIVE_STRAIN",
  # GitHub, Documentation links
  nav_item(
    tags$a(
      href = "https://github.com/aababc1/metafun",
      target = "_blank",
      class = "nav-link",
      tags$i(class = "fab fa-github", style = "font-size:24px;color:white;"),
      " GitHub"
    )
  ),
  nav_item(
    tags$a(
      href = "https://metafun-doc.readthedocs.io/",
      target = "_blank",
      class = "nav-link",
      tags$i(class = "fas fa-book", style = "font-size:24px;color:white;"),
      " Documentation"
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
      title = "Data Input",
      accordion(
        accordion_panel(
          title = "InStrain Input Directory",
          icon = bsicons::bs_icon("folder"),
          open = TRUE,

          # Main InStrain folder input
          div(
            style = "background-color: #e7f3ff; padding: 10px; margin-bottom: 15px; border-radius: 4px; border-left: 3px solid #2196F3;",
            tags$small(
              icon("info-circle"),
              " Provide the InStrain output directory containing profile results and comparison files."
            )
          ),

          textInput("instrain_base_path",
                   tags$span(icon("folder-open"), " InStrain Output Directory:"),
                   value = Sys.getenv("INSTRAIN_DATA_PATH", "/data/instrain/"),
                   width = "100%"),

          helpText("Directory should contain: *_instrain_profile_*/output/ folders"),

          br(),
          # Load data button
          actionButton("load_data", "Load InStrain Data",
                      icon = icon("upload"),
                      class = "btn-primary btn-block blink") |>
            tooltip("Load all InStrain output files from directory"),

          br(),
          actionButton("reset_data", "Reset Data",
                      icon = icon("redo"),
                      class = "btn-secondary btn-block") |>
            tooltip("Reset all loaded data")
        ),

        accordion_panel(
          title = "Additional Files (Optional)",
          icon = bsicons::bs_icon("file-earmark-plus"),
          open = FALSE,

          textInput("gtdb_file", "GTDB Metadata File:",
                   value = Sys.getenv("GTDB_METADATA_PATH", "/data/gtdb/bac_ar_combined_metadata.tsv"),
                   width = "100%"),

          textInput("instrain_file", "InStrain Compare File (auto-detected):",
                   value = "",
                   width = "100%"),
          helpText("Leave empty to auto-detect from directory"),

          textInput("metadata_file", "Sample Metadata File:",
                   value = Sys.getenv("SAMPLE_METADATA_PATH", "/data/metadata/sample_metadata.txt"),
                   width = "100%")
        ),

        accordion_panel(
          title = "Analysis Controls",
          icon = bsicons::bs_icon("play-circle"),
          open = FALSE,

          input_task_button(
            id = "master_run_analysis",
            label = "Run Master Analysis",
            icon = icon("play")
          ) |>
          tooltip("Run analysis across all modules")
        )
      )
    ),
    
    # Main navigation panels
    navset_underline(
      nav_panel(
        tags$span("Shared Taxa Analysis", class = "strain-analysis"),
        shared_taxa_UI("shared_taxa")
      ),
      nav_panel(
        tags$span("Nucleotide Diversity", class = "diversity-analysis"),
        nucleotide_diversity_UI("nucdiv")
      ),
      # Prevalence Analysis - Disabled
      # nav_panel(
      #   tags$span("Prevalence Analysis", class = "prevalence-analysis"),
      #   prevalence_UI("prevalence")
      # ),
      nav_panel(
        tags$span("pN/pS Analysis", class = "pnps-analysis"),
        navset_card_tab(
          height = "calc(100vh - 120px)",
          nav_panel(
            "Genome-wide pN/pS",
            pnps_UI("pnps")
          ),
          nav_panel(
            "Gene pN/pS",
            gene_level_pnps_UI("gene_pnps")
          ),
          nav_panel(
            "Diversity vs pN/pS",
            diversity_pnps_correlation_UI("div_pnps_cor")
          )
        )
      )
    )
  )
)

# =============================================================================
# Server Definition
# =============================================================================

server <- function(input, output, session) {

  # Reactive values for data management
  rv <- reactiveValues(
    gtdb_data = NULL,
    sample_metadata = NULL,
    integrated_data = NULL,
    integrated_path = NULL,
    diversity_data = NULL,
    pnps_data = NULL,
    data_loaded = FALSE,
    master_trigger = 0,
    analysis_completed = FALSE
  )
  
  # Initial UI state configuration
  shinyjs::disable("master_run_analysis")
  shinyjs::disable("reset_data")

  # ==========================================================================
  # AUTO-LOAD: Load cached data on app startup (for container deployment)
  # ==========================================================================
  observe({
    # Only run once on startup
    if (is.null(rv$data_loaded) || !rv$data_loaded) {
      rds_path <- file.path(app_dir, "data/integrated_microbiome_data.rds")
      if (file.exists(rds_path)) {
        message("Auto-loading cached data: ", rds_path)
        tryCatch({
          integrated_data <- readRDS(rds_path)

          rv$gtdb_data <- integrated_data$gtdb_data
          rv$sample_metadata <- integrated_data$metadata
          rv$integrated_data <- integrated_data$ani_data
          rv$integrated_path <- rds_path
          rv$diversity_data <- integrated_data$diversity_data
          rv$pnps_data <- integrated_data$pnps_data
          rv$data_loaded <- TRUE

          shinyjs::enable("master_run_analysis")
          shinyjs::enable("reset_data")
          shinyjs::removeClass(selector = "#load_data", class = "blink")

          showNotification("Data auto-loaded successfully!", type = "message", duration = 3)
        }, error = function(e) {
          message("Auto-load failed: ", e$message)
        })
      }
    }
  }) |> bindEvent(TRUE, once = TRUE)

  # Data loading event (manual button)
observeEvent(input$load_data, {
  tryCatch({
    # Load preprocessed file if it exists
    if(file.exists("data/integrated_microbiome_data.rds")) {
      show_modal_spinner(spin = "fading-circle", text = "Loading cached data...")
      showNotification("Loading preprocessed data...", type = "message")
      integrated_data <- readRDS("data/integrated_microbiome_data.rds")
      showNotification("Preprocessed data loaded successfully!", type = "message")
    } else {
      # Run full preprocessing if cached data doesn't exist
      show_modal_spinner(spin = "fading-circle", text = "Processing data with optimized approach...")
      showNotification("Running full preprocessing (this will take several minutes)...", type = "warning")
      integrated_data <- run_full_preprocessing()
    }

    # Assign data to reactive values
    rv$gtdb_data <- integrated_data$gtdb_data
    rv$sample_metadata <- integrated_data$metadata
    rv$integrated_data <- integrated_data$ani_data
    rv$integrated_path <- "data/integrated_microbiome_data.rds"
    rv$diversity_data <- integrated_data$diversity_data
    rv$pnps_data <- integrated_data$pnps_data

    rv$data_loaded <- TRUE
    rv$analysis_completed <- FALSE
    
    # Enable controls
    shinyjs::enable("master_run_analysis")
    shinyjs::enable("reset_data")
    shinyjs::removeClass(selector = "#load_data", class = "blink")
    
    # Close sidebar
    bslib::sidebar_toggle("master-sidebar", open = FALSE)

    remove_modal_spinner()
    showNotification("All data loading completed!", type = "message")
    
  }, error = function(e) {
    remove_modal_spinner()
    showNotification(paste("Error:", e$message), type = "error")
  })
})
  
  # Reset data event
  observeEvent(input$reset_data, {
    rv$gtdb_data <- NULL
    rv$sample_metadata <- NULL
    rv$integrated_data <- NULL
    rv$diversity_data <- NULL
    rv$pnps_data <- NULL
    rv$data_loaded <- FALSE
    rv$master_trigger <- 0
    rv$analysis_completed <- FALSE
    
    shinyjs::disable("master_run_analysis")
    shinyjs::disable("reset_data")
    shinyjs::addClass(selector = "#load_data", class = "blink")
    
    showNotification("All data has been reset.", type = "message")
  })
  
  # Master analysis event
  observeEvent(input$master_run_analysis, {
    req(rv$data_loaded)
    
    showNotification("Master analysis initiated! All modules will be updated.", type = "message")
    
    rv$master_trigger <- rv$master_trigger + 1
    rv$analysis_completed <- TRUE
  })
  
  # Enable/disable master analysis based on data loading status
  observe({
    if (rv$data_loaded) {
      shinyjs::enable("master_run_analysis")
    } else {
      shinyjs::disable("master_run_analysis")
    }
  })
  
  # Module servers
  shared_taxa_Server("shared_taxa", reactive({
    data <- rv$integrated_data
    if(!is.null(data)) {
      attr(data, "source_path") <- rv$integrated_path
    }
    data
  }))

  # Nucleotide diversity module server
  nucleotide_diversity_Server("nucdiv", reactive({
    # Create the full integrated data structure that the module expects
    full_data <- list(
      diversity_data = rv$diversity_data,
      ani_data = rv$integrated_data,
      gtdb_data = rv$gtdb_data,
      metadata = rv$sample_metadata
    )
    
    if(!is.null(full_data$diversity_data)) {
      attr(full_data, "source_path") <- rv$integrated_path
    }
    full_data
  }))
  
  # prevalence_Server("prevalence", reactive(rv$integrated_data))  # Disabled

  # pN/pS analysis module servers
  pnps_Server("pnps", reactive(rv$pnps_data), reactive(rv$diversity_data))
  gene_level_pnps_Server("gene_pnps",
                         gtdb_data = reactive(rv$gtdb_data),
                         sample_metadata = reactive(rv$sample_metadata))

  # Diversity vs pN/pS correlation module
  diversity_pnps_correlation_Server("div_pnps_cor",
                                    diversity_data = reactive(rv$diversity_data),
                                    pnps_data = reactive(rv$pnps_data))
}

# launch app
options(shiny.autoreload = TRUE)

shinyApp(ui = ui, server = server)

