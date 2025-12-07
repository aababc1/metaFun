# =============================================================================
# Shared Taxa Analysis Module - Fully revised version
# =============================================================================

# Metadata detection function - For Shared Taxa
get_shared_taxa_metadata_info <- function(data) {
  # Validate input
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    message("Shared Taxa: Data is NULL, not a data frame, or has 0 rows")
    return(list(categorical = character(0), numerical = character(0)))
  }

  message("Shared Taxa: Detecting metadata from data with ", nrow(data), " rows and ", ncol(data), " columns")

  # Find _sample1, _sample2 pattern
  sample1_cols <- grep("_sample1$", colnames(data), value = TRUE)
  sample2_cols <- grep("_sample2$", colnames(data), value = TRUE)

  message("Shared Taxa: Found ", length(sample1_cols), " _sample1 columns and ", length(sample2_cols), " _sample2 columns")
  
  if (length(sample1_cols) == 0 || length(sample2_cols) == 0) {
    return(list(categorical = character(0), numerical = character(0)))
  }
  
  # Extract common metadata variable names
  meta_vars <- intersect(
    gsub("_sample1$", "", sample1_cols),
    gsub("_sample2$", "", sample2_cols)
  )
  
  # Classify each metadata type
  categorical_vars <- character(0)
  numerical_vars <- character(0)
  
  for (var in meta_vars) {
    sample1_col <- paste0(var, "_sample1")
    sample2_col <- paste0(var, "_sample2")
    
    # Combine data from both columns
    combined_data <- c(data[[sample1_col]], data[[sample2_col]])
    combined_data <- combined_data[!is.na(combined_data)]
    
    if (length(combined_data) > 0) {
      if (is.factor(combined_data) || is.character(combined_data)) {
        categorical_vars <- c(categorical_vars, var)
      } else if (is.numeric(combined_data)) {
        numerical_vars <- c(numerical_vars, var)
      }
    }
  }
  
  return(list(categorical = categorical_vars, numerical = numerical_vars))
}

# UI Module
shared_taxa_UI <- function(id) {
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
          tags$li(class = "nav-item", actionLink(ns("tab_group"), label = "Group Comparison", class = "nav-link")),
          tags$li(class = "nav-item", actionLink(ns("tab_heatmap"), label = "Heatmap", class = "nav-link")),
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

            # Heatmap availability indicator
            div(
              id = ns("heatmap_indicator"),
              style = "margin: 10px 0; padding: 8px; border-radius: 4px; font-size: 0.9em;",
              htmlOutput(ns("heatmap_availability"))
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

            # Categorical variable options - Grouping Mode
            conditionalPanel(
              condition = "output.is_categorical_var == true",
              ns = ns,
              radioButtons(ns("grouping_mode"),
                          "Grouping Mode:",
                          choices = c("Inter/Intra Groups" = "inter_intra",
                                    "Original Categories" = "original"),
                          selected = "inter_intra"),
              helpText("Inter/Intra: Groups by within-group vs between-group comparisons.",
                      "Original: Shows all category combinations separately.")
            ),

            selectInput(ns("ani_type"),
                       "ANI Type:",
                       choices = c("Population ANI" = "popANI",
                                 "Consensus ANI" = "conANI"),
                       selected = "popANI"),
            
            numericInput(ns("ani_threshold"),
                        "ANI Threshold:",
                        value = 0.99999,
                        min = 0.90, 
                        max = 1.0, 
                        step = 0.00001),
            
                # Prevalence threshold is not used in shared event analysis
                # sliderInput(ns("min_prevalence"),
                #            "Minimum Prevalence (%):",
                #            min = 0, max = 100, value = 10, step = 5),
            
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
        # Message before data load
        conditionalPanel(
          condition = "output.data_loaded == false",
          ns = ns,
          div(
            style = "color: red; font-weight: bold; font-size: 20px; text-align: center; margin: 50px;",
            "Load integrated data first"
          )
        ),
        
        # Content after data load
        conditionalPanel(
          condition = "output.data_loaded == true",
          ns = ns,
          
              tabsetPanel(
                id = ns("main_tabs"),
                selected = "Overview",
                tabPanel("Overview",
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
                      title = "Total Samples",
                      value = textOutput(ns("overview_total_samples_text")),
                      showcase = bs_icon("stack"),
                      theme = "mint"
                    ),
                    value_box(
                      title = "Unique Sample Pairs",
                      value = textOutput(ns("overview_unique_pairs_text")),
                      showcase = bs_icon("diagram-2"),
                      theme = "mint"
                    ),
                    value_box(
                      title = "Total Shared Events",
                      value = textOutput(ns("overview_total_events_text")),
                      showcase = bs_icon("arrow-left-right"),
                      theme = "mint"
                    ),
                    value_box(
                      title = "Avg Events per Pair",
                      value = textOutput(ns("overview_events_per_pair_text")),
                      showcase = bs_icon("graph-up"),
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
                        create_plot_download_btn(ns, "overview_events")
                      ),
                      plotlyOutput(ns("overview_events_plot"), height = "400px")
                    ),
                    column(
                      width = 6,
                      div(
                        style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                        create_plot_download_btn(ns, "overview_taxonomy")
                      ),
                      plotlyOutput(ns("overview_taxonomy_plot"), height = "400px")
                    )
                  )
                ),
                
                tabPanel("Group Comparison",
                  div(
                    style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                    create_plot_download_btn(ns, "group_comparison")
                  ),
                  plotlyOutput(ns("group_comparison_plot"), height = "600px"),
                  br(),
                  # Statistical comparison section (only for categorical variables)
                  conditionalPanel(
                    condition = "output.is_categorical_var == true",
                    ns = ns,
                    card(
                      card_header("Statistical Comparison"),
                      card_body(
                        htmlOutput(ns("statistical_summary")),
                        br(),
                        DT::dataTableOutput(ns("statistical_results_table"))
                      )
                    ),
                    br()
                  ),
                  div(id = ns("clicked_data_section"),
                    card(
                      card_header("Selected Point Details"),
                      card_body(
                        downloadButton(ns("download_selected_details"), "Download CSV", class = "btn-secondary"),
                        br(),
                        DT::dataTableOutput(ns("clicked_point_table"))
                      )
                    )
                  ),
                  br(),
                  div(
                    style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
                    h4("Taxonomy Count Distribution", style = "margin: 0;"),
                    create_plot_download_btn(ns, "density")
                  ),
                  plotlyOutput(ns("density_plot"), height = "400px"),
                  br(),
                  conditionalPanel(
                    condition = "output.is_categorical_var == true",
                    ns = ns,
                    div(
                      style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
                      h4("Stacked Bar Chart by Category", style = "margin: 0;"),
                      create_plot_download_btn(ns, "stacked_bar")
                    ),
                    plotlyOutput(ns("stacked_bar_plot"), height = "400px")
                  )
                ),
                
                
                tabPanel("Heatmap",
                  div(
                    style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                    create_plot_download_btn(ns, "heatmap")
                  ),
                  plotlyOutput(ns("heatmap_plot"), height = "600px"),
                  br(),
                  # Click detail area below Heatmap
                  div(id = ns("clicked_data_section_hm"), 
                    card(
                      card_header("Selected Cell Details"),
                      card_body(
                        DT::dataTableOutput(ns("clicked_point_table_hm"))
                      )
                    )
                  ),
                  br()
                ),
                
                tabPanel("Data Table",
                  DT::dataTableOutput(ns("raw_table"))
                )
              )
            )
          )
        )
      )
    )
  )
}

# Function to load genome_info data from instrain folder
## Simple in-memory caches to avoid repeated disk reads and joins
gtdb_cache_env <- new.env(parent = emptyenv())
genome_cache_env <- new.env(parent = emptyenv())

get_gtdb_cached <- function(gtdb_file = "data/filtered_gtdb_metadata.rds") {
  if (!is.null(gtdb_cache_env$gtdb)) {
    return(gtdb_cache_env$gtdb)
  }
  if (!file.exists(gtdb_file)) {
    message("GTDB file not found: ", gtdb_file)
    gtdb_cache_env$gtdb <- data.frame()
    return(gtdb_cache_env$gtdb)
  }
  message("Loading GTDB metadata (cached): ", gtdb_file)
  gtdb_df <- readRDS(gtdb_file)
  # Normalize accession for robust joins (support both digits-only and GC[FA]_ formats)
  if (!"clean_accession" %in% colnames(gtdb_df) && "accession" %in% colnames(gtdb_df)) {
    gtdb_df$clean_accession <- gtdb_df$accession
  }
  # Extract GC[FA]_#########.# from any string
  extract_root <- function(x) {
    x <- as.character(x)
    m <- regexpr("GC[FA]_\\d+\\.\\d+", x, perl = TRUE)
    v <- regmatches(x, m)
    v[is.na(v) | v == ""] <- NA_character_
    return(v)
  }
  # If we have an 'accession' column, root is from it; otherwise try clean_accession
  gtdb_df$root_accession <- if ("accession" %in% colnames(gtdb_df)) extract_root(gtdb_df$accession) else extract_root(gtdb_df$clean_accession)
  # accession_core is digits-only; prefer existing clean_accession if it is already digits-only
  is_digits <- function(s) grepl("^\\d+\\.\\d+$", as.character(s))
  gtdb_df$accession_core <- ifelse(is_digits(gtdb_df$clean_accession),
                                   gtdb_df$clean_accession,
                                   sub("^GC[FA]_", "", as.character(gtdb_df$root_accession)))
  gtdb_cache_env$gtdb <- gtdb_df
  return(gtdb_cache_env$gtdb)
}

load_genome_info_for_samples <- function(instrain_base_path, sample_ids, target_taxon) {
  tryCatch({
    message("Loading genome_info for samples: ", paste(sample_ids, collapse = ", "))
    
    # Load GTDB data cache
    gtdb_data <- get_gtdb_cached()
    if (nrow(gtdb_data) == 0) {
      return(data.frame())
    }
    message("GTDB data available (cached): ", nrow(gtdb_data), " records")
    
    all_genome_data <- list()
    
    for(sample_id in sample_ids) {
      # Find instrain directory (with exact pattern)
      instrain_dirs <- list.files(instrain_base_path, 
                                 pattern = paste0("^", sample_id, "_fastp_hg38.sorted_instrain_profile_0.92$"), 
                                 full.names = FALSE)
      
      if(length(instrain_dirs) > 0) {
        dir_name <- instrain_dirs[1]  # Use first matching directory
        genome_file <- file.path(instrain_base_path, dir_name, "output", 
                               paste0(dir_name, "_genome_info.tsv"))
        
        if(file.exists(genome_file)) {
          message("Found genome file: ", genome_file)
          
          # Check cache: store genome_with_taxonomy per sample_id
          if (exists(sample_id, envir = genome_cache_env, inherits = FALSE)) {
            genome_with_taxonomy <- get(sample_id, envir = genome_cache_env, inherits = FALSE)
          } else {
            genome_data <- read.delim(genome_file, sep = "\t", header = TRUE, 
                                     stringsAsFactors = FALSE)
            
            message("Genome data loaded, rows: ", nrow(genome_data))
            
            if(nrow(genome_data) > 0) {
              # Extract root accession from filename robustly (e.g., GCA_004000625.1_ASM... or ..._genomic.fna.gz)
              extract_root <- function(x) {
                x <- as.character(x)
                m <- regexpr("GC[FA]_\\d+\\.\\d+", x, perl = TRUE)
                v <- regmatches(x, m)
                v[is.na(v) | v == ""] <- NA_character_
                return(v)
              }
              genome_data$root_accession <- extract_root(basename(genome_data$genome))
              # Backward-compatible fields
              genome_data$clean_genome <- gsub("_genomic\\.fna(\\.gz)?$", "", genome_data$genome)
              genome_data$clean_accession <- ifelse(is.na(genome_data$root_accession) | genome_data$root_accession == "",
                                                    genome_data$clean_genome,
                                                    genome_data$root_accession)
              genome_data$accession_core <- sub("^GC[FA]_", "", as.character(genome_data$clean_accession))
              
              message("Sample genome accessions: ", paste(head(genome_data$clean_accession, 3), collapse = ", "))
              message("GTDB clean_accession format: ", paste(head(gtdb_data$clean_accession, 3), collapse = ", "))
              
              # Join with GTDB data to get species information
              # 1) join by digits-only core accession
genome_with_taxonomy <- genome_data %>%
  left_join(gtdb_data %>% select(accession_core, domain, phylum, class, order, family, genus, species),
            by = c("accession_core" = "accession_core"))

# If still NA, match clean_accession directly with GTDB accession
na_rows <- is.na(genome_with_taxonomy$species)
if (any(na_rows)) {
  genome_with_taxonomy[na_rows, ] <- genome_with_taxonomy[na_rows, ] %>%
    select(-domain, -phylum, -class, -order, -family, -genus, -species) %>%
    left_join(gtdb_data %>% 
                mutate(accession_clean = gsub("^RS_|^GB_", "", accession)) %>%
                select(accession_clean, domain, phylum, class, order, family, genus, species),
              by = c("clean_accession" = "accession_clean"))
}

# If finally still NA, try with root_accession
na_rows <- is.na(genome_with_taxonomy$species)
if (any(na_rows) && "root_accession" %in% colnames(gtdb_data)) {
  genome_with_taxonomy[na_rows, ] <- genome_with_taxonomy[na_rows, ] %>%
    select(-domain, -phylum, -class, -order, -family, -genus, -species) %>%
    left_join(gtdb_data %>% select(root_accession, domain, phylum, class, order, family, genus, species),
              by = c("clean_accession" = "root_accession"))
}
                # Create taxonomy columns if they don't exist and handle safely
                for(colnm in c("domain","phylum","class","order","family","genus","species")){
                  if(!colnm %in% names(genome_with_taxonomy)) genome_with_taxonomy[[colnm]] <- NA_character_
                }
                genome_with_taxonomy <- genome_with_taxonomy %>%
                  mutate(
                    domain = ifelse(is.na(domain), "Unknown", domain),
                    phylum = ifelse(is.na(phylum), "Unknown", phylum),
                    class = ifelse(is.na(class), "Unknown", class),
                    order = ifelse(is.na(order), "Unknown", order),
                    family = ifelse(is.na(family), "Unknown", family),
                    genus = ifelse(is.na(genus), "Unknown", genus),
                    species = ifelse(is.na(species), "Unknown", species)
                  )
              
              assign(sample_id, genome_with_taxonomy, envir = genome_cache_env)
            } else {
              genome_with_taxonomy <- data.frame()
            }
          }
          
          if (nrow(genome_with_taxonomy) > 0) {
            message("After GTDB join (cached) - rows with taxonomy: ", sum(!is.na(genome_with_taxonomy$species)))
            message("Sample species found: ", paste(unique(genome_with_taxonomy$species[!is.na(genome_with_taxonomy$species)])[1:3], collapse = ", "))
            
            # Filter genomes matching target_taxon
            filtered_genomes <- genome_with_taxonomy %>%
              filter(
                grepl(target_taxon, species, ignore.case = TRUE) |
                grepl(target_taxon, genus, ignore.case = TRUE) |
                grepl(target_taxon, family, ignore.case = TRUE)
              )
            
            if(nrow(filtered_genomes) > 0) {
              filtered_data <- filtered_genomes %>%
                mutate(
                  sample_id = sample_id,
                  # Keep GTDB tokens as-is (d__/p__/...; suffixes preserved)
                  domain = ifelse(is.na(domain) | domain == "", "Unknown", domain),
                  phylum = ifelse(is.na(phylum) | phylum == "", "Unknown", phylum),
                  class = ifelse(is.na(class) | class == "", "Unknown", class),
                  order = ifelse(is.na(order) | order == "", "Unknown", order),
                  family = ifelse(is.na(family) | family == "", "Unknown", family),
                  genus = ifelse(is.na(genus) | genus == "", "Unknown", genus),
                  species = ifelse(is.na(species) | species == "", "Unknown", species)
                ) %>%
                mutate(`Full Taxonomy` = paste(domain, phylum, class, order, family, genus, species, sep = " | ")) %>%
                select(
                  `Sample ID` = sample_id,
                  `Genome` = genome,
                  `Accession` = clean_genome,
                  `Clean Accession` = clean_accession,
                  `Full Taxonomy`,
                  `Domain` = domain,
                  `Phylum` = phylum,
                  `Class` = class,
                  `Order` = order,
                  `Family` = family,
                  `Genus` = genus,
                  `Species` = species,
                  `Coverage` = coverage,
                  `Breadth` = breadth,
                  `Nucleotide Diversity` = nucl_diversity,
                  `SNV Count` = SNV_count,
                  `Genome Size` = length,
                  `Coverage Median` = coverage_median,
                  `Coverage Std` = coverage_std,
                  `True Scaffolds` = true_scaffolds,
                  `Detected Scaffolds` = detected_scaffolds
                )
              
              all_genome_data[[sample_id]] <- filtered_data
              message("Added genome data for sample: ", sample_id, " (", nrow(filtered_data), " records)")
              message("Species found: ", paste(unique(filtered_data$Species), collapse = ", "))
            } else {
              message("No genomes found matching taxon: ", target_taxon, " for sample: ", sample_id)
            }
          }
        } else {
          message("Genome file not found: ", genome_file)
        }
      } else {
        message("No instrain directory found for sample: ", sample_id)
      }
    }
    
    if(length(all_genome_data) > 0) {
      combined_data <- do.call(rbind, all_genome_data)
      message("Loaded genome_info data: ", nrow(combined_data), " records")
      return(combined_data)
    } else {
      message("No genome_info data found for the samples")
      return(data.frame())
    }
    
  }, error = function(e) {
    message("Error loading genome_info: ", e$message)
    return(data.frame())
  })
}

# Server Module
shared_taxa_Server <- function(id, integrated_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive value for clicked data
    clicked_data_rv <- reactiveValues(data = NULL)
    
    # Error handling function
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

        # Check if data is a data frame with rows (robust validation)
        if(is.data.frame(data) && nrow(data) > 0) {
          message("Shared Taxa: Detecting metadata from data frame with ", nrow(data), " rows")
          get_shared_taxa_metadata_info(data)
        } else if(is.list(data) && "ani_data" %in% names(data) && !is.null(data$ani_data)) {
          # If data is a list structure, extract ani_data
          ani_data <- data$ani_data
          message("Shared Taxa: Detecting metadata from ani_data, rows: ", nrow(ani_data))
          get_shared_taxa_metadata_info(ani_data)
        } else {
          message("Shared Taxa: No valid data found for metadata detection")
          return(list(categorical = character(0), numerical = character(0)))
        }
      }, error = function(e) {
        showNotification(paste("Metadata detection error:", e$message), type = "warning")
        return(list(categorical = character(0), numerical = character(0)))
      })
    })
    
    # Heatmap click -> populate details (improved version)
    observeEvent(event_data("plotly_click", source = "heatmap_plot"), {
      click <- event_data("plotly_click", source = "heatmap_plot")
      req(click)

      tryCatch({
        cat("=== HEATMAP CLICK EVENT ===\n")
        cat("Click data:\n")
        print(click)

        data <- processed_data()
        req(data)

        # Extract clicked information
        clicked_taxon <- as.character(click$x)
        clicked_y_index <- click$y  # 0-based index

        cat("Clicked taxon:", clicked_taxon, "\n")
        cat("Clicked Y index:", clicked_y_index, "\n")

        # Map Y-axis labels from current_heatmap_info
        if(exists("current_heatmap_info") && !is.null(current_heatmap_info)) {
          available_sharing_types <- current_heatmap_info$sharing_types
          cat("Available sharing types in heatmap:", paste(available_sharing_types, collapse = " | "), "\n")

          # Convert Y index to actual sharing_type
          clicked_sharing <- if(clicked_y_index >= 0 && clicked_y_index < length(available_sharing_types)) {
            available_sharing_types[clicked_y_index + 1]  # R is 1-based
          } else {
            # Fallback: use customdata or first sharing_type
            if(!is.null(click$customdata)) {
              as.character(click$customdata)
            } else {
              available_sharing_types[1]
            }
          }
        } else {
          # Fallback: extract directly from data$shared_events
          unique_sharing_types <- sort(unique(data$shared_events$sharing_type))
          cat("Using fallback sharing types:", paste(unique_sharing_types, collapse = " | "), "\n")

          clicked_sharing <- if(clicked_y_index >= 0 && clicked_y_index < length(unique_sharing_types)) {
            unique_sharing_types[clicked_y_index + 1]
          } else {
            # direct value fallback
            as.character(click$y)
          }
        }

        cat("Mapped to sharing type:", clicked_sharing, "\n")

        # Filter shared events
        shared_events_data <- data$shared_events %>%
          filter(
            taxon == clicked_taxon,
            sharing_type == clicked_sharing
          ) %>%
          distinct(sample1, sample2, species, .keep_all = TRUE)

        cat("Found", nrow(shared_events_data), "events for cell (", clicked_taxon, ",", clicked_sharing, ")\n")

        if(nrow(shared_events_data) == 0) {
          showNotification("No events for the selected cell", type = "warning")
          clicked_data_rv$data <- data.frame()
          shinyjs::hide("clicked_data_section_hm")
          return()
        }

        # Load genome info
        instrain_base_path <- Sys.getenv("INSTRAIN_DATA_PATH", "/data/instrain/")

        # Generate coverage summary
        coverage_list <- lapply(seq_len(nrow(shared_events_data)), function(i) {
          row <- shared_events_data[i, ]

          sample1_genome_data <- load_genome_info_for_samples(
            instrain_base_path, row$sample1, row$species
          )
          sample2_genome_data <- load_genome_info_for_samples(
            instrain_base_path, row$sample2, row$species
          )

          sample1_info <- if(nrow(sample1_genome_data) > 0) {
            sample1_genome_data %>%
              slice_max(`Coverage`, n = 1) %>%
              transmute(
                `Full Taxonomy_Sample1` = `Full Taxonomy`,
                `Coverage_Sample1` = Coverage,
                `Breadth_Sample1` = Breadth,
                `Nucleotide Diversity_Sample1` = `Nucleotide Diversity`,
                `SNV Count_Sample1` = `SNV Count`,
                `Genome Size_Sample1` = `Genome Size`
              )
          } else {
            tibble::tibble(
              `Full Taxonomy_Sample1` = NA_character_,
              `Coverage_Sample1` = NA_real_,
              `Breadth_Sample1` = NA_real_,
              `Nucleotide Diversity_Sample1` = NA_real_,
              `SNV Count_Sample1` = NA_real_,
              `Genome Size_Sample1` = NA_real_
            )
          }

          sample2_info <- if(nrow(sample2_genome_data) > 0) {
            sample2_genome_data %>%
              slice_max(`Coverage`, n = 1) %>%
              transmute(
                `Full Taxonomy_Sample2` = `Full Taxonomy`,
                `Coverage_Sample2` = Coverage,
                `Breadth_Sample2` = Breadth,
                `Nucleotide Diversity_Sample2` = `Nucleotide Diversity`,
                `SNV Count_Sample2` = `SNV Count`,
                `Genome Size_Sample2` = `Genome Size`
              )
          } else {
            tibble::tibble(
              `Full Taxonomy_Sample2` = NA_character_,
              `Coverage_Sample2` = NA_real_,
              `Breadth_Sample2` = NA_real_,
              `Nucleotide Diversity_Sample2` = NA_real_,
              `SNV Count_Sample2` = NA_real_,
              `Genome Size_Sample2` = NA_real_
            )
          }

          join_key <- paste(row$sample1, row$sample2, row$species, sep = "|")

          tibble::tibble(join_key) %>%
            bind_cols(sample1_info) %>%
            bind_cols(sample2_info)
        })

        coverage_summary <- if(length(coverage_list) > 0) {
          dplyr::bind_rows(coverage_list)
        } else {
          tibble::tibble(join_key = character(0))
        }

        # Generate detailed data
        detailed_data <- shared_events_data %>%
          mutate(join_key = paste(sample1, sample2, species, sep = "|")) %>%
          left_join(coverage_summary, by = "join_key") %>%
          select(
            `Sample 1` = sample1,
            `Sample 2` = sample2,
            `Taxon` = taxon,
            `Species` = species,
            `Group 1` = group1,
            `Group 2` = group2,
            `Sharing Type` = sharing_type,
            `Combination` = metadata_combination,
            starts_with("Full"),
            ends_with("Sample1"),
            ends_with("Sample2")
          )

        # Store in reactive value
        clicked_data_rv$data <- detailed_data

        # Render table
        output$clicked_point_table_hm <- DT::renderDataTable({
          DT::datatable(
            detailed_data,
            options = list(
              pageLength = 10,
              scrollX = TRUE,
              dom = 'Bfrtip',
              buttons = c('copy', 'csv', 'excel')
            ),
            extensions = 'Buttons',
            caption = paste0("Shared Events for ", clicked_taxon, " (", clicked_sharing, ") - ",
                            nrow(detailed_data), " events"),
            rownames = FALSE
          )
        })

        shinyjs::show("clicked_data_section_hm")
        showNotification(
          paste0("Loaded ", nrow(detailed_data), " events for ", clicked_taxon, " (", clicked_sharing, ")"),
          type = "message"
        )

      }, error = function(e) {
        cat("=== HEATMAP CLICK ERROR ===\n")
        cat("Error:", e$message, "\n")
        showNotification(paste("Heatmap click error:", e$message), type = "error")
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
      if (input$metadata_var %in% meta_info$categorical) {
        return("categorical")
      } else if (input$metadata_var %in% meta_info$numerical) {
        return("numerical")
      }
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

    # Heatmap availability indicator
    output$heatmap_availability <- renderText({
      if(is.null(input$metadata_var) || input$metadata_var == "") {
        return(HTML("<div style='color: #6c757d; background: #f8f9fa; padding: 10px; border-radius: 4px;'>
                     <strong>🔥 Heatmap:</strong> Select a variable above to check availability
                   </div>"))
      }

      tryCatch({
        var_type <- current_var_type()
        if(var_type == "categorical") {
          return(HTML("<div style='color: #155724; background: #d4edda; border: 1px solid #c3e6cb; padding: 10px; border-radius: 4px;'>
                       <strong>✅ Heatmap Available:</strong> Categorical variable selected - heatmap will work!
                     </div>"))
        } else if(var_type == "numerical") {
          return(HTML("<div style='color: #721c24; background: #f8d7da; border: 1px solid #f5c6cb; padding: 10px; border-radius: 4px;'>
                       <strong>❌ Heatmap Unavailable:</strong> Numerical variables not supported in heatmap<br>
                       <small>Please select a categorical variable (e.g., disease_group, country, host_sex)</small>
                     </div>"))
        } else {
          return(HTML("<div style='color: #856404; background: #fff3cd; border: 1px solid #ffeaa7; padding: 10px; border-radius: 4px;'>
                       <strong>⚠️ Heatmap Status:</strong> Unknown variable type
                     </div>"))
        }
      }, error = function(e) {
        return(HTML("<div style='color: #721c24; background: #f8d7da; padding: 10px; border-radius: 4px;'>
                     <strong>❌ Error:</strong> Could not determine heatmap availability
                   </div>"))
      })
    })

    # UI control - when data is available
    observe({
      req(integrated_data())

      # Additional validation: ensure we have valid metadata info
      meta_info <- metadata_info()
      if(is.null(meta_info)) {
        return()
      }

      shinyjs::hide("analysis_message")

      all_vars <- c(meta_info$categorical, meta_info$numerical)
      
      if (length(all_vars) > 0) {
        # Update metadata selection
        choices <- setNames(all_vars, tools::toTitleCase(gsub("_", " ", all_vars)))
        
        current_selection <- input$metadata_var
        if (is.null(current_selection) || !current_selection %in% all_vars) {
          # If first categorical exists use that, otherwise use first variable
          default_selection <- if(length(meta_info$categorical) > 0) meta_info$categorical[1] else all_vars[1]
        } else {
          default_selection <- current_selection
        }
        
        updateSelectInput(session, "metadata_var",
                         choices = choices,
                         selected = default_selection)
        
        # Enable buttons
        shinyjs::enable("update_analysis")
        shinyjs::enable("metadata_var")
      } else {
        shinyjs::disable("metadata_var")
        shinyjs::disable("update_analysis")
        showNotification("No valid metadata variables found", type = "warning")
      }
    })

    # Header tab click handlers to control main tabset (no ns() here; id matches directly)
    observeEvent(input$tab_overview, ignoreInit = TRUE, {
      updateTabsetPanel(session, "main_tabs", selected = "Overview")
    })
    observeEvent(input$tab_group, ignoreInit = TRUE, {
      updateTabsetPanel(session, "main_tabs", selected = "Group Comparison")
    })
    observeEvent(input$tab_heatmap, ignoreInit = TRUE, {
      cat("=== HEATMAP DEBUG: Tab heatmap clicked ===\n")
      cat("Current metadata_var:", input$metadata_var, "\n")
      cat("Current tax_level:", input$tax_level, "\n")
      if (exists("processed_data") && !is.null(tryCatch(processed_data(), error = function(e) NULL))) {
        data <- processed_data()
        cat("Processed data exists, var_type:", ifelse(is.null(data$var_type), "NULL", data$var_type), "\n")
      } else {
        cat("No processed data available\n")
      }
      updateTabsetPanel(session, "main_tabs", selected = "Heatmap")
    })
    observeEvent(input$tab_data, ignoreInit = TRUE, {
      updateTabsetPanel(session, "main_tabs", selected = "Data Table")
    })
    
    # Update taxon selection
    observe({
      req(processed_data())
      
      data <- processed_data()
      if(!is.null(data)) {
        # Extract unique taxon list from analyzed data
        unique_taxa <- unique(data$prevalence_data$taxon)
        unique_taxa <- unique_taxa[!is.na(unique_taxa) & unique_taxa != ""]
        
        if(length(unique_taxa) > 0) {
          # Add "All" option at the front
          choices <- c("All" = "all", setNames(unique_taxa, unique_taxa))
          
          updateSelectizeInput(
            session = session,
            inputId = "selected_taxon",
            choices = choices,
            selected = "all",  # Select "All" by default
            options = list(
              placeholder = paste('Select', input$tax_level, 'or All'),
              maxOptions = 500
            )
          )
        }
      }
    })
    
    # Display selected taxon information
    output$selected_taxon_info <- renderText({
      req(input$selected_taxon, processed_data())
      
      data <- processed_data()
      if(is.null(data) || is.null(input$selected_taxon)) {
        return(HTML("<span style='color: #6c757d;'>No taxon selected</span>"))
      }
      
      tryCatch({
        if(input$selected_taxon == "all") {
          # Display all data information
          total_taxa <- length(unique(data$prevalence_data$taxon))
          total_groups <- length(unique(data$prevalence_data$group))
          avg_prevalence <- mean(data$prevalence_data$prevalence, na.rm = TRUE)
          
          info_text <- sprintf(
            "<strong>All Taxa</strong><br>
            Total Taxa: %d | Groups: %d | Avg Prevalence: %.1f%%",
            total_taxa,
            total_groups,
            avg_prevalence
          )
          
          HTML(info_text)
        } else {
          # Extract information for selected specific taxon
          taxon_data <- data$prevalence_data %>%
            filter(taxon == input$selected_taxon)
          
          if(nrow(taxon_data) > 0) {
            total_groups <- length(unique(taxon_data$group))
            avg_prevalence <- mean(taxon_data$prevalence, na.rm = TRUE)
            total_species <- sum(taxon_data$unique_species, na.rm = TRUE)
            
            info_text <- sprintf(
              "<strong>%s</strong><br>
              Groups: %d | Avg Prevalence: %.1f%% | Unique Species: %d",
              input$selected_taxon,
              total_groups,
              avg_prevalence,
              total_species
            )
            
            HTML(info_text)
          } else {
            HTML("<span style='color: #dc3545;'>Taxon not found in analysis results</span>")
          }
        }
      }, error = function(e) {
        HTML("<span style='color: #dc3545;'>Error loading taxon information</span>")
      })
    })
    
    # UI control - when data is not available
    observe({
      if (is.null(integrated_data())) {
        shinyjs::show("analysis_message")
        shinyjs::disable("metadata_var")
        shinyjs::disable("update_analysis")
      } else {
        # Move to Overview tab when data is loaded
        updateTabsetPanel(session, "main_tabs", selected = "Overview")
      }
    })
    
    # Process analysis data
    processed_data <- eventReactive(input$update_analysis, {
      req(integrated_data(), input$metadata_var)
      
      tryCatch({
        message("Starting analysis...")
        data <- integrated_data()
        message("Data loaded, rows: ", nrow(data))
        
        # ANI filtering - use data generated by preprocessing_scripts.R
        ani_col <- input$ani_type
        message("ANI column: ", ani_col)
        message("ANI threshold: ", input$ani_threshold)
        
        # Use ani_data - data generated by preprocessing_scripts.R
        if("ani_data" %in% names(data) && !is.null(data$ani_data)) {
          filtered_data <- data$ani_data[!is.na(data$ani_data[[ani_col]]) & data$ani_data[[ani_col]] >= input$ani_threshold, ]
          message("Using ani_data, rows after filtering: ", nrow(filtered_data))
        } else {
        filtered_data <- data[!is.na(data[[ani_col]]) & data[[ani_col]] >= input$ani_threshold, ]
          message("Using main data, rows after filtering: ", nrow(filtered_data))
        }
        
        # Selected taxonomic level and metadata variable
        tax_col <- input$tax_level
        
        # Debugging: check available columns
        message("=== AVAILABLE COLUMNS DEBUGGING ===")
        message("All columns in data: ", paste(colnames(data), collapse = ", "))
        message("Selected metadata variable: ", input$metadata_var)
        
        # Select columns matching data structure generated by preprocessing_scripts.R
        # ani_data already has metadata joined with _sample1, _sample2 suffix
        # Verify that sample-specific metadata columns for selected variable actually exist
        meta_col_1 <- paste0(input$metadata_var, "_sample1")
        meta_col_2 <- paste0(input$metadata_var, "_sample2")
        if(!all(c(meta_col_1, meta_col_2) %in% colnames(filtered_data))) {
          # disease_group is a commonly used default variable, so check additionally
          if(input$metadata_var == "disease_group" &&
             all(c("disease_group_sample1", "disease_group_sample2") %in% colnames(filtered_data))) {
            meta_col_1 <- "disease_group_sample1"
            meta_col_2 <- "disease_group_sample2"
          } else {
            message("Missing metadata columns for variable '", input$metadata_var,
                    "': expected ", meta_col_1, ", ", meta_col_2)
            showNotification(paste0("Missing metadata columns for variable '",
                                    input$metadata_var, "'"), type = "error")
            validate(need(FALSE, paste0("Missing metadata columns for variable '",
                                        input$metadata_var, "'")))
          }
        }
        message("Using sample-specific metadata columns: ", meta_col_1, ", ", meta_col_2)
        
        # Numerical variable processing - use Original only
        if (current_var_type() == "numerical" && !is.null(input$numerical_analysis_type)) {
          if (input$numerical_analysis_type == "original") {
            # Keep original numerical data (no transformation)
            message("Using original numerical data without transformation")
          }
        }
        
        # Event-based analysis: each row = one shared event (sample1 ↔ sample2)
        message("Tax column: ", tax_col)
        message("Meta columns: ", meta_col_1, ", ", meta_col_2)
        
        # Debugging: check original data structure
        message("=== ORIGINAL DATA STRUCTURE DEBUGGING ===")
        message("Filtered data columns: ", paste(colnames(filtered_data), collapse = ", "))
        message("Sample filtered data:")
        print(head(filtered_data %>% select(clean_name1, clean_name2, all_of(meta_col_1), all_of(meta_col_2)), 5))
        message("Unique values in meta_col_1: ", paste(unique(filtered_data[[meta_col_1]]), collapse = ", "))
        message("Unique values in meta_col_2: ", paste(unique(filtered_data[[meta_col_2]]), collapse = ", "))
        
        # Check data structure generated by preprocessing_scripts.R
        message("=== PREPROCESSING DATA STRUCTURE DEBUGGING ===")
        if("diversity_data" %in% names(data)) {
          message("Diversity data available: ", nrow(data$diversity_data), " records")
          message("Diversity data columns: ", paste(colnames(data$diversity_data), collapse = ", "))
        }
        if("ani_data" %in% names(data)) {
          message("ANI data available: ", nrow(data$ani_data), " records")
          message("ANI data columns: ", paste(colnames(data$ani_data), collapse = ", "))
        }
        
        # Generate shared event data - use directly from ani_data
        message("=== USING ANI DATA WITH PRE-JOINED METADATA ===")
        message("ANI data columns: ", paste(colnames(filtered_data), collapse = ", "))
        message("Using metadata columns: ", meta_col_1, ", ", meta_col_2)
        
        # Get grouping mode early for raw_display (will be used again for shared_events)
        grouping_mode <- if(!is.null(input$grouping_mode)) input$grouping_mode else "inter_intra"

        # Generate raw display data same as Data Table (save for reuse)
        raw_display <- filtered_data %>%
          mutate(
            group1 = .data[[meta_col_1]],
            group2 = .data[[meta_col_2]],
            taxon = .data[[tax_col]]
          ) %>%
          mutate(
            group1_chr = trimws(as.character(group1)),
            group2_chr = trimws(as.character(group2)),
            category_combination = paste0(
              pmin(group1_chr, group2_chr), " ↔ ",
              pmax(group1_chr, group2_chr)
            )
          )

        # Apply grouping logic based on mode
        if(grouping_mode == "inter_intra") {
          raw_display <- raw_display %>%
            mutate(
              sharing_type = if_else(group1_chr == group2_chr, "Intra-group", "Inter-group"),
              sharing_label = if_else(
                group1_chr == group2_chr,
                paste0("Intra-group (", group1_chr, ")"),
                paste0("Inter-group (", pmin(group1_chr, group2_chr), " ↔ ", pmax(group1_chr, group2_chr), ")")
              ),
              metadata_combination = category_combination
            )
        } else {
          raw_display <- raw_display %>%
            mutate(
              sharing_type = category_combination,
              sharing_label = category_combination,
              metadata_combination = category_combination
            )
        }

        raw_display <- raw_display %>%
          select(-group1_chr, -group2_chr)

        front_cols_display <- c("clean_name1", "clean_name2", "species", tax_col,
                                 "popANI", "conANI", "sharing_type", "sharing_label", "metadata_combination")
        front_cols_display <- intersect(front_cols_display, colnames(raw_display))
        other_cols_display <- setdiff(colnames(raw_display), front_cols_display)
        raw_display <- raw_display[, c(front_cols_display, other_cols_display), drop = FALSE]

        # Use directly from ani_data (metadata already joined)
        shared_events <- filtered_data %>%
          filter(!is.na(.data[[tax_col]]) & 
                 !is.na(.data[[meta_col_1]]) & 
                 !is.na(.data[[meta_col_2]])) %>%
          select(
            sample1 = clean_name1,
            sample2 = clean_name2,
            taxon = all_of(tax_col),
            group1 = all_of(meta_col_1),
            group2 = all_of(meta_col_2),
            species = species
          ) %>%
          mutate(
            group1 = trimws(as.character(group1)),
            group2 = trimws(as.character(group2))
          )

        # Grouping mode was already set above for raw_display
        message("Grouping mode: ", grouping_mode)

        shared_events <- shared_events %>%
          # Create metadata combination category - process differently based on mode
          mutate(
            # Create sorted category combination (collapse symmetric pairs)
            category_combination = paste0(
              pmin(group1, group2), " ↔ ",
              pmax(group1, group2)
            )
          )

        # Apply different grouping logic based on mode
        if(grouping_mode == "inter_intra") {
          # Original behavior: Group as Inter/Intra
          shared_events <- shared_events %>%
            mutate(
              # If selected metadata values are same: Intra-group, different: Inter-group
              sharing_type = if_else(group1 == group2, "Intra-group", "Inter-group"),
              # Label with values specified
              sharing_label = if_else(
                group1 == group2,
                paste0("Intra-group (", group1, ")"),
                paste0("Inter-group (", pmin(group1, group2), " ↔ ", pmax(group1, group2), ")")
              ),
              # Stable sorted combination label for matching and display
              metadata_combination = category_combination
            )
        } else {
          # Original categories mode: Show each category combination
          shared_events <- shared_events %>%
            mutate(
              # Use the actual category combination as sharing_type
              sharing_type = category_combination,
              sharing_label = category_combination,
              metadata_combination = category_combination
            )
        }

        shared_events <- shared_events %>%
          # Add logs for debugging
          mutate(
            debug_group1 = group1,
            debug_group2 = group2,
            debug_group_equal = as.character(group1) == as.character(group2),
            debug_sharing_type = sharing_type,
            debug_metadata_combination = metadata_combination,
            debug_grouping_mode = grouping_mode
          ) %>%
          # Remove duplicates so sample1, sample2, species combination is unique
          distinct(sample1, sample2, species, .keep_all = TRUE)
        
        message("Shared events rows: ", nrow(shared_events))
        message("Unique metadata combinations: ", n_distinct(shared_events$metadata_combination))
        
        # Debugging: check actual data
        message("=== DEBUGGING: Actual Data Check ===")
        sample_data <- shared_events %>% 
          select(debug_group1, debug_group2, debug_group_equal, debug_sharing_type, debug_metadata_combination) %>%
          head(10)
        print(sample_data)
        
        # Check Within-group and Cross-group distribution
        sharing_type_counts <- table(shared_events$debug_sharing_type)
        message("Sharing type distribution:")
        print(sharing_type_counts)
        
        # Check original group1, group2 values
        message("=== ORIGINAL GROUP VALUES DEBUGGING ===")
        message("Unique group1 values: ", paste(unique(shared_events$group1), collapse = ", "))
        message("Unique group2 values: ", paste(unique(shared_events$group2), collapse = ", "))
        message("Group comparison results:")
        group_comparison <- shared_events %>%
          select(sample1, sample2, group1, group2, debug_group_equal) %>%
          head(10)
        print(group_comparison)
        
        # Debugging: check shared_events data structure
        message("=== SHARED EVENTS STRUCTURE DEBUGGING ===")
        message("Sample shared_events data:")
        print(head(shared_events %>% select(sample1, sample2, group1, group2, taxon, sharing_type, metadata_combination), 10))
        message("Unique sharing_types: ", paste(unique(shared_events$sharing_type), collapse = ", "))
        message("Unique metadata_combinations: ", paste(unique(shared_events$metadata_combination), collapse = ", "))
        
        # Special check for Control-Control cases
        control_control_cases <- shared_events %>% 
          filter(group1 == "Control" & group2 == "Control")
        message("Control-Control cases: ", nrow(control_control_cases))
        if(nrow(control_control_cases) > 0) {
          message("Control-Control sharing_types: ", paste(unique(control_control_cases$sharing_type), collapse = ", "))
          message("Control-Control metadata_combinations: ", paste(unique(control_control_cases$metadata_combination), collapse = ", "))
        }
        
        # Check CRC-Control cases
        crc_control_cases <- shared_events %>% 
          filter((group1 == "CRC" & group2 == "Control") | (group1 == "Control" & group2 == "CRC"))
        message("CRC-Control cases: ", nrow(crc_control_cases))
        if(nrow(crc_control_cases) > 0) {
          message("CRC-Control sharing_types: ", paste(unique(crc_control_cases$sharing_type), collapse = ", "))
          message("CRC-Control metadata_combinations: ", paste(unique(crc_control_cases$metadata_combination), collapse = ", "))
        }
        
        # Intra vs Inter group analysis
        intra_events <- shared_events %>% filter(group1 == group2)
        inter_events <- shared_events %>% filter(group1 != group2)
        
        message("Intra-group events: ", nrow(intra_events))
        message("Inter-group events: ", nrow(inter_events))
        
        if(nrow(intra_events) > 0) {
          message("Intra-group breakdown:")
          print(table(intra_events$sharing_type))
        }
        
        if(nrow(inter_events) > 0) {
          message("Inter-group breakdown:")
          print(table(inter_events$sharing_type))
        }
        
        # Calculate Event-based prevalence
        # 1. Calculate total comparison events (by metadata combination)
        total_events_by_combination <- shared_events %>%
          group_by(metadata_combination, sharing_type) %>%
          summarise(
            total_events = n(),
            unique_sample_pairs = n_distinct(paste0(sample1, "-", sample2)),
            .groups = "drop"
          )
        
        message("Total events by combination:")
        print(total_events_by_combination)
        
        # 2. Calculate Event-based prevalence - by metadata combination
        # Group by taxon, but calculate species belonging to each taxon individually
        event_prevalence_data_raw <- shared_events %>%
          group_by(metadata_combination, sharing_type, taxon) %>%
          summarise(
            shared_events_count = n(),  # Number of events where this taxon was shared
            unique_species = n_distinct(species, na.rm = TRUE),  # Number of unique species in this taxon
            unique_sample_pairs = n_distinct(paste0(sample1, "-", sample2)),  # Number of related sample pairs
            species_list = paste(unique(species), collapse = ", "),  # List of species in this taxon
            .groups = "drop"
          ) %>%
          left_join(total_events_by_combination, by = c("metadata_combination", "sharing_type")) %>%
          mutate(
            event_prevalence = ifelse(total_events > 0, (shared_events_count / total_events) * 100, 0)
          )
        
        # Normalize sharing_type from metadata_combination to avoid mislabeling
        # e.g., ensure "Control ↔ Control" is always Intra-group
        if (nrow(event_prevalence_data_raw) > 0) {
          parts <- strsplit(event_prevalence_data_raw$metadata_combination, "↔", fixed = TRUE)
          left_vals <- trimws(vapply(parts, function(x) if(length(x)>=1) x[1] else NA_character_, character(1)))
          right_vals <- trimws(vapply(parts, function(x) if(length(x)>=2) x[2] else NA_character_, character(1)))
          event_prevalence_data_raw$sharing_type <- ifelse(!is.na(left_vals) & left_vals == right_vals, "Intra-group", "Inter-group")
        }
        
        message("All event prevalence data (before filtering):")
        print(event_prevalence_data_raw %>% select(sharing_type, taxon, event_prevalence) %>% arrange(sharing_type, desc(event_prevalence)))
        
        # Do not apply prevalence threshold in Event-based analysis
        event_prevalence_data <- event_prevalence_data_raw
        
        message("Total event prevalence data rows: ", nrow(event_prevalence_data))
        message("Event occurrence range: ", round(min(event_prevalence_data$event_prevalence, na.rm = TRUE), 2), "% - ", 
                round(max(event_prevalence_data$event_prevalence, na.rm = TRUE), 2), "%")
        
        # Debugging: check event_prevalence_data structure
        message("=== EVENT PREVALENCE DATA STRUCTURE DEBUGGING ===")
        message("Sample event_prevalence_data:")
        print(head(event_prevalence_data %>% select(taxon, sharing_type, metadata_combination, event_prevalence), 10))
        message("Unique sharing_types in event_prevalence_data: ", paste(unique(event_prevalence_data$sharing_type), collapse = ", "))
        message("Unique metadata_combinations in event_prevalence_data: ", paste(unique(event_prevalence_data$metadata_combination), collapse = ", "))
        
        message("Event occurrence data rows: ", nrow(event_prevalence_data))
        
        # Handle case when no data is available
        if(nrow(event_prevalence_data) == 0) {
          showNotification("No shared events found. Try lowering the ANI threshold or check your data.",
                          type = "warning", duration = 5)
          return(list(
            event_data = data.frame(),
            event_prevalence_data = data.frame(),
            total_events_by_combination = data.frame(),
            shared_events = shared_events,
            filtered_data = filtered_data,
            tax_level = tax_col,
            meta_var = input$metadata_var,
            ani_type = ani_col,
            var_type = current_var_type(),
            grouping_mode = input$grouping_mode
          ))
        }
        
        showNotification("Analysis completed successfully!", type = "message", duration = 3)
        
        return(list(
          event_data = shared_events,
          event_prevalence_data = event_prevalence_data,
          total_events_by_combination = total_events_by_combination,
          shared_events = shared_events,
          filtered_data = filtered_data,
          raw_display = raw_display,
          tax_level = tax_col,
          summary_stats = list(
            total_taxa = n_distinct(event_prevalence_data$taxon),
            total_combinations = n_distinct(event_prevalence_data$metadata_combination),
            total_events = nrow(shared_events),
            avg_event_prevalence = round(mean(event_prevalence_data$event_prevalence, na.rm = TRUE), 2)
          ),
          meta_var = input$metadata_var,
          ani_type = ani_col,
          var_type = current_var_type(),
          grouping_mode = input$grouping_mode
        ))
        
      }, error = function(e) {
        global_error_handler(e, session, output, 
          code_snippet = "processed_data reactive function",
          input_list = list(
            metadata_var = input$metadata_var,
            tax_level = input$tax_level,
            ani_type = input$ani_type,
            ani_threshold = input$ani_threshold
          )
        )
        return(NULL)
      })
    })

    # Statistical comparison for categorical variables
    statistical_comparison <- reactive({
      data <- processed_data()
      req(data)

      # Only compute for categorical variables
      if(data$var_type != "categorical") {
        return(NULL)
      }

      tryCatch({
        message("Computing statistical comparison for categorical variable")

        # Get grouping mode
        grouping_mode <- if(!is.null(input$grouping_mode)) input$grouping_mode else "inter_intra"
        message("Statistical comparison using grouping mode: ", grouping_mode)

        # Prepare data from shared_events - count by taxon and sharing_type
        # Use the sharing_type that was already calculated in processed_data based on grouping mode
        counts_by_taxon <- data$shared_events %>%
          group_by(sharing_type, taxon) %>%
          summarise(
            shared_events_count = n(),
            unique_species = n_distinct(species, na.rm = TRUE),
            unique_sample_pairs = n_distinct(paste0(sample1, "-", sample2)),
            .groups = "drop"
          )

        message("Unique sharing types in data: ", paste(unique(counts_by_taxon$sharing_type), collapse = ", "))

        # Count number of unique groups
        n_groups <- n_distinct(counts_by_taxon$sharing_type)
        message("Number of groups for statistical testing: ", n_groups)

        # Adaptive statistical testing based on number of groups
        if(n_groups < 2) {
          message("Insufficient groups for statistical comparison")
          return(list(
            has_comparison = FALSE,
            message = "Insufficient groups for statistical comparison (need at least 2 groups)"
          ))
        } else if(n_groups == 2) {
          # Two groups: Use Wilcoxon test
          message("Using Wilcoxon rank sum test (2 groups)")

          group_names <- unique(counts_by_taxon$sharing_type)
          group1_data <- counts_by_taxon %>% filter(sharing_type == group_names[1])
          group2_data <- counts_by_taxon %>% filter(sharing_type == group_names[2])

          message("Group 1 (", group_names[1], ") data points: ", nrow(group1_data))
          message("Group 2 (", group_names[2], ") data points: ", nrow(group2_data))

          # Check if we have enough data
          if(nrow(group1_data) < 2 || nrow(group2_data) < 2) {
            message("Insufficient data for statistical comparison")
            return(list(
              has_comparison = FALSE,
              message = "Insufficient data for statistical comparison (need at least 2 data points in each group)"
            ))
          }

          # Wilcoxon rank sum test
          test_result <- wilcox.test(
            group1_data$shared_events_count,
            group2_data$shared_events_count,
            alternative = "two.sided"
          )

          # Calculate effect size (rank-biserial correlation)
          n1 <- nrow(group1_data)
          n2 <- nrow(group2_data)
          U <- test_result$statistic
          effect_size <- 1 - (2 * U) / (n1 * n2)

          # Calculate descriptive statistics
          group_stats <- list()
          for(gname in group_names) {
            gdata <- counts_by_taxon %>% filter(sharing_type == gname)
            group_stats[[gname]] <- list(
              n = nrow(gdata),
              median = median(gdata$shared_events_count, na.rm = TRUE),
              mean = mean(gdata$shared_events_count, na.rm = TRUE),
              sd = sd(gdata$shared_events_count, na.rm = TRUE),
              min = min(gdata$shared_events_count, na.rm = TRUE),
              max = max(gdata$shared_events_count, na.rm = TRUE)
            )
          }

          message("Wilcoxon test completed successfully")
          message("P-value: ", test_result$p.value)
          message("Effect size: ", effect_size)

          return(list(
            has_comparison = TRUE,
            test_type = "wilcoxon",
            test_result = test_result,
            effect_size = effect_size,
            group_stats = group_stats,
            n_groups = n_groups,
            data_for_plot = counts_by_taxon
          ))

        } else {
          # Three or more groups: Use Kruskal-Wallis test
          message("Using Kruskal-Wallis test (", n_groups, " groups)")

          # Check data availability for each group
          for(gname in unique(counts_by_taxon$sharing_type)) {
            gdata <- counts_by_taxon %>% filter(sharing_type == gname)
            message("Group '", gname, "' data points: ", nrow(gdata))
            if(nrow(gdata) < 2) {
              return(list(
                has_comparison = FALSE,
                message = paste0("Insufficient data in group '", gname, "' (need at least 2 data points per group)")
              ))
            }
          }

          # Kruskal-Wallis test
          test_result <- kruskal.test(shared_events_count ~ sharing_type, data = counts_by_taxon)

          # Post-hoc pairwise comparisons using Wilcoxon with BH correction
          pairwise_result <- pairwise.wilcox.test(
            counts_by_taxon$shared_events_count,
            counts_by_taxon$sharing_type,
            p.adjust.method = "BH"
          )

          # Calculate descriptive statistics for each group
          group_stats <- list()
          for(gname in unique(counts_by_taxon$sharing_type)) {
            gdata <- counts_by_taxon %>% filter(sharing_type == gname)
            group_stats[[gname]] <- list(
              n = nrow(gdata),
              median = median(gdata$shared_events_count, na.rm = TRUE),
              mean = mean(gdata$shared_events_count, na.rm = TRUE),
              sd = sd(gdata$shared_events_count, na.rm = TRUE),
              min = min(gdata$shared_events_count, na.rm = TRUE),
              max = max(gdata$shared_events_count, na.rm = TRUE)
            )
          }

          # Calculate effect size (eta-squared)
          # eta^2 = (H - k + 1) / (n - k) where H is the test statistic, k is number of groups, n is total sample size
          H <- test_result$statistic
          k <- n_groups
          n_total <- nrow(counts_by_taxon)
          eta_squared <- (H - k + 1) / (n_total - k)

          message("Kruskal-Wallis test completed successfully")
          message("P-value: ", test_result$p.value)
          message("Eta-squared: ", eta_squared)

          return(list(
            has_comparison = TRUE,
            test_type = "kruskal",
            test_result = test_result,
            pairwise_result = pairwise_result,
            effect_size = eta_squared,
            group_stats = group_stats,
            n_groups = n_groups,
            data_for_plot = counts_by_taxon
          ))
        }

      }, error = function(e) {
        message("Error in statistical comparison: ", e$message)
        return(list(
          has_comparison = FALSE,
          message = paste("Error:", e$message)
        ))
      })
    })

    # Statistical summary output
    output$statistical_summary <- renderUI({
      stats <- statistical_comparison()
      req(stats)

      if(!stats$has_comparison) {
        return(tags$div(
          style = "padding: 10px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
          tags$strong("Note: "),
          stats$message %||% "Statistical comparison not available for this data."
        ))
      }

      # Format p-value
      p_val <- stats$test_result$p.value
      p_text <- if(p_val < 0.001) {
        "< 0.001"
      } else {
        sprintf("= %.4f", p_val)
      }

      # Determine significance
      significance <- if(p_val < 0.001) {
        "***"
      } else if(p_val < 0.01) {
        "**"
      } else if(p_val < 0.05) {
        "*"
      } else {
        "ns"
      }

      # Effect size interpretation
      effect_interp <- if(abs(stats$effect_size) < 0.1) {
        "negligible"
      } else if(abs(stats$effect_size) < 0.3) {
        "small"
      } else if(abs(stats$effect_size) < 0.5) {
        "medium"
      } else {
        "large"
      }

      # Background color based on significance
      bg_color <- if(p_val < 0.05) "#d4edda" else "#f8d7da"
      border_color <- if(p_val < 0.05) "#28a745" else "#dc3545"

      # Build UI based on test type
      if(stats$test_type == "wilcoxon") {
        # Wilcoxon test (2 groups)
        test_title <- sprintf("Wilcoxon Rank Sum Test (%d groups)", stats$n_groups)
        test_stat_label <- "Test Statistic (W)"
        effect_label <- "Effect Size (r)"

        # Get group names
        group_names <- names(stats$group_stats)

        tags$div(
          style = sprintf("padding: 15px; background-color: %s; border: 2px solid %s; border-radius: 5px;", bg_color, border_color),
          tags$h5(style = "margin-top: 0;", test_title),
          tags$p(
            style = "margin-bottom: 5px;",
            tags$strong(test_stat_label, ": "), sprintf("%.2f", stats$test_result$statistic)
          ),
          tags$p(
            style = "margin-bottom: 5px;",
            tags$strong("P-value: "), p_text, " ", tags$span(style = "font-size: 1.2em;", significance)
          ),
          tags$p(
            style = "margin-bottom: 5px;",
            tags$strong(effect_label, ": "), sprintf("%.3f", stats$effect_size), " (", effect_interp, ")"
          ),
          tags$hr(style = "margin: 10px 0;"),
          tags$div(
            style = "display: flex; justify-content: space-around;",
            lapply(group_names, function(gname) {
              gstats <- stats$group_stats[[gname]]
              tags$div(
                style = "text-align: center;",
                tags$strong(gname),
                tags$br(),
                sprintf("n = %d", gstats$n),
                tags$br(),
                sprintf("Median = %.1f", gstats$median),
                tags$br(),
                sprintf("Mean ± SD = %.1f ± %.1f", gstats$mean, gstats$sd)
              )
            })
          ),
          tags$p(
            style = "margin-top: 10px; margin-bottom: 0; font-size: 0.9em; color: #666;",
            "*** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant"
          )
        )

      } else {
        # Kruskal-Wallis test (3+ groups)
        test_title <- sprintf("Kruskal-Wallis Test (%d groups)", stats$n_groups)
        test_stat_label <- "Test Statistic (H)"
        effect_label <- "Effect Size (η²)"

        # Get group names
        group_names <- names(stats$group_stats)

        # Build pairwise comparison text
        pairwise_text <- NULL
        if(!is.null(stats$pairwise_result)) {
          pairwise_matrix <- stats$pairwise_result$p.value
          if(!is.null(pairwise_matrix)) {
            sig_pairs <- c()
            for(i in 2:nrow(pairwise_matrix)) {
              for(j in 1:(i-1)) {
                p_adj <- pairwise_matrix[i, j]
                if(!is.na(p_adj) && p_adj < 0.05) {
                  sig <- if(p_adj < 0.001) "***" else if(p_adj < 0.01) "**" else "*"
                  sig_pairs <- c(sig_pairs, sprintf("%s vs %s: p = %.4f %s",
                                                    rownames(pairwise_matrix)[i],
                                                    colnames(pairwise_matrix)[j],
                                                    p_adj, sig))
                }
              }
            }
            if(length(sig_pairs) > 0) {
              pairwise_text <- paste(sig_pairs, collapse = "; ")
            } else {
              pairwise_text <- "No significant pairwise differences (BH-corrected)"
            }
          }
        }

        tags$div(
          style = sprintf("padding: 15px; background-color: %s; border: 2px solid %s; border-radius: 5px;", bg_color, border_color),
          tags$h5(style = "margin-top: 0;", test_title),
          tags$p(
            style = "margin-bottom: 5px;",
            tags$strong(test_stat_label, ": "), sprintf("%.2f", stats$test_result$statistic)
          ),
          tags$p(
            style = "margin-bottom: 5px;",
            tags$strong("P-value: "), p_text, " ", tags$span(style = "font-size: 1.2em;", significance)
          ),
          tags$p(
            style = "margin-bottom: 5px;",
            tags$strong(effect_label, ": "), sprintf("%.3f", stats$effect_size), " (", effect_interp, ")"
          ),
          if(!is.null(pairwise_text)) {
            tags$p(
              style = "margin-bottom: 5px;",
              tags$strong("Post-hoc Pairwise Comparisons (BH-adjusted): "),
              tags$br(),
              tags$span(style = "font-size: 0.9em;", pairwise_text)
            )
          },
          tags$hr(style = "margin: 10px 0;"),
          tags$div(
            style = "display: flex; flex-wrap: wrap; justify-content: space-around;",
            lapply(group_names, function(gname) {
              gstats <- stats$group_stats[[gname]]
              tags$div(
                style = "text-align: center; margin: 5px;",
                tags$strong(gname),
                tags$br(),
                sprintf("n = %d", gstats$n),
                tags$br(),
                sprintf("Median = %.1f", gstats$median),
                tags$br(),
                sprintf("Mean ± SD = %.1f ± %.1f", gstats$mean, gstats$sd)
              )
            })
          ),
          tags$p(
            style = "margin-top: 10px; margin-bottom: 0; font-size: 0.9em; color: #666;",
            "*** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant"
          )
        )
      }
    })

    # Statistical results table
    output$statistical_results_table <- DT::renderDataTable({
      stats <- statistical_comparison()
      req(stats)

      if(!stats$has_comparison) {
        return(DT::datatable(data.frame(
          Message = "No statistical comparison available"
        )))
      }

      # Create detailed comparison table by sharing type
      comparison_table <- stats$data_for_plot %>%
        group_by(sharing_type) %>%
        summarise(
          N = n(),
          `Total Events` = sum(shared_events_count),
          `Mean Events` = mean(shared_events_count),
          `Median Events` = median(shared_events_count),
          `SD` = sd(shared_events_count),
          `Min` = min(shared_events_count),
          `Max` = max(shared_events_count),
          .groups = "drop"
        ) %>%
        arrange(sharing_type)

      DT::datatable(
        comparison_table,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        extensions = 'Buttons',
        caption = "Descriptive Statistics by Metadata Combination and Sharing Type",
        rownames = FALSE
      ) %>%
        DT::formatRound(columns = c("Mean Events", "Median Events", "SD"), digits = 2)
    })

    # Handle Group Comparison Plot point click event
    observeEvent(event_data("plotly_click", source = "group_comparison"), {
      tryCatch({
        message("=== CLICK EVENT START ===")
        click_data <- event_data("plotly_click", source = "group_comparison")
        message("Click data: ", toString(click_data))
        # Safe default values: initialize first so they can be referenced in any branch
        detailed_data <- data.frame()
        clicked_point_data <- data.frame()
        
        if(!is.null(click_data)) {
          data <- processed_data()
          req(data)
          message("Data loaded successfully for click event")
          
          message("=== PROCESSING CLICK DATA ===")
          # Extract clicked point information
          clicked_x <- click_data$x
          clicked_y <- click_data$y
          message("Clicked coordinates: x=", clicked_x, ", y=", clicked_y)
          message("Variable type: ", data$var_type)
          
          # Check full shared_events data
          message("=== DEBUGGING: Data Source Comparison ===")
          message("Total shared_events in processed_data: ", nrow(data$shared_events))
          message("Total event_prevalence_data in processed_data: ", nrow(data$event_prevalence_data))
          
          message("Sample shared_events entries:")
          if(nrow(data$shared_events) > 0) {
            sample_events <- head(data$shared_events, 3)
            for(i in 1:nrow(sample_events)) {
              message("  Event ", i, ": ", sample_events$sample1[i], " vs ", sample_events$sample2[i], 
                     " | Species: ", sample_events$species[i], " | Taxon: ", sample_events$taxon[i])
            }
          }
          
          message("Sample event_prevalence_data entries:")
          if(nrow(data$event_prevalence_data) > 0) {
            sample_prevalence <- head(data$event_prevalence_data, 3)
            for(i in 1:nrow(sample_prevalence)) {
              message("  Prevalence ", i, ": Taxon: ", sample_prevalence$taxon[i], 
                     " | Sharing: ", sample_prevalence$sharing_type[i], 
                     " | Prevalence: ", sample_prevalence$event_prevalence[i])
            }
          }
          
          # Compare taxon lists
          shared_taxons <- unique(data$shared_events$taxon)
          prevalence_taxons <- unique(data$event_prevalence_data$taxon)
          message("Unique taxons in shared_events: ", length(shared_taxons))
          message("Unique taxons in event_prevalence_data: ", length(prevalence_taxons))
          message("Taxons only in shared_events: ", paste(setdiff(shared_taxons, prevalence_taxons), collapse = ", "))
          message("Taxons only in event_prevalence_data: ", paste(setdiff(prevalence_taxons, shared_taxons), collapse = ", "))
          
          if(data$var_type == "categorical") {
            # Handle categorical variable click - restructure to match count-based plot
            message("Processing categorical variable click")

            # 1) Regenerate aggregated data same as plot (count-based)
            counts_by_taxon <- data$shared_events %>%
              group_by(metadata_combination, sharing_type, taxon) %>%
              summarise(
                shared_events_count = n(),
                unique_species = n_distinct(species, na.rm = TRUE),
                unique_sample_pairs = n_distinct(paste0(sample1, "-", sample2)),
                .groups = "drop"
              )
            totals_by_comb <- data$shared_events %>%
              group_by(metadata_combination, sharing_type) %>%
              summarise(total_events = n(), .groups = "drop")
            plot_data_counts <- counts_by_taxon %>%
              left_join(totals_by_comb, by = c("metadata_combination", "sharing_type"))

            # 2) Extract taxon/combination/type from click text or key
            clicked_text <- if(!is.null(click_data$text)) as.character(click_data$text) else ""
            clicked_key  <- if(!is.null(click_data$key)) as.character(click_data$key) else ""
            parse_field <- function(txt, label) {
              m <- regmatches(txt, regexpr(paste0(label, " [^<]+"), txt))
              if(length(m) == 1) sub(paste0(label, " "), "", m) else NA_character_
            }
            # Try tooltip first
            clicked_taxon  <- parse_field(clicked_text, "Taxon:")
            clicked_combo  <- parse_field(clicked_text, "Combination:")
            clicked_type   <- parse_field(clicked_text, "Type:")
            # Fallback to key: taxon||combination||type
            if((is.na(clicked_taxon) || clicked_taxon == "") && grepl("\\|\\|", clicked_key)) {
              parts <- strsplit(clicked_key, "\\|\\|", fixed = FALSE)[[1]]
              if(length(parts) >= 3) {
                clicked_taxon <- parts[1]
                clicked_combo <- parts[2]
                clicked_type  <- parts[3]
              }
            }

            # Fallback to pointNumber index if still missing
            if((is.na(clicked_taxon) || clicked_taxon == "" || is.na(clicked_combo) || clicked_combo == "" || is.na(clicked_type) || clicked_type == "") &&
               !is.null(click_data$pointNumber) && NROW(plot_data_counts) > 0) {
              idx <- as.integer(click_data$pointNumber) + 1L
              if(!is.na(idx) && idx >= 1 && idx <= NROW(plot_data_counts)) {
                fallback_row <- plot_data_counts[idx, ]
                clicked_taxon <- ifelse(is.na(clicked_taxon) || clicked_taxon == "", as.character(fallback_row$taxon), clicked_taxon)
                clicked_combo <- ifelse(is.na(clicked_combo) || clicked_combo == "", as.character(fallback_row$metadata_combination), clicked_combo)
                clicked_type  <- ifelse(is.na(clicked_type)  || clicked_type  == "", as.character(fallback_row$sharing_type), clicked_type)
              }
            }

            # 3) Correction: if text is empty, correct from coordinates
            if(is.na(clicked_type) || clicked_type == "") {
              clicked_type <- if(is.character(clicked_x)) clicked_x else as.character(clicked_x)
              # numeric fallback -> map to category index
              if(!clicked_type %in% c("Inter-group", "Intra-group")) {
                levels <- c("Inter-group", "Intra-group")
                idx <- round(as.numeric(clicked_type))
                if(is.finite(idx) && idx >= 1 && idx <= length(levels)) clicked_type <- levels[idx] else clicked_type <- levels[1]
              }
            }

            message("=== CATEGORICAL CLICK DATA ===")
            message("Target taxon: ", clicked_taxon)
            message("Target sharing_type: ", clicked_type)
            message("Target metadata_combination: ", clicked_combo)

            # 4) Detailed data filter: Use pre-computed sharing_type from processed_data
            # This respects the grouping mode (inter_intra or original) selected by the user
            shared_events_data <- data$shared_events %>%
              filter(
                taxon == clicked_taxon,
                # Use the pre-existing sharing_type and metadata_combination columns
                # These were computed in processed_data() based on the selected grouping mode
                (is.na(clicked_type)  | clicked_type  == "" | sharing_type == clicked_type),
                (is.na(clicked_combo) | clicked_combo == "" | metadata_combination == clicked_combo)
              ) %>%
              distinct(sample1, sample2, species, .keep_all = TRUE)

            message("Filtered shared_events data (after calc-based filters): ", nrow(shared_events_data), " unique events")
            if(nrow(shared_events_data) == 0) {
              showNotification("No events for the selected point", type = "warning")
              clicked_data_rv$data <- data.frame()
              return()
            }

           # Defer building of detailed_data to the unified block below where
           # coverage_summary is computed. This avoids referencing variables
           # before they are created and keeps categorical/numerical flows consistent.
            
          } else if(data$var_type == "numerical") {
            # Numerical variable processing - use direct method
            message("Processing numerical variable click")
            
            # Use actual data information from clicked point in Plotly
            # click_data may contain pointNumber or curveNumber information
            message("Click data details: ", paste(names(click_data), "=", click_data, collapse = ", "))
            
            # Find data that exactly matches clicked coordinates
            shared_events <- data$shared_events
            
            # Regenerate summary_data (same structure as plot)
            # Load GTDB data and map phylum (same logic as plot)
            filtered_gtdb_path <- "data/filtered_gtdb_metadata.rds"
            message("Looking for GTDB file at: ", filtered_gtdb_path)
            if(file.exists(filtered_gtdb_path)) {
              message("GTDB file found, loading...")
              gtdb_data <- readRDS(filtered_gtdb_path)
              message("GTDB data loaded: ", nrow(gtdb_data), " records")
              
              current_tax_level <- input$tax_level
              message("Current tax level: ", current_tax_level)
              
              event_data <- shared_events %>%
                mutate(
                  x_difference = abs(as.numeric(as.character(group1)) - as.numeric(as.character(group2))),
                  y_value = as.numeric(as.character(group1)),
                  pair_id = paste0(pmin(sample1, sample2), "_", pmax(sample1, sample2))
                ) %>%
                left_join(gtdb_data %>% select(!!sym(current_tax_level), phylum), 
                         by = setNames(current_tax_level, "taxon")) %>%
                mutate(phylum = if_else(is.na(phylum), "Unknown", phylum))
              
              message("Event data after GTDB join: ", nrow(event_data), " rows")
              
              # Generate summary data (same as plot)
              message("Creating summary data...")
              summary_data <- tryCatch({
                event_data %>%
                  group_by(x_difference, sharing_type, sample1, sample2) %>%
                  summarise(
                    unique_pairs = n_distinct(paste(sample1, sample2)),
                    unique_species = n_distinct(species),
                    total_events = n(),
                    dominant_phylum = {
                      phylum_clean <- as.character(phylum[!is.na(phylum) & phylum != ""])
                      if(length(phylum_clean) > 0) {
                        phylum_counts <- table(phylum_clean)
                        as.character(names(phylum_counts)[which.max(phylum_counts)])
                      } else {
                        "Unknown"
                      }
                    },
                    y_value = first(y_value),
                    .groups = 'drop'
                  )
              }, error = function(e) {
                message("Error creating summary data: ", e$message)
                return(data.frame())
              })
              
              message("Summary data created: ", nrow(summary_data), " rows")
            } else {
              message("GTDB file not found, cannot process numerical click")
              return()
            }
            
            # Find point that exactly matches clicked coordinates (based on summary_data)
            tolerance <- 2.0  # increased tolerance
            matched_point <- summary_data %>%
              filter(
                abs(x_difference - clicked_x) <= tolerance,
                abs(y_value - clicked_y) <= tolerance
              ) %>%
              slice(1)
            
            message("Found ", nrow(matched_point), " matching summary points")
            
            if(nrow(matched_point) > 0) {
              # Extract sample pair information from matched point
              target_sample1 <- matched_point$sample1
              target_sample2 <- matched_point$sample2
              target_sharing_type <- matched_point$sharing_type
              
              message("=== SUMMARY MATCH FOUND ===")
              message("Sample Pair: ", target_sample1, " vs ", target_sample2)
              
              # Find all shared strains from this sample pair
              all_shared_species <- shared_events %>%
                filter(
                  sample1 == target_sample1,
                  sample2 == target_sample2,
                  sharing_type == target_sharing_type
                )
              
              message("Found ", nrow(all_shared_species), " shared species for this pair")
              
              # Use first species as representative, show all species later
              clicked_point_data <- all_shared_species %>% slice(1)
              target_sample_pair <- paste(target_sample1, target_sample2, sep = " vs ")
              
            } else {
              message("=== NO MATCH FOUND ===")
              message("Available summary points: ", nrow(summary_data))
              if(nrow(summary_data) > 0) {
                message("X range: ", round(min(summary_data$x_difference), 2), " - ", round(max(summary_data$x_difference), 2))
                message("Y range: ", round(min(summary_data$y_value), 2), " - ", round(max(summary_data$y_value), 2))
                message("Clicked: X=", round(clicked_x, 2), ", Y=", round(clicked_y, 2))
              }
              return()
              
            }
          } else {
            clicked_point_data <- data.frame()
          }
          
          # Display only data related to clicked point
          if(data$var_type == "numerical" && exists("all_shared_species") && nrow(all_shared_species) > 0) {
              # Numerical variable: show all shared strains for this sample pair
              message("=== SHOWING ALL SHARED SPECIES FOR PAIR ===")
              message("Sample Pair: ", target_sample1, " vs ", target_sample2)
              message("Total shared species: ", nrow(all_shared_species))
              
              # Load genome info for all shared strains
              shared_events_data <- all_shared_species %>%
                distinct(sample1, sample2, species, .keep_all = TRUE)
              
              message("Processing ", nrow(shared_events_data), " shared species for genome info")
          } else if(data$var_type == "categorical") {
              # Categorical variable: Use already-filtered shared_events_data from the earlier block
              # The first categorical block (lines 1568-1657) already parsed and filtered the data
              # We should NOT filter again here - just use what was already created

              message("=== CATEGORICAL: Using pre-filtered data ===")

              # Check if shared_events_data was created in the earlier categorical block
              if(!exists("shared_events_data") || is.null(shared_events_data) || nrow(shared_events_data) == 0) {
                message("WARNING: shared_events_data not found from earlier block, filtering now as fallback")

                target_taxon <- clicked_taxon
                target_sharing_type <- clicked_type
                target_metadata_combination <- clicked_combo

                message("=== CATEGORICAL CLICK DATA (FALLBACK) ===")
                message("Target taxon: ", target_taxon)
                message("Target sharing_type: ", target_sharing_type)
                message("Target metadata_combination: ", target_metadata_combination)

                # Use pre-computed sharing_type from processed_data (fallback path)
                # This respects the grouping mode selected by the user
                shared_events_data <- data$shared_events %>%
                  filter(
                    taxon == target_taxon,
                    # IMPORTANT: Use pre-existing sharing_type and metadata_combination
                    # These were computed in processed_data() based on the selected grouping mode
                    (is.na(target_sharing_type)  | target_sharing_type  == "" | sharing_type == target_sharing_type),
                    (is.na(target_metadata_combination) | target_metadata_combination == "" | metadata_combination == target_metadata_combination)
                  ) %>%
                  distinct(sample1, sample2, species, .keep_all = TRUE)

                message("Filtered shared_events data (fallback): ", nrow(shared_events_data), " unique events")
              } else {
                message("Using shared_events_data from earlier block: ", nrow(shared_events_data), " unique events")
              }

              if(nrow(shared_events_data) > 0) {
                message("Filtered shared_events sample:")
                print(head(shared_events_data %>% select(sample1, sample2, taxon, sharing_type, metadata_combination, species), 5))
              }

              # Extract related samples
              relevant_samples <- unique(c(shared_events_data$sample1, shared_events_data$sample2))
              message("Relevant samples: ", paste(relevant_samples, collapse = ", "))
              # Downstream sections expect clicked_point_data to exist
              clicked_point_data <- if(nrow(shared_events_data) > 0) shared_events_data %>% slice(1) else data.frame()
            }
            
            # Try to load genome_info data
            tryCatch({
              instrain_base_path <- Sys.getenv("INSTRAIN_DATA_PATH", "/data/instrain/")

              fd_raw <- data$filtered_data %>%
                mutate(join_key = paste(clean_name1, clean_name2, species, sep = "|"))

              coverage_list <- lapply(seq_len(nrow(shared_events_data)), function(i) {
                row <- shared_events_data[i, ]
                join_key <- paste(row$sample1, row$sample2, row$species, sep = "|")

                sample1_genome_data <- load_genome_info_for_samples(instrain_base_path, row$sample1, row$species)
                sample2_genome_data <- load_genome_info_for_samples(instrain_base_path, row$sample2, row$species)

                sample1_info <- if(nrow(sample1_genome_data) > 0) {
                  sample1_genome_data %>%
                    slice_max(`Coverage`, n = 1) %>%
                    transmute(
                      `Full Taxonomy_Sample1` = `Full Taxonomy`,
                      `Coverage_Sample1` = Coverage,
                      `Breadth_Sample1` = Breadth,
                      `Nucleotide Diversity_Sample1` = `Nucleotide Diversity`,
                      `SNV Count_Sample1` = `SNV Count`,
                      `Genome Size_Sample1` = `Genome Size`
                    )
                } else {
                  tibble::tibble(
                    `Full Taxonomy_Sample1` = NA_character_,
                    `Coverage_Sample1` = NA_real_,
                    `Breadth_Sample1` = NA_real_,
                    `Nucleotide Diversity_Sample1` = NA_real_,
                    `SNV Count_Sample1` = NA_real_,
                    `Genome Size_Sample1` = NA_real_
                  )
                }

                sample2_info <- if(nrow(sample2_genome_data) > 0) {
                  sample2_genome_data %>%
                    slice_max(`Coverage`, n = 1) %>%
                    transmute(
                      `Full Taxonomy_Sample2` = `Full Taxonomy`,
                      `Coverage_Sample2` = Coverage,
                      `Breadth_Sample2` = Breadth,
                      `Nucleotide Diversity_Sample2` = `Nucleotide Diversity`,
                      `SNV Count_Sample2` = `SNV Count`,
                      `Genome Size_Sample2` = `Genome Size`
                    )
                } else {
                  tibble::tibble(
                    `Full Taxonomy_Sample2` = NA_character_,
                    `Coverage_Sample2` = NA_real_,
                    `Breadth_Sample2` = NA_real_,
                    `Nucleotide Diversity_Sample2` = NA_real_,
                    `SNV Count_Sample2` = NA_real_,
                    `Genome Size_Sample2` = NA_real_
                  )
                }

                tibble::tibble(join_key) %>%
                  bind_cols(sample1_info) %>%
                  bind_cols(sample2_info)
              })

              coverage_summary <- if(length(coverage_list) > 0) dplyr::bind_rows(coverage_list) else tibble::tibble(join_key = character(0))
              message("Coverage summary rows: ", nrow(coverage_summary))

              selected_keys <- paste(shared_events_data$sample1, shared_events_data$sample2, shared_events_data$species, sep = "|")

              detailed_data <- data$raw_display %>%
                mutate(join_key = paste(clean_name1, clean_name2, species, sep = "|")) %>%
                filter(join_key %in% selected_keys) %>%
                left_join(coverage_summary, by = "join_key") %>%
                # Explicitly select only needed columns with select() (removed everything())
                select(
                  `Sample 1` = clean_name1,
                  `Sample 2` = clean_name2,
                  `Taxon` = !!sym(data$tax_level),
                  `Species` = species,
                  `Group 1` = group1,
                  `Group 2` = group2,
                  `Sharing Type` = sharing_type,
                  `Combination` = metadata_combination,
                  popANI, conANI,
                  starts_with("Full"), ends_with("Sample1"), ends_with("Sample2"),
                  contains("Coverage"), contains("Breadth"),
                  contains("Nucleotide"), contains("SNV"), contains("Genome Size")
                )

              tax_col <- data$tax_level
              preferred_front <- c(
                "Sample 1","Sample 2","Species",tax_col,"popANI","conANI",
                "Group 1","Group 2","Sharing Type","Combination",
                "Full Taxonomy","Full Taxonomy_Sample1","Coverage_Sample1","Breadth_Sample1",
                "Nucleotide Diversity_Sample1","SNV Count_Sample1","Genome Size_Sample1",
                "Full Taxonomy_Sample2","Coverage_Sample2","Breadth_Sample2",
                "Nucleotide Diversity_Sample2","SNV Count_Sample2","Genome Size_Sample2"
              )
              existing_front <- intersect(preferred_front, colnames(detailed_data))
              other_cols <- setdiff(colnames(detailed_data), existing_front)
              detailed_data <- detailed_data[, c(existing_front, other_cols), drop = FALSE]
              message("Showing combined data with genome_info: ", nrow(detailed_data), " records")
              
            }, error = function(e) {
              message("Error loading genome_info: ", e$message)
              # Even on error, join with original filtered_data to provide as many columns as possible
              fd <- data$filtered_data
              detailed_data <- shared_events_data %>%
                select(
                  sample1, sample2, taxon, species, group1, group2, sharing_type, metadata_combination
                ) %>%
                left_join(
                  fd,
                  by = c("sample1" = "clean_name1", "sample2" = "clean_name2", "species" = "species")
                ) %>%
                # Change column names for display and organize order
                rename(
                  `Sample 1` = sample1,
                  `Sample 2` = sample2,
                  `Taxon` = taxon,
                  `Species` = species,
                  `Group 1` = group1,
                  `Group 2` = group2,
                  `Sharing Type` = sharing_type,
                  `Combination` = metadata_combination
                )
              # Place main columns at the front
              tax_col <- data$tax_level
              front_cols <- c("Sample 1","Sample 2","Species",tax_col,"popANI","conANI","Group 1","Group 2","Sharing Type","Combination")
              front_cols <- intersect(front_cols, colnames(detailed_data))
              other_cols <- setdiff(colnames(detailed_data), front_cols)
              detailed_data <- detailed_data[, c(front_cols, other_cols), drop = FALSE]
            })
            
          # Store data in reactive value (safe even if empty)
          clicked_data_rv$data <- detailed_data
          # Force re-render: ensure immediate table update even for numeric clicks
          output$clicked_point_table <- DT::renderDataTable({
            req(clicked_data_rv$data)
            DT::datatable(
              clicked_data_rv$data,
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'colvis')
              ),
              extensions = 'Buttons',
              caption = paste0("Shared Events Data (", nrow(clicked_data_rv$data), " events)")
            )
          })
            
          # Show table section
          message("Showing data table section")
          shinyjs::show("clicked_data_section")
          message("Data table section should now be visible")
          
          showNotification(
            paste("Data table updated with", nrow(detailed_data), "events"),
            type = "message",
            duration = 3
          )
          
        }
      }, error = function(e) {
        message("=== CLICK EVENT ERROR ===")
        message("Error details: ", e$message)
        message("Error call: ", deparse(e$call))
        if(exists("traceback")) {
          message("Traceback: ", paste(traceback(), collapse = "\n"))
        }
        showNotification(paste("Click event error:", e$message), type = "error")
      })
    })
    
    # Clicked point data table
    output$clicked_point_table <- DT::renderDataTable({
      req(clicked_data_rv$data)
      DT::datatable(
        clicked_data_rv$data,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'colvis')
        ),
        extensions = 'Buttons',
        caption = paste0("Shared Events Data (", nrow(clicked_data_rv$data), " events)")
      )
    })

    # Download handler for Selected Point Details (CSV)
    output$download_selected_details <- downloadHandler(
      filename = function() {
        paste0("selected_point_details_", Sys.Date(), ".csv")
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
    
    # Value boxes
    output$total_taxa <- renderValueBox({
      data <- processed_data()
      if(is.null(data)) return(valueBox(0, "Total Taxa", icon = icon("bacteria"), color = "blue"))
      
      n_taxa <- length(unique(data$shared_events$taxon))
      valueBox(
        value = n_taxa,
        subtitle = "Total Taxa",
        icon = icon("bacteria"),
        color = "blue"
      )
    })
    
    output$shared_taxa <- renderValueBox({
      data <- processed_data()
      if(is.null(data)) return(valueBox(0, "Shared Taxa", icon = icon("share-alt"), color = "green"))
      
      # Total shared events count
      total_events <- nrow(data$shared_events)
      
      valueBox(
        value = total_events,
        subtitle = "Total Shared Events",
        icon = icon("share-alt"),
        color = "green"
      )
    })
    
    output$total_comparisons <- renderValueBox({
      data <- processed_data()
      if(is.null(data)) return(valueBox(0, "Event Combinations", icon = icon("chart-line"), color = "yellow"))
      
      n_combinations <- nrow(data$total_events_by_combination)
      
      valueBox(
        value = n_combinations,
        subtitle = "Event Combinations",
        icon = icon("chart-line"),
        color = "yellow"
      )
    })
    
    output$filtered_comparisons <- renderValueBox({
      req(integrated_data())
      data <- integrated_data()
      
      valueBox(
        value = nrow(data),
        subtitle = "Total Comparisons",
        icon = icon("database"),
        color = "purple"
      )
    })
    
    # Overview metrics with text outputs (like Available Variables)
    output$overview_total_samples_text <- renderText({
      data <- processed_data()
      if(is.null(data) || is.null(data$shared_events) || nrow(data$shared_events) == 0) {
        return("0")
      }
      all_samples <- unique(c(data$shared_events$sample1, data$shared_events$sample2))
      return(as.character(length(all_samples)))
    })

    output$overview_unique_pairs_text <- renderText({
      data <- processed_data()
      if(is.null(data) || is.null(data$shared_events) || nrow(data$shared_events) == 0) {
        return("0")
      }
      pair_count <- data$shared_events %>%
        mutate(pair_id = paste(sample1, sample2, sep = " ↔ ")) %>%
        summarise(unique_pairs = n_distinct(pair_id)) %>%
        pull(unique_pairs)
      return(as.character(pair_count))
    })

    output$overview_total_events_text <- renderText({
      data <- processed_data()
      if(is.null(data) || is.null(data$shared_events)) {
        return("0")
      }
      return(as.character(nrow(data$shared_events)))
    })

    output$overview_events_per_pair_text <- renderText({
      data <- processed_data()
      if(is.null(data) || is.null(data$shared_events) || nrow(data$shared_events) == 0) {
        return("0")
      }
      pair_stats <- data$shared_events %>%
        mutate(pair_id = paste(sample1, sample2, sep = " ↔ ")) %>%
        count(pair_id, name = "event_count")
      avg_events <- mean(pair_stats$event_count)
      return(as.character(round(avg_events, 2)))
    })

    # Overview plots
    output$overview_events_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      shared_events <- data$shared_events
      if(is.null(shared_events) || nrow(shared_events) == 0) {
        return(plot_ly() %>%
                 add_text(x = 0.5, y = 0.5, text = "No shared events available",
                          textfont = list(size = 16, color = "#6c757d")) %>%
                 layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      # Count events per pair and create histogram
      pair_data <- shared_events %>%
        mutate(pair_id = paste(sample1, sample2, sep = " ↔ ")) %>%
        count(pair_id, name = "event_count")
      
      # Create histogram of event counts
      event_counts <- table(pair_data$event_count)
      hist_data <- data.frame(
        events_in_pair = as.numeric(names(event_counts)),
        frequency = as.numeric(event_counts)
      )

      p <- ggplot(hist_data, aes(x = events_in_pair, y = frequency,
                                 text = paste0("Events in Pair: ", events_in_pair, "<br>Frequency: ", frequency))) +
        geom_col(fill = "#4e79a7", alpha = 0.8) +
        labs(x = "Number of Shared Events in a Pair", 
             y = "Event Occurrence Frequency") +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())

      ggplotly(p, tooltip = "text")
    })

    output$overview_taxonomy_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      shared_events <- data$shared_events
      if(is.null(shared_events) || nrow(shared_events) == 0) {
        return(plot_ly() %>%
                 add_text(x = 0.5, y = 0.5, text = "No taxonomy data available",
                          textfont = list(size = 16, color = "#6c757d")) %>%
                 layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      tax_summary <- shared_events %>%
        filter(!is.na(taxon) & taxon != "") %>%
        count(taxon, name = "event_count") %>%
        arrange(desc(event_count)) %>%
        slice_head(n = 20)

      if(nrow(tax_summary) == 0) {
        return(plot_ly() %>%
                 add_text(x = 0.5, y = 0.5, text = "No taxonomy data available",
                          textfont = list(size = 16, color = "#6c757d")) %>%
                 layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      tax_summary$taxon <- factor(tax_summary$taxon, levels = rev(tax_summary$taxon))

      p <- ggplot(tax_summary, aes(x = taxon, y = event_count,
                                   text = paste0("Taxon: ", taxon, "<br>Events: ", event_count))) +
        geom_col(fill = "#f28e2b") +
        coord_flip() +
        labs(x = "Taxon", y = "Shared Events") +
        theme_minimal(base_size = 13)

      ggplotly(p, tooltip = "text") %>%
        layout(margin = list(l = 140))
    })
    
    # Event-based Group comparison plot
    output$group_comparison_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      
      tryCatch({
        # Check variable type - categorical vs numerical
        var_type <- data$var_type
        
        # Treat Quartile or binned data as categorical
        group1_values <- data$shared_events$group1[!is.na(data$shared_events$group1)]
        is_binned_data <- any(grepl("Q[1-4]|Bin_", group1_values))
        
        # Automatically detect categorical data
        is_categorical_data <- FALSE
        if(length(group1_values) > 0) {
          # Try converting to number
          numeric_test <- tryCatch({
            as.numeric(as.character(group1_values))
          }, error = function(e) {
            return(rep(NA, length(group1_values)))
          })
          
          # If many NA, judge as categorical
          na_ratio <- sum(is.na(numeric_test)) / length(numeric_test)
          is_categorical_data <- na_ratio > 0.5  # Categorical if 50% or more are NA
        }
        
        message("Variable type: ", var_type)
        message("Is binned data: ", is_binned_data)
        message("Is categorical data: ", is_categorical_data)
        message("Sample group1 values: ", paste(head(group1_values), collapse = ", "))
        
        # For categorical variables, visualize distribution by metadata combination
        if(var_type == "categorical" || is_binned_data || is_categorical_data) {
          message("Processing categorical variable")
          # Categorical variable: visualize distribution by metadata combination

          # Use the sharing_type that was already calculated in processed_data based on grouping mode
          # This respects the user's choice between Inter/Intra mode or Original Categories mode
          counts_by_taxon <- data$shared_events %>%
            group_by(sharing_type, taxon) %>%
            summarise(
              shared_events_count = n(),
              unique_species = n_distinct(species, na.rm = TRUE),
              unique_sample_pairs = n_distinct(paste0(sample1, "-", sample2)),
              metadata_combination = first(metadata_combination),
              .groups = "drop"
            )
          totals_by_sharing <- data$shared_events %>%
            group_by(sharing_type) %>%
            summarise(total_events = n(), .groups = "drop")
          plot_data <- counts_by_taxon %>%
            left_join(totals_by_sharing, by = "sharing_type")
          plot_title <- "Shared Event Counts by Grouping"
          
          # Handle case when no data is available
          if(nrow(plot_data) == 0) {
            return(plot_ly() %>% 
              add_text(x = 0.5, y = 0.5, text = "No data available for plotting", 
                       textfont = list(size = 20, color = "red")) %>%
              layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
          }
          
          # Get color palette safely
          tryCatch({
            color_values <- get_palette(input$color_palette)
          }, error = function(e) {
            color_values <- rainbow(length(unique(plot_data$sharing_type)))
          })
          
          # Debugging: check plot_data structure
          message("=== CATEGORICAL PLOT DATA DEBUGGING ===")
          message("Plot data rows: ", nrow(plot_data))
          message("Unique sharing_types: ", paste(unique(plot_data$sharing_type), collapse = ", "))
          message("Unique metadata_combinations: ", paste(unique(plot_data$metadata_combination), collapse = ", "))
          message("Sample plot data (counts):")
          print(head(plot_data %>% select(taxon, sharing_type, metadata_combination, shared_events_count), 10))
          
          # Categorical plot: shared event count by metadata combination (box + jitter)
          p <- ggplot(plot_data, aes(x = sharing_type, y = shared_events_count, fill = sharing_type)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.2) +
            geom_jitter(
              width = 0.2,
              alpha = input$point_alpha,
              size = input$point_size,
              mapping = aes(
                text = paste0(
                  "Taxon: ", taxon, "<br>",
                  "Combination: ", metadata_combination, "<br>",
                  "Type: ", sharing_type, "<br>",
                  "Shared Events: ", shared_events_count
                ),
                key = interaction(taxon, metadata_combination, sharing_type, sep = "||", drop = TRUE)
              )
            ) +
            guides(fill = "none") +
            theme(legend.position = "none") +
            labs(title = plot_title,
                 x = "Sharing Type",
                 y = "Shared Event Count") +
            scale_fill_manual(values = color_values) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        } else {
          # Numerical variable: show points and lines per pair (without cross-group grouping)
          
          # Shared events data for selected taxon
          if(!is.null(input$selected_taxon) && input$selected_taxon != "" && input$selected_taxon != "all") {
            event_data <- data$shared_events %>% filter(taxon == input$selected_taxon)
            plot_title <- paste("Shared Events for", input$selected_taxon, "- Sample Pairs")
          } else {
            event_data <- data$shared_events
            plot_title <- "All Shared Events - Sample Value Pairs"
          }
          
          # Handle case when no data is available
          if(nrow(event_data) == 0) {
            return(plot_ly() %>% 
              add_text(x = 0.5, y = 0.5, text = "No data available for plotting", 
                       textfont = list(size = 20, color = "red")) %>%
              layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
          }
          
          # Numerical variable: directly display original values for each pair
          message("Processing numerical variable - direct pair values")
          
          # Convert each pair to two points (Sample1 value, Sample2 value)
          pair_data <- event_data %>%
            mutate(
              # Sample1 data
              sample_type = "Sample 1",
              sample_id = sample1,
              group_value = as.numeric(as.character(group1)),
              taxon = taxon
            ) %>%
            bind_rows(
              event_data %>%
                mutate(
                  # Sample2 data
                  sample_type = "Sample 2",
                  sample_id = sample2,
                  group_value = as.numeric(as.character(group2)),
                  taxon = taxon
                )
            ) %>%
            # Add ID for pair identification
            mutate(
              pair_id = paste(sample1, sample2, sep = "-"),
              # X-axis: absolute value of numeric difference per pair (|group1 - group2|)
              x_difference = abs(as.numeric(as.character(group1)) - as.numeric(as.character(group2))),
              # Y-axis: original numeric value
              y_value = group_value
            )
          
          message("Pair data created: ", nrow(pair_data), " rows")
          
          # Numerical variable: use categorical palette (per taxon)
          message("Using categorical palette for taxon")
          
          # Generate colors per taxon
          unique_taxa <- unique(pair_data$taxon)
          n_taxa <- length(unique_taxa)
          message("Number of unique taxa: ", n_taxa)
          
          # Generate sufficient colors
          if(n_taxa <= 8) {
            color_values <- RColorBrewer::brewer.pal(max(3, n_taxa), "Set3")
          } else if(n_taxa <= 12) {
            color_values <- RColorBrewer::brewer.pal(n_taxa, "Set3")
          } else {
            color_values <- rainbow(n_taxa)
          }
          
          # Map colors per taxon
          names(color_values) <- unique_taxa
          message("Categorical palette created: ", length(color_values), " colors")
          
          # Generate summary information for each point
          summary_data <- pair_data %>%
            group_by(x_difference, sharing_type, sample1, sample2) %>%
            summarise(
              unique_pairs = n_distinct(paste(sample1, sample2)),
              unique_species = n_distinct(species),
              total_events = n(),
              dominant_taxon = {
                taxon_clean <- as.character(taxon[!is.na(taxon) & taxon != ""])
                if(length(taxon_clean) > 0) {
                  taxon_counts <- table(taxon_clean)
                  as.character(names(taxon_counts)[which.max(taxon_counts)])
                } else {
                  "Unknown"
                }
              },
              .groups = 'drop'
            ) %>%
            # Representative y value for each point (use first sample's value)
            left_join(
              pair_data %>% 
                filter(sample_type == "Sample 1") %>%
                select(x_difference, sample1, sample2, sharing_type, y_value) %>%
                distinct(),
              by = c("x_difference", "sample1", "sample2", "sharing_type")
            )
          
          # Numerical plot: use point size for event count, simple tooltip
          p <- ggplot(summary_data, aes(x = x_difference, y = y_value, color = dominant_taxon,
                                       size = total_events,
                                       text = paste0("Sample Pair: ", sample1, " ↔ ", sample2, "<br>",
                                                   "Unique Pairs: ", unique_pairs, "<br>",
                                                   "Unique Species: ", unique_species, "<br>",
                                                   "Total Events: ", total_events, "<br>",
                                                   "Dominant ", tools::toTitleCase(input$tax_level), ": ", dominant_taxon, "<br>",
                                                   "Value Difference: ", round(x_difference, 2)))) +
            geom_point(alpha = input$point_alpha) +
            scale_size_continuous(range = c(2, 10), name = "Event Count") +
            labs(title = plot_title,
                 x = paste("Value Difference (", data$meta_var, ")"),
                 y = paste("Age Values (", data$meta_var, ")")) +
            scale_color_manual(values = color_values, name = paste("Dominant", tools::toTitleCase(input$tax_level))) +
            theme_minimal()
          
          # Additional plot: strain sharing event count (log scale)
          count_plot <- event_data %>%
            group_by(taxon) %>%
            summarise(
              count = n(),
              .groups = 'drop'
            ) %>%
            ggplot(aes(x = taxon, y = count, fill = taxon)) +
            geom_bar(stat = "identity", alpha = 0.7) +
            scale_y_log10() +
            labs(title = "Strain Sharing Event Count (log10 scale)",
                 x = "Taxon",
                 y = "Count (log10)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
          
          # Additional plot: density plot by taxonomy unit
          density_plot <- ggplot(pair_data, aes(x = x_difference, fill = phylum, color = phylum)) +
            geom_density(alpha = 0.3) +
            labs(title = paste("Distribution of Age Differences by", data$tax_column),
                 x = paste("Age Difference (", data$meta_var, ")"),
                 y = "Density") +
            scale_fill_manual(values = color_values) +
            scale_color_manual(values = color_values) +
            theme_minimal() +
            theme(legend.position = "bottom")
        }
        
        # Apply plot style (safely)
        tryCatch({
        p <- apply_plot_style(p, input)
        }, error = function(e) {
          message("Error applying plot style: ", e$message)
        })
        
        # Return only the existing Group Comparison plot
          tryCatch({
          plotly::event_register(
            ggplotly(p, tooltip = "text", source = "group_comparison") %>%
              layout(
                height = input$plot_height,
                width = input$plot_width
              ),
            "plotly_click"
          )
          }, error = function(e) {
          message("Error creating plotly: ", e$message)
          # Return default plotly object
          plot_ly() %>% 
            add_text(x = 0.5, y = 0.5, text = "Error generating plot", 
                     textfont = list(size = 16, color = "red")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE))
        })
        
        
      }, error = function(e) {
        showNotification(paste("Plot error:", e$message), type = "error")
        return(plot_ly() %>% 
          add_text(x = 0.5, y = 0.5, text = paste("Error:", e$message), 
                   textfont = list(size = 16, color = "red")) %>%
          layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      })
    })
    
    
        # Taxon-Category Heatmap (categorical variables only)
    output$heatmap_plot <- renderPlotly({
      cat("=== HEATMAP DEBUG: Starting heatmap rendering ===\n")

      # Check if analysis has been run
      if (is.null(input$update_analysis) || !shiny::isTruthy(input$metadata_var)) {
        cat("=== HEATMAP DEBUG: No analysis run yet ===\n")
        return(plot_ly() %>%
          add_text(x = 0.5, y = 0.5, text = "Please run analysis first:\n1. Select a metadata variable\n2. Click 'Update Analysis'",
                   textfont = list(size = 14, color = "orange")) %>%
          layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      data <- tryCatch({
        processed_data()
      }, error = function(e) {
        cat("=== HEATMAP DEBUG: Error getting processed_data:", e$message, "===\n")
        return(NULL)
      })

      req(data)

      cat("=== HEATMAP DEBUG: Data loaded ===\n")
      cat("Var type:", data$var_type, "\n")

      # Only process categorical variables
      if(data$var_type != "categorical") {
        cat("=== HEATMAP DEBUG: Not categorical ===\n")
        meta_info <- metadata_info()
        available_categorical <- if(length(meta_info$categorical) > 0) {
          paste("•", paste(meta_info$categorical, collapse = "\n• "))
        } else {
          "No categorical variables available"
        }

        return(plot_ly() %>%
          add_text(x = 0.5, y = 0.5,
                   text = paste("❌ Heatmap requires CATEGORICAL variables\n\n",
                               "Current selection:", input$metadata_var, "(", data$var_type, ")\n\n",
                               "📋 Available categorical variables:\n", available_categorical),
                   textfont = list(size = 12, color = "#dc3545")) %>%
          layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      tryCatch({
        cat("=== HEATMAP DEBUG: Processing categorical data ===\n")

        # Validate shared_events
        if (is.null(data$shared_events) || nrow(data$shared_events) == 0) {
          cat("=== HEATMAP DEBUG: No shared events ===\n")
          return(plot_ly() %>%
            add_text(x = 0.5, y = 0.5,
                     text = "No shared events to display.\n\nTry:\n- Lowering ANI threshold\n- Different taxonomic level",
                     textfont = list(size = 14, color = "#6c757d")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }

        # Check required columns
        required_cols <- c("sharing_type", "taxon")
        missing_cols <- setdiff(required_cols, colnames(data$shared_events))
        if (length(missing_cols) > 0) {
          cat("=== HEATMAP DEBUG: Missing columns:", paste(missing_cols, collapse = ", "), "===\n")
          return(plot_ly() %>%
            add_text(x = 0.5, y = 0.5,
                     text = paste("Data error: Missing", paste(missing_cols, collapse = ", ")),
                     textfont = list(size = 14, color = "red")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }

        # Check grouping mode
        grouping_mode <- if(is.null(data$grouping_mode) || is.na(data$grouping_mode)) {
          "inter_intra"
        } else {
          data$grouping_mode
        }
        cat("=== HEATMAP DEBUG: Grouping mode =", grouping_mode, "===\n")

        # Check sharing_type in actual data
        unique_sharing_types <- unique(data$shared_events$sharing_type)
        cat("Unique sharing types:", paste(unique_sharing_types, collapse = ", "), "\n")

        # Organize sharing_type into clear categories
        all_sharing_types <- sort(unique(data$shared_events$sharing_type))
        all_taxa <- sort(unique(data$shared_events$taxon))

        cat("Total sharing types:", length(all_sharing_types), "\n")
        cat("Total taxa:", length(all_taxa), "\n")

        # Calculate actual count
        heatmap_data <- data$shared_events %>%
          count(sharing_type, taxon, name = "count") %>%
          complete(
            sharing_type = all_sharing_types,
            taxon = all_taxa,
            fill = list(count = 0)
          ) %>%
          pivot_wider(
            names_from = taxon,
            values_from = count,
            values_fill = 0
          )

        cat("Heatmap data dimensions:", nrow(heatmap_data), "x", ncol(heatmap_data), "\n")

        if(nrow(heatmap_data) == 0) {
          return(plot_ly() %>%
            add_text(x = 0.5, y = 0.5, text = "No data for heatmap",
                     textfont = list(size = 16, color = "red")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }

        # Generate matrix - keep sharing_type as clear string
        sharing_types_labels <- heatmap_data$sharing_type
        heatmap_matrix <- as.matrix(heatmap_data[, -1, drop = FALSE])
        rownames(heatmap_matrix) <- sharing_types_labels

        cat("Matrix dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix), "\n")
        cat("Matrix range:", min(heatmap_matrix), "-", max(heatmap_matrix), "\n")
        cat("Row names (sharing types):", paste(rownames(heatmap_matrix), collapse = ", "), "\n")

        # Remove empty rows/columns
        row_sums <- rowSums(heatmap_matrix)
        col_sums <- colSums(heatmap_matrix)

        cat("Rows with data:", sum(row_sums > 0), "/", length(row_sums), "\n")
        cat("Cols with data:", sum(col_sums > 0), "/", length(col_sums), "\n")

        if(sum(row_sums > 0) == 0 || sum(col_sums > 0) == 0) {
          return(plot_ly() %>%
            add_text(x = 0.5, y = 0.5,
                     text = "All values are zero.\n\nCheck:\n- ANI threshold is not too high\n- Selected metadata has shared events",
                     textfont = list(size = 14, color = "orange")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }

        heatmap_matrix <- heatmap_matrix[row_sums > 0, col_sums > 0, drop = FALSE]

        # Display only top 50 taxa (for readability)
        if(ncol(heatmap_matrix) > 50) {
          col_order <- order(colSums(heatmap_matrix), decreasing = TRUE)
          heatmap_matrix <- heatmap_matrix[, head(col_order, 50), drop = FALSE]
          cat("Limited to top 50 taxa for readability\n")
        }

        cat("Final matrix dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix), "\n")
        cat("Final row names:", paste(rownames(heatmap_matrix), collapse = ", "), "\n")

        # Prepare Y-axis labels (clear category names)
        y_labels <- rownames(heatmap_matrix)

        # Calculate Y-axis height (sufficient space for each category)
        num_categories <- length(y_labels)
        dynamic_height <- max(400, num_categories * 80)  # 80px per category

        cat("Y-axis labels:", paste(y_labels, collapse = " | "), "\n")
        cat("Dynamic height:", dynamic_height, "\n")

        # Store current heatmap info for click events
        current_heatmap_info <<- list(
          sharing_types = y_labels,
          taxa = colnames(heatmap_matrix),
          matrix = heatmap_matrix
        )

        # Generate Plotly heatmap
        p <- plot_ly(
          z = heatmap_matrix,
          x = colnames(heatmap_matrix),
          y = y_labels,  # Use clear category names
          type = "heatmap",
          colorscale = "Viridis",
          hovertemplate = paste0(
            "<b>Taxon:</b> %{x}<br>",
            "<b>Sharing Type:</b> %{y}<br>",
            "<b>Event Count:</b> %{z}<br>",
            "<extra></extra>"
          ),
          source = "heatmap_plot",
          xgap = 2,
          ygap = 2,
          # Add customdata for click functionality
          customdata = array(
            rep(y_labels, each = ncol(heatmap_matrix)),
            dim = dim(heatmap_matrix)
          )
        ) %>%
        layout(
          title = list(
            text = paste0(
              "<b>Sharing Type vs Taxon Heatmap</b><br>",
              "<sub style='font-size:11px'>Grouping Mode: ",
              ifelse(grouping_mode == "inter_intra", "Inter/Intra Groups", "Original Categories"),
              " | Click cells to view details</sub>"
            ),
            font = list(size = 16)
          ),
          xaxis = list(
            title = paste("<b>Taxon (", data$tax_level, ")</b>", sep = ""),
            tickangle = 45,
            tickfont = list(size = 10),
            side = "bottom",
            automargin = TRUE
          ),
          yaxis = list(
            title = "<b>Sharing Type</b>",
            tickfont = list(size = 12),
            automargin = TRUE,
            type = "category",  # Specify Y-axis as category type
            categoryorder = "array",
            categoryarray = y_labels,  # Specify exact order
            tickmode = "array",
            tickvals = seq(0, length(y_labels) - 1),  # Start from 0
            ticktext = y_labels,  # Clear text labels
            dtick = 1  # Display for each category
          ),
          width = input$plot_width,
          height = dynamic_height,  # Dynamic height
          margin = list(
            l = 200,  # Increase left margin (Y-axis label space)
            r = 50,
            t = 100,
            b = 150   # Increase bottom margin (X-axis label space)
          ),
          hoverlabel = list(
            bgcolor = "white",
            font = list(size = 12)
          )
        )

        cat("=== HEATMAP DEBUG: Successfully created heatmap ===\n")
        return(plotly::event_register(p, "plotly_click"))

      }, error = function(e) {
        cat("=== HEATMAP DEBUG: ERROR ===\n")
        cat("Error:", e$message, "\n")
        cat("Call:", deparse(e$call), "\n")

        showNotification(paste("Heatmap error:", e$message), type = "error")

        return(plot_ly() %>%
          add_text(x = 0.5, y = 0.5,
                   text = paste0("Error generating heatmap:\n", e$message, "\n\nCheck console for details"),
                   textfont = list(size = 12, color = "red")) %>%
          layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      })
    })
    
    # Summary table
    output$summary_table <- DT::renderDataTable({
      data <- processed_data()
      req(data)
      
      tryCatch({
        # Filter only when specific taxon is selected and not "all"
        if(!is.null(input$selected_taxon) && input$selected_taxon != "" && input$selected_taxon != "all") {
          display_data <- data$event_prevalence_data %>%
            filter(taxon == input$selected_taxon) %>%
          mutate(
              event_prevalence = round(event_prevalence, 2),
              shared_events_count = as.integer(shared_events_count),
              total_events = as.integer(total_events),
              unique_species = as.integer(unique_species)
            ) %>%
            arrange(desc(event_prevalence))
        } else {
          # Show all data when "All" is selected or by default
          display_data <- data$event_prevalence_data %>%
            mutate(
              event_prevalence = round(event_prevalence, 2),
              shared_events_count = as.integer(shared_events_count),
              total_events = as.integer(total_events),
              unique_species = as.integer(unique_species)
            ) %>%
            arrange(desc(event_prevalence))
        }
        
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
        ) %>%
        DT::formatStyle(
          'event_prevalence',
          background = DT::styleColorBar(range(display_data$event_prevalence), '#2FA4E7'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
        
      }, error = function(e) {
        showNotification(paste("Table error:", e$message), type = "error")
        return(DT::datatable(data.frame(Error = "Failed to generate table")))
      })
    }, options = list(pageLength = 15, scrollX = TRUE))
    
    # Raw ANI-filtered table (show near-original rows with minimal derived columns)
    output$raw_table <- DT::renderDataTable({
      data <- processed_data()
      req(data)
      tryCatch({
        fd <- data$filtered_data
        validate(need(!is.null(fd) && nrow(fd) > 0, "No ANI-filtered records to display"))
        
        tax_col <- data$tax_level
        base_var <- data$meta_var
        meta_col_1 <- paste0(base_var, "_sample1")
        meta_col_2 <- paste0(base_var, "_sample2")
        if(!all(c(meta_col_1, meta_col_2) %in% colnames(fd))) {
          if(base_var == "disease_group" && all(c("disease_group_sample1", "disease_group_sample2") %in% colnames(fd))) {
            meta_col_1 <- "disease_group_sample1"
            meta_col_2 <- "disease_group_sample2"
          }
        }
        
        # Enrich with taxonomy if missing using GTDB cache (by accession)
        if(any(c("genus","species","domain","phylum","class","order","family") %in% colnames(fd)) == FALSE ||
           ("genus" %in% colnames(fd) && anyNA(fd$genus)) || ("species" %in% colnames(fd) && anyNA(fd$species))) {
          gtdb_data <- get_gtdb_cached()
          # attempt to enrich using accession columns if available
          acc_cols <- intersect(colnames(fd), c("clean_accession","accession","clean_genome","genome","genome_file"))
          if(length(acc_cols) > 0) {
            # derive root accession from any available accession-like column
            extract_root <- function(x) {
              x <- as.character(x)
              m <- regexpr("GC[FA]_\\d+\\.\\d+", x, perl = TRUE)
              v <- regmatches(x, m)
              v[is.na(v) | v == ""] <- NA_character_
              return(v)
            }
            fd$root_accession_tmp <- extract_root(do.call(coalesce, lapply(fd[acc_cols], as.character)))
            fd$accession_core_tmp <- sub("^GC[FA]_", "", fd$root_accession_tmp)
            fd <- fd %>% left_join(gtdb_data %>% select(accession_core, domain, phylum, class, order, family, genus, species),
                                   by = c("accession_core_tmp" = "accession_core")) %>%
              mutate(
                genus = coalesce(genus.x, genus.y),
                species = coalesce(species.x, species.y),
                domain = coalesce(domain.x, domain.y),
                phylum = coalesce(phylum.x, phylum.y),
                class = coalesce(class.x, class.y),
                order = coalesce(order.x, order.y),
                family = coalesce(family.x, family.y)
              ) %>%
              select(-ends_with(".x"), -ends_with(".y"), -root_accession_tmp, -accession_core_tmp)
          }
        }

        # Add minimal derived fields for grouping clarity
        display_data <- fd %>%
          mutate(
            group1 = .data[[meta_col_1]],
            group2 = .data[[meta_col_2]],
            taxon = .data[[tax_col]]
          ) %>%
          mutate(
            group1_chr = trimws(as.character(group1)),
            group2_chr = trimws(as.character(group2)),
            sharing_type = if_else(group1_chr == group2_chr, "Intra-group", "Inter-group"),
            sharing_label = if_else(
              group1_chr == group2_chr,
              paste0("Intra-group (", group1_chr, ")"),
              paste0("Inter-group (", sort(c(group1_chr, group2_chr))[1], " ↔ ", sort(c(group1_chr, group2_chr))[2], ")")
            ),
            metadata_combination = paste0(
              sort(c(group1_chr, group2_chr))[1], " ↔ ", sort(c(group1_chr, group2_chr))[2]
            )
          ) %>%
          select(-group1_chr, -group2_chr)
        
        # Reorder columns for readability while preserving raw columns
        front_cols <- c("clean_name1", "clean_name2", "species", tax_col, "popANI", "conANI", "sharing_type", "sharing_label", "metadata_combination")
        front_cols <- intersect(front_cols, colnames(display_data))
        other_cols <- setdiff(colnames(display_data), front_cols)
        display_data <- display_data[, c(front_cols, other_cols), drop = FALSE]
        
        DT::datatable(
          display_data,
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            searchHighlight = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel', 'colvis')
          ),
          rownames = FALSE,
          extensions = 'Buttons'
        )
      }, error = function(e) {
        showNotification(paste("Raw table error:", e$message), type = "error")
        return(DT::datatable(data.frame(Error = "Failed to generate raw table")))
      })
    })

    # ========== Download Modal Handlers ==========
    # Store reactive ggplot objects for each plot

    # Overview Events Plot ggplot
    overview_events_ggplot <- reactive({
      data <- processed_data()
      req(data)

      event_data <- data$event_prevalence_data
      req(nrow(event_data) > 0)

      # Create summary by sharing type
      summary_data <- event_data %>%
        group_by(sharing_type) %>%
        summarise(
          total_events = sum(shared_events_count, na.rm = TRUE),
          unique_taxa = n_distinct(taxon),
          .groups = "drop"
        )

      p <- ggplot(summary_data, aes(x = sharing_type, y = total_events, fill = sharing_type)) +
        geom_col(alpha = 0.8) +
        geom_text(aes(label = total_events), vjust = -0.3, size = 4) +
        labs(
          title = "Sharing Events by Type",
          x = "Sharing Type",
          y = "Total Events"
        ) +
        theme_minimal() +
        theme(legend.position = "none")

      apply_plot_style(p, input)
    })

    # Overview Taxonomy Plot ggplot
    overview_taxonomy_ggplot <- reactive({
      data <- processed_data()
      req(data)

      event_data <- data$event_prevalence_data
      req(nrow(event_data) > 0)

      # Top 10 taxa by event count
      top_taxa <- event_data %>%
        group_by(taxon) %>%
        summarise(total_events = sum(shared_events_count, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(total_events)) %>%
        slice_head(n = 10)

      p <- ggplot(top_taxa, aes(x = reorder(taxon, total_events), y = total_events, fill = taxon)) +
        geom_col(alpha = 0.8) +
        coord_flip() +
        labs(
          title = paste("Top 10 Taxa (", data$tax_level, ")", sep = ""),
          x = "Taxon",
          y = "Total Sharing Events"
        ) +
        theme_minimal() +
        theme(legend.position = "none")

      apply_plot_style(p, input)
    })

    # Group Comparison Plot ggplot
    group_comparison_ggplot <- reactive({
      data <- processed_data()
      req(data)

      event_data <- data$event_prevalence_data
      req(nrow(event_data) > 0)

      var_type <- data$var_type

      if (var_type == "categorical") {
        # For categorical: box plot or bar chart
        p <- ggplot(event_data, aes(x = sharing_type, y = event_prevalence, fill = sharing_type)) +
          geom_boxplot(alpha = 0.7) +
          geom_jitter(width = 0.2, alpha = 0.3, size = input$point_size) +
          labs(
            title = "Event Prevalence by Sharing Type",
            x = "Sharing Type",
            y = "Event Prevalence (%)"
          ) +
          theme_minimal()
      } else {
        # For numerical
        p <- ggplot(event_data, aes(x = taxon, y = event_prevalence, fill = sharing_type)) +
          geom_col(position = "dodge", alpha = 0.8) +
          coord_flip() +
          labs(
            title = "Event Prevalence by Taxon",
            x = "Taxon",
            y = "Event Prevalence (%)"
          ) +
          theme_minimal()
      }

      apply_plot_style(p, input)
    })

    # Density Plot ggplot
    density_ggplot <- reactive({
      data <- processed_data()
      req(data)

      event_data <- data$shared_events
      req(nrow(event_data) > 0)

      count_data <- event_data %>%
        count(taxon) %>%
        arrange(desc(n)) %>%
        slice_head(n = 20)

      p <- ggplot(count_data, aes(x = reorder(taxon, n), y = n, fill = taxon)) +
        geom_col(alpha = 0.7) +
        coord_flip() +
        labs(
          title = "Taxonomy Count Distribution",
          x = "Taxon",
          y = "Count"
        ) +
        theme_minimal() +
        theme(legend.position = "none")

      apply_plot_style(p, input)
    })

    # Stacked Bar Plot ggplot
    stacked_bar_ggplot <- reactive({
      data <- processed_data()
      req(data)

      event_data <- data$shared_events
      req(nrow(event_data) > 0)
      req(data$var_type == "categorical")

      # Get top taxa
      top_taxa <- event_data %>%
        count(taxon) %>%
        arrange(desc(n)) %>%
        slice_head(n = 15) %>%
        pull(taxon)

      stacked_data <- event_data %>%
        filter(taxon %in% top_taxa) %>%
        count(taxon, metadata_combination)

      p <- ggplot(stacked_data, aes(x = reorder(taxon, n), y = n, fill = metadata_combination)) +
        geom_col(alpha = 0.8) +
        coord_flip() +
        labs(
          title = "Taxa Distribution by Group",
          x = "Taxon",
          y = "Count",
          fill = "Group"
        ) +
        theme_minimal() +
        scale_fill_manual(values = get_palette(input$color_palette))

      apply_plot_style(p, input)
    })

    # Heatmap ggplot - matches the Plotly heatmap structure
    heatmap_ggplot <- reactive({
      data <- processed_data()
      req(data)
      req(data$var_type == "categorical")

      shared_events <- data$shared_events
      req(!is.null(shared_events) && nrow(shared_events) > 0)

      # Count by sharing_type and taxon (same as Plotly version)
      all_sharing_types <- sort(unique(shared_events$sharing_type))
      all_taxa <- sort(unique(shared_events$taxon))

      heatmap_data <- shared_events %>%
        count(sharing_type, taxon, name = "count") %>%
        tidyr::complete(
          sharing_type = all_sharing_types,
          taxon = all_taxa,
          fill = list(count = 0)
        )

      # Filter out empty rows and columns
      non_zero_taxa <- heatmap_data %>%
        group_by(taxon) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        filter(total > 0) %>%
        pull(taxon)

      non_zero_sharing <- heatmap_data %>%
        group_by(sharing_type) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        filter(total > 0) %>%
        pull(sharing_type)

      heatmap_data <- heatmap_data %>%
        filter(taxon %in% non_zero_taxa, sharing_type %in% non_zero_sharing)

      # Limit to top 50 taxa by total count (same as Plotly)
      if (length(unique(heatmap_data$taxon)) > 50) {
        top_taxa <- heatmap_data %>%
          group_by(taxon) %>%
          summarise(total = sum(count), .groups = "drop") %>%
          arrange(desc(total)) %>%
          slice_head(n = 50) %>%
          pull(taxon)

        heatmap_data <- heatmap_data %>%
          filter(taxon %in% top_taxa)
      }

      # Order taxa by total count
      taxa_order <- heatmap_data %>%
        group_by(taxon) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        arrange(desc(total)) %>%
        pull(taxon)

      heatmap_data$taxon <- factor(heatmap_data$taxon, levels = taxa_order)

      # Get grouping mode for title
      grouping_mode <- if (is.null(data$grouping_mode) || is.na(data$grouping_mode)) {
        "inter_intra"
      } else {
        data$grouping_mode
      }
      mode_label <- ifelse(grouping_mode == "inter_intra", "Inter/Intra Groups", "Original Categories")

      p <- ggplot(heatmap_data, aes(x = taxon, y = sharing_type, fill = count)) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_viridis_c(option = "viridis", name = "Event Count") +
        labs(
          title = paste0("Sharing Type vs Taxon Heatmap\n(", mode_label, ")"),
          x = paste("Taxon (", data$tax_level, ")", sep = ""),
          y = "Sharing Type"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )

      apply_plot_style(p, input)
    })

    # Setup download modal handlers for each plot
    setup_download_modal_handler(input, output, session, "overview_events", overview_events_ggplot, "shared_taxa_overview_events")
    setup_download_modal_handler(input, output, session, "overview_taxonomy", overview_taxonomy_ggplot, "shared_taxa_overview_taxonomy")
    setup_download_modal_handler(input, output, session, "group_comparison", group_comparison_ggplot, "shared_taxa_group_comparison")
    setup_download_modal_handler(input, output, session, "density", density_ggplot, "shared_taxa_density")
    setup_download_modal_handler(input, output, session, "stacked_bar", stacked_bar_ggplot, "shared_taxa_stacked_bar")
    setup_download_modal_handler(input, output, session, "heatmap", heatmap_ggplot, "shared_taxa_heatmap")

    # Density plot (newly added)
    output$density_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      
      tryCatch({
        # Check variable type
        var_type <- data$var_type
        event_data <- data$shared_events
        
        if(var_type == "categorical") {
          # Categorical: Count plot (sorted by number count) + Stacked bar chart

          # 1. Count plot (sorted)
          count_data <- event_data %>%
            count(taxon) %>%
            arrange(desc(n))
          
          count_plot <- count_data %>%
            ggplot(aes(x = reorder(taxon, n), y = n, fill = taxon)) +
            geom_col(alpha = 0.7) +
            labs(title = "Taxonomy Count Distribution (Sorted)",
                 x = "Taxonomy",
                 y = "Count") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none") +
            coord_flip()
          
        } else {
          # Numerical: Density plot by taxonomy
          
          # Use all 2 values per pair (include both group1 and group2)
          density_data <- event_data %>%
            select(taxon, group1, group2) %>%
            pivot_longer(cols = c(group1, group2), 
                        names_to = "position", 
                        values_to = "value") %>%
            mutate(value = as.numeric(as.character(value))) %>%
            filter(!is.na(value))
          
          # Load GTDB data to add phylum information (using data parsed from preprocessing_scripts.R)
          tryCatch({
            gtdb_file <- "data/filtered_gtdb_metadata.rds"
            if(file.exists(gtdb_file)) {
              gtdb_data <- readRDS(gtdb_file)
              message("GTDB data loaded: ", nrow(gtdb_data), " records")
              message("Available GTDB columns: ", paste(colnames(gtdb_data), collapse = ", "))
              
              # Create mapping table based on selected taxonomy level
              tax_level <- data$tax_level
              message("Current taxonomy level: ", tax_level)
              
              # Create phylum mapping table for each taxonomy level
              tax_to_phylum_mapping <- switch(tax_level,
                "species" = gtdb_data %>% 
                  select(species, phylum) %>% 
                  distinct() %>% 
                  filter(!is.na(species) & !is.na(phylum)),
                "genus" = gtdb_data %>% 
                  select(genus, phylum) %>% 
                  distinct() %>% 
                  filter(!is.na(genus) & !is.na(phylum)),
                "family" = gtdb_data %>% 
                  select(family, phylum) %>% 
                  distinct() %>% 
                  filter(!is.na(family) & !is.na(phylum)),
                "order" = gtdb_data %>% 
                  select(order, phylum) %>% 
                  distinct() %>% 
                  filter(!is.na(order) & !is.na(phylum)),
                "class" = gtdb_data %>% 
                  select(class, phylum) %>% 
                  distinct() %>% 
                  filter(!is.na(class) & !is.na(phylum)),
                # Default: species mapping
                gtdb_data %>% 
                  select(species, phylum) %>% 
                  distinct() %>% 
                  filter(!is.na(species) & !is.na(phylum))
              )
              
              # Standardize mapping table column names
              colnames(tax_to_phylum_mapping) <- c("taxon_name", "phylum")
              
              message("Mapping table created: ", nrow(tax_to_phylum_mapping), " entries")
              message("Sample mapping entries: ", paste(head(tax_to_phylum_mapping$taxon_name, 5), collapse = ", "))
              
              # Check taxon names before mapping
              unique_taxons <- unique(density_data$taxon)
              message("Unique taxons in data: ", paste(head(unique_taxons, 10), collapse = ", "))
              message("Total unique taxons: ", length(unique_taxons))
              
              # Check matching status
              matched_taxons <- intersect(unique_taxons, tax_to_phylum_mapping$taxon_name)
              message("Matched taxons: ", length(matched_taxons), " out of ", length(unique_taxons))
              if(length(matched_taxons) > 0) {
                message("Examples of matched: ", paste(head(matched_taxons, 5), collapse = ", "))
              }
              
              unmatched_taxons <- setdiff(unique_taxons, tax_to_phylum_mapping$taxon_name)
              if(length(unmatched_taxons) > 0) {
                message("Examples of unmatched: ", paste(head(unmatched_taxons, 5), collapse = ", "))
              }
              
              # Map with taxon
              density_data <- density_data %>%
                left_join(tax_to_phylum_mapping, by = c("taxon" = "taxon_name")) %>%
                mutate(phylum = ifelse(is.na(phylum), "Unknown", phylum))
              
              message("Phylum mapping completed. Unique phylums: ", paste(unique(density_data$phylum), collapse = ", "))
              message("Mapping success rate: ", round(100 * sum(!is.na(density_data$phylum) & density_data$phylum != "Unknown") / nrow(density_data), 1), "%")
              
            } else {
              message("GTDB file not found: ", gtdb_file)
              density_data$phylum <- "Unknown"
            }
          }, error = function(e) {
            message("Error loading GTDB data: ", e$message)
            density_data$phylum <- "Unknown"
          })
          
          # Determine legend name based on selected taxonomy level
          selected_tax_level <- data$tax_level
          
          # Create plot_data using phylum-mapped density_data
          # First create taxon-phylum mapping table
          taxon_phylum_map <- density_data %>%
            group_by(taxon) %>%
            slice(1) %>%  # Keep only one phylum info per taxon
            ungroup() %>%
            select(taxon, phylum) %>%
            filter(!is.na(taxon))
          
          message("Taxon-Phylum mapping table: ", nrow(taxon_phylum_map), " entries")
          message("Sample taxon-phylum pairs: ")
          if(nrow(taxon_phylum_map) > 0) {
            sample_pairs <- head(taxon_phylum_map, 5)
            for(i in 1:nrow(sample_pairs)) {
              message("  - ", sample_pairs$taxon[i], " -> ", sample_pairs$phylum[i])
            }
          }
          
          # Check taxons in event_data
          event_taxons <- unique(event_data$taxon)
          message("Event data taxons: ", paste(head(event_taxons, 5), collapse = ", "))
          
          # Check matching status
          matched_event_taxons <- intersect(event_taxons, taxon_phylum_map$taxon)
          message("Matched event taxons: ", length(matched_event_taxons), " out of ", length(event_taxons))
          
          plot_data <- event_data %>%
            left_join(taxon_phylum_map, by = "taxon") %>%
            mutate(
              value_diff = abs(as.numeric(as.character(group1)) - as.numeric(as.character(group2))),
              # Additional info for click events
              pair_info = paste0("Sample1: ", sample1, " (", group1, ")", "<br>",
                               "Sample2: ", sample2, " (", group2, ")", "<br>",
                               "Species: ", species, "<br>",
                               "Phylum: ", ifelse(is.na(phylum), "Unknown", phylum), "<br>",
                               "Difference: ", round(value_diff, 2)),
              phylum = ifelse(is.na(phylum), "Unknown", phylum)
            ) %>%
            filter(!is.na(value_diff))
          
          message("Plot data prepared: ", nrow(plot_data), " events")
          message("Unique phylums in plot data: ", paste(unique(plot_data$phylum), collapse = ", "))
          message("Phylum distribution in plot data:")
          phylum_counts <- table(plot_data$phylum)
          for(i in 1:length(phylum_counts)) {
            message("  - ", names(phylum_counts)[i], ": ", phylum_counts[i], " events")
          }
          
          count_plot <- plot_data %>%
            ggplot(aes(x = value_diff, fill = phylum, 
                      text = pair_info)) +
            geom_histogram(alpha = 0.7, bins = 30, position = "stack") +
            labs(title = paste("Strain Sharing Events by", data$meta_var, "Difference"),
                 x = paste(data$meta_var, "Difference (Absolute)"),
                 y = "Count (Events)",
                 fill = "Phylum") +
            theme_minimal() +
            theme(legend.position = "right") +
            guides(fill = guide_legend(override.aes = list(alpha = 1)))
        }
        
        tryCatch({
          count_plot <- apply_plot_style(count_plot, input)
        }, error = function(e) {
          message("Error applying plot style to density plot: ", e$message)
        })
        
        tryCatch({
          plotly::event_register(
            ggplotly(count_plot, tooltip = "text", source = "density_plot") %>%
              layout(
                height = 400,
                width = input$plot_width
              ),
            "plotly_click"
          )
        }, error = function(e) {
          message("Error creating density plotly: ", e$message)
          plot_ly() %>% 
            add_text(x = 0.5, y = 0.5, text = "Error generating density plot", 
                     textfont = list(size = 16, color = "red")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE))
        })
        
      }, error = function(e) {
        showNotification(paste("Density plot error:", e$message), type = "error")
        return(plot_ly() %>% 
          add_text(x = 0.5, y = 0.5, text = paste("Error:", e$message), 
                   textfont = list(size = 16, color = "red")) %>%
          layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      })
    })
    
    # Handle density plot click event
    observeEvent(event_data("plotly_click", source = "density_plot"), {
      click_data <- event_data("plotly_click", source = "density_plot")
      req(click_data)
      
      tryCatch({
        data <- processed_data()
        req(data)
        
        if(data$var_type == "numerical") {
          # Find value_diff range of clicked point
          clicked_x <- click_data$x
          
          # Set range (clicked point ± certain range)
          range_tolerance <- 2.0  # adjustable
          min_x <- clicked_x - range_tolerance
          max_x <- clicked_x + range_tolerance
          
          # Find shared events in this range
          event_data <- data$shared_events
          
          # Calculate value_diff and filter
          relevant_events <- event_data %>%
            mutate(value_diff = abs(as.numeric(as.character(group1)) - as.numeric(as.character(group2)))) %>%
            filter(value_diff >= min_x & value_diff <= max_x) %>%
            select(sample1, sample2, species, group1, group2, value_diff, sharing_type) %>%
            arrange(desc(value_diff))
          
          message("Density plot clicked at x=", clicked_x, ", found ", nrow(relevant_events), " events in range [", min_x, ", ", max_x, "]")
          
          if(nrow(relevant_events) > 0) {
            # Store data in reactive value
            clicked_data_rv$data <- relevant_events
            
            # Show table section
            shinyjs::show("clicked_data_section")
            
            showNotification(
              paste("Found", nrow(relevant_events), "shared events in the selected range"),
              type = "message",
              duration = 3
            )
          } else {
            showNotification("No shared events found in the selected range", type = "warning")
          }
        }
        
      }, error = function(e) {
        showNotification(paste("Density plot click error:", e$message), type = "error")
      })
    })
    
    # Stacked bar chart (categorical variables only)
    output$stacked_bar_plot <- renderPlotly({
      data <- processed_data()
      req(data)
      
      tryCatch({
        # Execute only for categorical variables
        if(data$var_type != "categorical") {
          return(plot_ly() %>% 
            add_text(x = 0.5, y = 0.5, text = "Stacked bar chart is only available for categorical variables", 
                     textfont = list(size = 16, color = "gray")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }
        
        event_data <- data$shared_events
        
        # Stacked bar chart: stacked by category metadata per taxon
        stacked_data <- event_data %>%
          count(taxon, sharing_type) %>%
          arrange(desc(n))
        
        stacked_plot <- stacked_data %>%
          ggplot(aes(x = reorder(taxon, n), y = n, fill = sharing_type)) +
          geom_col(position = "stack", alpha = 0.8) +
          labs(title = "Taxonomy Count by Category (Stacked)",
               x = "Taxonomy",
               y = "Count",
               fill = "Category") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          coord_flip()
        
        tryCatch({
          stacked_plot <- apply_plot_style(stacked_plot, input)
        }, error = function(e) {
          message("Error applying plot style to stacked plot: ", e$message)
        })
        
        tryCatch({
          ggplotly(stacked_plot, tooltip = "text") %>% 
            layout(
              height = 400,
              width = input$plot_width
            )
        }, error = function(e) {
          message("Error creating stacked plotly: ", e$message)
          plot_ly() %>% 
            add_text(x = 0.5, y = 0.5, text = "Error generating stacked plot", 
                     textfont = list(size = 16, color = "red")) %>%
            layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE))
        })
        
      }, error = function(e) {
        showNotification(paste("Stacked plot error:", e$message), type = "error")
        return(plot_ly() %>% 
          add_text(x = 0.5, y = 0.5, text = paste("Error:", e$message), 
                   textfont = list(size = 16, color = "red")) %>%
          layout(showlegend = FALSE, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      })
    })
  })
}