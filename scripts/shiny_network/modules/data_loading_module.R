# =============================================================================
# Data Loading and Overview Module
# Displays data overview and quality metrics
# =============================================================================

library(shiny)
library(shinyjs)
library(bslib)
library(plotly)
library(DT)
library(ggplot2)
library(phyloseq)
library(dplyr)

#' Data Loading UI
data_loading_UI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),

    card(
      full_screen = TRUE,
      class = "main-card",
      card_header(
        style = "background-color: #2FA4E7; color: white;",
        "Data Overview"
      ),
      card_body(
        # Message before data is loaded
        div(
          id = ns("load_message"),
          style = "color: red; font-weight: bold; margin-bottom: 20px;",
          icon("exclamation-triangle"),
          " Please load data using the Master Controller (click the arrow in upper left)"
        ),

        # Data summary cards
        layout_column_wrap(
          width = "250px",
          fill = FALSE,
          value_box(
            title = "Total Samples",
            value = textOutput(ns("n_samples")),
            showcase = bs_icon("diagram-3"),
            theme = "success"
          ),
          value_box(
            title = "Total Taxa",
            value = textOutput(ns("n_taxa")),
            showcase = bs_icon("grid-3x3"),
            theme = "info"
          ),
          value_box(
            title = "Metadata Columns",
            value = textOutput(ns("n_metadata")),
            showcase = bs_icon("table"),
            theme = "warning"
          )
        ),

        br(),

        # Taxonomic prevalence vs abundance plot with rank selector
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header(
            style = "background-color: #2FA4E7; color: white;",
            "Taxonomic Prevalence vs Mean Abundance"
          ),
          card_body(
            style = "height: 600px;",
            # Taxonomic rank selector
            div(
              style = "margin-bottom: 15px;",
              selectInput(
                ns("tax_rank_selector"),
                "Select Taxonomic Rank:",
                choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                selected = "Species",
                width = "250px"
              )
            ),
            plotlyOutput(ns("tax_prevalence_abundance"), height = "100%")
          )
        ),

        br(),

        # Taxonomy overview table
        card(
          full_screen = TRUE,
          fill = FALSE,
          card_header(
            style = "background-color: #2FA4E7; color: white;",
            "Top Taxa Overview"
          ),
          card_body(
            DTOutput(ns("tax_table"), height = "900px")
          )
        )
      )
    )
  )
}

#' Data Loading Server
data_loading_Server <- function(id, phyloseq_reactive, data_loaded_reactive, tax_rank_reactive = reactive("Species")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Taxonomic rank hierarchy (lower index = higher level)
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    # Show/hide load message
    observe({
      if (data_loaded_reactive()) {
        shinyjs::hide("load_message")
      } else {
        shinyjs::show("load_message")
      }
    })

    # Update taxonomic rank selector based on Master Controller selection
    observe({
      req(tax_rank_reactive())
      master_rank <- tax_rank_reactive()

      # Get allowed ranks (master rank and above only)
      idx <- which(tax_ranks == master_rank)
      if (length(idx) == 0) idx <- length(tax_ranks)
      allowed_ranks <- tax_ranks[1:idx]

      # Update tax_rank_selector selectInput
      updateSelectInput(session, "tax_rank_selector",
                       choices = allowed_ranks,
                       selected = master_rank)
    })

    # Number of samples
    output$n_samples <- renderText({
      req(phyloseq_reactive())
      ps <- phyloseq_reactive()
      as.character(nsamples(ps))
    })

    # Number of taxa
    output$n_taxa <- renderText({
      req(phyloseq_reactive())
      ps <- phyloseq_reactive()
      as.character(ntaxa(ps))
    })

    # Number of metadata columns
    output$n_metadata <- renderText({
      req(phyloseq_reactive())
      ps <- phyloseq_reactive()

      # Get sample data
      sample_df <- as(sample_data(ps), "data.frame")

      # Return number of metadata columns
      as.character(ncol(sample_df))
    })

    # Taxonomic prevalence vs abundance plot (responsive to rank selector)
    output$tax_prevalence_abundance <- renderPlotly({
      req(phyloseq_reactive(), input$tax_rank_selector)
      ps <- phyloseq_reactive()
      selected_rank <- input$tax_rank_selector

      # Get OTU table
      otu_mat <- as(otu_table(ps), "matrix")
      if (!taxa_are_rows(ps)) {
        otu_mat <- t(otu_mat)
      }

      # Get taxonomy table
      tax_mat <- as(tax_table(ps), "matrix")

      # Find the selected taxonomic rank column
      rank_col <- NULL
      for (col_name in colnames(tax_mat)) {
        if (tolower(col_name) == tolower(selected_rank)) {
          rank_col <- tax_mat[, col_name]
          break
        }
      }

      # If not found by exact match, try partial match
      if (is.null(rank_col)) {
        rank_idx <- grep(selected_rank, colnames(tax_mat), ignore.case = TRUE)
        if (length(rank_idx) > 0) {
          rank_col <- tax_mat[, rank_idx[1]]
        } else {
          return(plotly_empty() %>%
                  layout(title = paste("No", selected_rank, "column found in taxonomy table")))
        }
      }

      # Create data frame with taxon assignments
      df <- data.frame(
        taxon = rank_col,
        otu_mat,
        stringsAsFactors = FALSE
      )

      # Remove NA taxa
      df <- df[!is.na(df$taxon) & df$taxon != "", ]

      # Aggregate abundance to the selected taxonomic rank
      # Sum all ASVs/OTUs belonging to the same taxon across all samples
      # Also count how many ASVs/OTUs belong to each taxon
      tax_summary <- df %>%
        group_by(taxon) %>%
        summarise(
          across(where(is.numeric), sum),
          n_taxa = n(),  # Count number of ASVs/OTUs in this taxon
          .groups = "drop"
        )

      # Now calculate prevalence and mean abundance at the aggregated level
      abundance_matrix <- as.matrix(tax_summary[, !(names(tax_summary) %in% c("taxon", "n_taxa"))])
      rownames(abundance_matrix) <- tax_summary$taxon

      tax_summary <- data.frame(
        taxon = tax_summary$taxon,
        mean_abundance = rowMeans(abundance_matrix),
        prevalence = rowSums(abundance_matrix > 0) / ncol(abundance_matrix) * 100,  # As percentage
        n_taxa = tax_summary$n_taxa,  # Number of ASVs/OTUs aggregated
        stringsAsFactors = FALSE
      )

      # Filter for positive values (required for log scale)
      tax_summary <- tax_summary %>%
        filter(mean_abundance > 0, prevalence > 0)

      # Check if we have any data left
      if (nrow(tax_summary) == 0) {
        return(plotly_empty() %>%
                layout(title = paste("No data available for", selected_rank, "level")))
      }

      # Add small offset to avoid log10(0)
      tax_summary$mean_abundance <- tax_summary$mean_abundance + 1e-10
      tax_summary$prevalence <- tax_summary$prevalence + 0.1

      # Sort by mean abundance and keep top 10 for legend
      tax_summary <- tax_summary %>%
        arrange(desc(mean_abundance))

      # Create a factor for top 10 + "Other"
      if (nrow(tax_summary) > 10) {
        top_10_taxa <- tax_summary$taxon[1:10]
        tax_summary$display_taxon <- ifelse(
          tax_summary$taxon %in% top_10_taxa,
          as.character(tax_summary$taxon),
          "Other"
        )
        # Make "Other" points smaller and grey
        tax_summary$display_color <- tax_summary$display_taxon
      } else {
        tax_summary$display_taxon <- as.character(tax_summary$taxon)
        tax_summary$display_color <- tax_summary$display_taxon
      }

      # Create tooltip text
      tax_summary$tooltip_text <- paste0(
        "Taxon: ", tax_summary$taxon, "\n",
        "Prevalence: ", round(tax_summary$prevalence, 2), "%\n",
        "Mean Abundance: ", format(tax_summary$mean_abundance, scientific = TRUE, digits = 3), "\n",
        "Number of Taxa: ", tax_summary$n_taxa
      )

      # Create plot - simplified for plotly compatibility
      p <- ggplot(tax_summary, aes(x = prevalence, y = mean_abundance,
                                    color = display_color, size = n_taxa)) +
        geom_point(alpha = 0.7) +
        scale_x_log10(
          breaks = c(0.1, 1, 10, 50, 100),
          labels = c("0", "1", "10", "50", "100")
        ) +
        scale_y_log10() +
        labs(
          title = paste(selected_rank, "Level: Prevalence vs Mean Abundance"),
          x = "Prevalence (% of samples, log10 scale)",
          y = "Mean Abundance (log10 scale)",
          size = "Number of Taxa",
          color = selected_rank
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "right",
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10)
        )

      # Convert to plotly with error handling
      tryCatch({
        # Convert ggplot to plotly
        plotly_obj <- ggplotly(p, tooltip = c("x", "y", "color", "size"))

        # Add custom tooltips manually
        for (i in seq_along(plotly_obj$x$data)) {
          if (!is.null(plotly_obj$x$data[[i]]$text)) {
            # Map to our tooltip text
            plotly_obj$x$data[[i]]$text <- tax_summary$tooltip_text[i]
          }
        }

        plotly_obj %>%
          layout(
            height = 900,
            legend = list(
              orientation = "v",
              x = 1.02,
              y = 1,
              font = list(size = 10),
              title = list(font = list(size = 11))
            )
          )
      }, error = function(e) {
        # If conversion fails, create a simple plotly plot directly
        plot_ly(tax_summary,
                x = ~prevalence,
                y = ~mean_abundance,
                color = ~display_color,
                size = ~n_taxa,
                text = ~tooltip_text,
                type = "scatter",
                mode = "markers",
                marker = list(opacity = 0.7)) %>%
          layout(
            height = 650,
            title = paste(selected_rank, "Level: Prevalence vs Mean Abundance"),
            xaxis = list(
              title = "Prevalence (% of samples, log10 scale)",
              type = "log",
              tickvals = c(0.1, 1, 10, 50, 100),
              ticktext = c("0", "1", "10", "50", "100")
            ),
            yaxis = list(
              title = "Mean Abundance (log10 scale)",
              type = "log"
            ),
            legend = list(
              orientation = "v",
              x = 1.02,
              y = 1,
              font = list(size = 9)
            )
          )
      })
    })

    # Taxonomy table
    output$tax_table <- renderDT({
      req(phyloseq_reactive())
      ps <- phyloseq_reactive()

      # Get top 100 taxa by abundance
      otu_mat <- as(otu_table(ps), "matrix")
      if (!taxa_are_rows(ps)) {
        otu_mat <- t(otu_mat)
      }

      mean_abundance <- rowMeans(otu_mat)
      top_taxa <- names(sort(mean_abundance, decreasing = TRUE)[1:min(100, length(mean_abundance))])

      # Get taxonomy
      tax_mat <- as(tax_table(ps), "matrix")
      tax_df <- as.data.frame(tax_mat[top_taxa, , drop = FALSE])
      tax_df$MeanAbundance <- mean_abundance[top_taxa]
      tax_df$Prevalence <- rowSums(otu_mat[top_taxa, , drop = FALSE] > 0) / ncol(otu_mat) * 100

      # Reorder columns
      col_order <- c(colnames(tax_df)[!colnames(tax_df) %in% c("MeanAbundance", "Prevalence")],
                    "MeanAbundance", "Prevalence")
      tax_df <- tax_df[, col_order]

      datatable(
        tax_df,
        options = list(
          pageLength = 50,  # Increased from 20
          scrollX = TRUE,
          autoWidth = TRUE,
          order = list(list(ncol(tax_df) - 1, 'desc'))  # Sort by MeanAbundance
        ),
        rownames = TRUE,
        filter = 'top'
      ) %>%
        formatRound(columns = c("MeanAbundance", "Prevalence"), digits = 2)
    })
  })
}
