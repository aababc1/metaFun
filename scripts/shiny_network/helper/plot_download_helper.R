# =============================================================================
# Plot Download Helper Functions with Live Preview
# Provides modal dialogs for customizing and downloading plots
# =============================================================================

library(shiny)
library(ggplot2)
library(plotly)
library(RColorBrewer)

# Available themes
theme_options <- list(
  "Classic" = theme_classic,
  "Black & White" = theme_bw,
  "Minimal" = theme_minimal,
  "Gray" = theme_gray,
  "Light" = theme_light,
  "Dark" = theme_dark,
  "Void" = theme_void
)

# Available color palettes
color_palettes <- list(
  "Default" = NULL,
  "Set1" = "Set1",
  "Set2" = "Set2",
  "Set3" = "Set3",
  "Paired" = "Paired",
  "Dark2" = "Dark2",
  "Accent" = "Accent",
  "Pastel1" = "Pastel1",
  "Pastel2" = "Pastel2",
  "Spectral" = "Spectral",
  "RdYlBu" = "RdYlBu",
  "RdYlGn" = "RdYlGn",
  "PuOr" = "PuOr",
  "BrBG" = "BrBG",
  "Viridis" = "viridis",
  "Plasma" = "plasma",
  "Inferno" = "inferno",
  "Magma" = "magma"
)

#' Create a download button for individual plots
#' @param ns Namespace function
#' @param plot_id ID of the plot
#' @param btn_label Label for the button
#' @param btn_class Optional CSS class for the button
#' @return Action button for opening download modal
create_plot_download_btn <- function(ns, plot_id, btn_label = "Download",
                                     btn_class = "btn-sm btn-outline-primary") {
  actionButton(
    ns(paste0("download_btn_", plot_id)),
    label = tagList(icon("download"), btn_label),
    class = btn_class,
    style = "margin-left: 10px; float: right;"
  )
}

#' Create download modal dialog with full customization and live preview
#' @param ns Namespace function
#' @param plot_id ID of the plot
#' @param title Modal title
#' @return Modal dialog UI
create_download_modal_full <- function(ns, plot_id, title = "Download Plot") {
  modalDialog(
    title = tagList(icon("download"), " ", title),
    size = "xl",
    easyClose = TRUE,

    fluidRow(
      # Left panel - Settings
      column(4,
        div(
          style = "max-height: 600px; overflow-y: auto; padding-right: 10px;",

          # File Settings
          wellPanel(
            style = "background: #f8f9fa; padding: 12px; margin-bottom: 10px;",
            h5(icon("file-image"), " File Settings", style = "margin-top: 0; margin-bottom: 10px;"),

            selectInput(
              ns(paste0("modal_format_", plot_id)),
              "Format:",
              choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"),
              selected = "png",
              width = "100%"
            ),

            fluidRow(
              column(6,
                numericInput(ns(paste0("modal_width_", plot_id)), "Width (in):",
                            value = 10, min = 4, max = 20, step = 0.5, width = "100%")
              ),
              column(6,
                numericInput(ns(paste0("modal_height_", plot_id)), "Height (in):",
                            value = 8, min = 3, max = 15, step = 0.5, width = "100%")
              )
            ),

            sliderInput(ns(paste0("modal_dpi_", plot_id)), "DPI:",
                       min = 72, max = 600, value = 300, step = 24, width = "100%")
          ),

          # Theme & Style
          wellPanel(
            style = "background: #e8f4f8; padding: 12px; margin-bottom: 10px;",
            h5(icon("palette"), " Theme & Style", style = "margin-top: 0; margin-bottom: 10px;"),

            selectInput(
              ns(paste0("modal_theme_", plot_id)),
              "Theme:",
              choices = names(theme_options),
              selected = "Classic",
              width = "100%"
            ),

            selectInput(
              ns(paste0("modal_palette_", plot_id)),
              "Color Palette:",
              choices = names(color_palettes),
              selected = "Default",
              width = "100%"
            ),

            sliderInput(ns(paste0("modal_base_size_", plot_id)), "Base Font Size:",
                       min = 8, max = 24, value = 14, step = 1, width = "100%"),

            checkboxInput(ns(paste0("modal_show_grid_", plot_id)), "Show Grid Lines", value = TRUE)
          ),

          # Point & Line Settings
          wellPanel(
            style = "background: #f0f8e8; padding: 12px; margin-bottom: 10px;",
            h5(icon("circle"), " Points & Lines", style = "margin-top: 0; margin-bottom: 10px;"),

            sliderInput(ns(paste0("modal_point_size_", plot_id)), "Point Size:",
                       min = 0.5, max = 5, value = 2, step = 0.25, width = "100%"),

            sliderInput(ns(paste0("modal_point_alpha_", plot_id)), "Point Alpha:",
                       min = 0.1, max = 1, value = 0.7, step = 0.1, width = "100%"),

            sliderInput(ns(paste0("modal_line_size_", plot_id)), "Line Size:",
                       min = 0.25, max = 3, value = 1, step = 0.25, width = "100%")
          ),

          # Legend Settings
          wellPanel(
            style = "background: #fff8e8; padding: 12px; margin-bottom: 10px;",
            h5(icon("list"), " Legend", style = "margin-top: 0; margin-bottom: 10px;"),

            selectInput(
              ns(paste0("modal_legend_pos_", plot_id)),
              "Position:",
              choices = c("Right" = "right", "Left" = "left",
                         "Top" = "top", "Bottom" = "bottom", "None" = "none"),
              selected = "right",
              width = "100%"
            ),

            sliderInput(ns(paste0("modal_legend_size_", plot_id)), "Legend Text Size:",
                       min = 6, max = 16, value = 10, step = 1, width = "100%")
          ),

          # Action Buttons
          div(
            style = "margin-top: 15px;",
            actionButton(
              ns(paste0("preview_btn_", plot_id)),
              label = tagList(icon("eye"), " Update Preview"),
              class = "btn-info btn-block",
              style = "margin-bottom: 8px;"
            ),

            downloadButton(
              ns(paste0("modal_download_", plot_id)),
              label = tagList(icon("download"), " Download Plot"),
              class = "btn-primary btn-block"
            )
          )
        )
      ),

      # Right panel - Preview
      column(8,
        div(
          style = "border: 1px solid #ddd; border-radius: 4px; padding: 15px; background: white; min-height: 550px;",
          h5(icon("image"), " Preview", style = "margin-top: 0; color: #666; margin-bottom: 15px;"),
          div(
            style = "overflow: auto; max-height: 520px; text-align: center;",
            plotOutput(ns(paste0("modal_preview_", plot_id)), height = "auto")
          )
        )
      )
    ),

    footer = tagList(
      modalButton("Close")
    )
  )
}

#' Apply custom styling to a ggplot based on modal inputs
#' @param p ggplot object
#' @param input Shiny input
#' @param plot_id Plot ID for namespacing
#' @return Modified ggplot object
apply_modal_style <- function(p, input, plot_id) {
  # Get settings with defaults
  theme_name <- input[[paste0("modal_theme_", plot_id)]]
  base_size <- input[[paste0("modal_base_size_", plot_id)]]
  show_grid <- input[[paste0("modal_show_grid_", plot_id)]]
  legend_pos <- input[[paste0("modal_legend_pos_", plot_id)]]
  legend_size <- input[[paste0("modal_legend_size_", plot_id)]]
  palette_name <- input[[paste0("modal_palette_", plot_id)]]

  if (is.null(theme_name)) theme_name <- "Classic"
  if (is.null(base_size)) base_size <- 14
  if (is.null(show_grid)) show_grid <- TRUE
  if (is.null(legend_pos)) legend_pos <- "right"
  if (is.null(legend_size)) legend_size <- 10
  if (is.null(palette_name)) palette_name <- "Default"

  # Apply theme
  theme_fn <- theme_options[[theme_name]]
  if (!is.null(theme_fn)) {
    p <- p + theme_fn(base_size = base_size)
  }

  # Apply color palette
  palette_value <- color_palettes[[palette_name]]
  if (!is.null(palette_value)) {
    tryCatch({
      if (palette_value %in% c("viridis", "plasma", "inferno", "magma")) {
        # Viridis-style palettes
        p <- p + scale_fill_viridis_d(option = palette_value) +
                 scale_color_viridis_d(option = palette_value)
      } else {
        # Brewer palettes
        p <- p + scale_fill_brewer(palette = palette_value) +
                 scale_color_brewer(palette = palette_value)
      }
    }, error = function(e) {
      # If palette doesn't work (e.g., continuous data), try continuous versions
      tryCatch({
        if (palette_value %in% c("viridis", "plasma", "inferno", "magma")) {
          p <- p + scale_fill_viridis_c(option = palette_value) +
                   scale_color_viridis_c(option = palette_value)
        } else {
          p <- p + scale_fill_distiller(palette = palette_value) +
                   scale_color_distiller(palette = palette_value)
        }
      }, error = function(e2) NULL)
    })
  }

  # Apply additional theme modifications
  p <- p + theme(
    legend.position = legend_pos,
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size + 2)
  )

  # Remove grid if requested
  if (!show_grid) {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  return(p)
}

#' Server-side handler for full download modal with preview
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param plot_id ID of the plot
#' @param plot_fn Function that creates the base ggplot (receives input for customization)
#' @param filename_prefix Prefix for the downloaded file name
#' @param modal_title Title for the modal
setup_full_download_handler <- function(input, output, session, plot_id, plot_fn,
                                        filename_prefix = "plot", modal_title = "Download Plot") {
  ns <- session$ns

  # Open modal when download button clicked
  observeEvent(input[[paste0("download_btn_", plot_id)]], {
    showModal(create_download_modal_full(ns, plot_id, modal_title))

    # Trigger initial preview after modal opens
    shinyjs::delay(200, {
      shinyjs::click(paste0("preview_btn_", plot_id))
    })
  })

  # Create styled plot
  styled_plot <- reactive({
    # Re-evaluate when preview button is clicked
    input[[paste0("preview_btn_", plot_id)]]

    # Get base plot
    p <- tryCatch({
      plot_fn()
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(p)) return(NULL)

    # Get point/line settings
    point_size <- input[[paste0("modal_point_size_", plot_id)]]
    point_alpha <- input[[paste0("modal_point_alpha_", plot_id)]]
    line_size <- input[[paste0("modal_line_size_", plot_id)]]

    if (is.null(point_size)) point_size <- 2
    if (is.null(point_alpha)) point_alpha <- 0.7
    if (is.null(line_size)) line_size <- 1

    # Try to update geom defaults (this works for new geom calls)
    # For existing plots, we need to modify the layers
    tryCatch({
      # Modify existing point layers
      for (i in seq_along(p$layers)) {
        layer <- p$layers[[i]]
        if (inherits(layer$geom, "GeomPoint")) {
          p$layers[[i]]$aes_params$size <- point_size
          p$layers[[i]]$aes_params$alpha <- point_alpha
        }
        if (inherits(layer$geom, "GeomLine")) {
          p$layers[[i]]$aes_params$linewidth <- line_size
        }
      }
    }, error = function(e) NULL)

    # Apply theme and style
    p <- apply_modal_style(p, input, plot_id)

    return(p)
  })

  # Render preview
  output[[paste0("modal_preview_", plot_id)]] <- renderPlot({
    p <- styled_plot()
    req(p)
    p
  }, width = function() {
    w <- input[[paste0("modal_width_", plot_id)]]
    if (is.null(w)) w <- 10
    min(w * 55, 800)  # Scale for preview, max 800px
  }, height = function() {
    h <- input[[paste0("modal_height_", plot_id)]]
    if (is.null(h)) h <- 8
    min(h * 55, 600)
  })

  # Download handler
  output[[paste0("modal_download_", plot_id)]] <- downloadHandler(
    filename = function() {
      format <- input[[paste0("modal_format_", plot_id)]]
      if (is.null(format)) format <- "png"
      paste0(filename_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", format)
    },
    content = function(file) {
      p <- styled_plot()
      req(p)

      width <- input[[paste0("modal_width_", plot_id)]]
      height <- input[[paste0("modal_height_", plot_id)]]
      dpi <- input[[paste0("modal_dpi_", plot_id)]]
      format <- input[[paste0("modal_format_", plot_id)]]

      if (is.null(width)) width <- 10
      if (is.null(height)) height <- 8
      if (is.null(dpi)) dpi <- 300
      if (is.null(format)) format <- "png"

      ggsave(
        file,
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        device = format
      )
    }
  )

  # Return the styled plot reactive for potential reuse
  return(styled_plot)
}

#' Simplified setup for plotly downloads (HTML/PNG export)
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param plot_id ID of the plot
#' @param plotly_fn Function that creates the plotly object
#' @param data_fn Optional function that returns the underlying data
#' @param filename_prefix Prefix for the downloaded file name
setup_plotly_download_handler <- function(input, output, session, plot_id, plotly_fn,
                                          data_fn = NULL, filename_prefix = "plot") {
  ns <- session$ns

  # Open modal
  observeEvent(input[[paste0("download_btn_", plot_id)]], {
    showModal(modalDialog(
      title = tagList(icon("download"), " Download Plot"),
      size = "m",
      easyClose = TRUE,

      wellPanel(
        style = "background: #f8f9fa; padding: 15px;",

        selectInput(
          ns(paste0("plotly_format_", plot_id)),
          "Export Format:",
          choices = c(
            "PNG (static)" = "png",
            "HTML (interactive)" = "html",
            "CSV (data only)" = "csv"
          ),
          selected = "png"
        ),

        conditionalPanel(
          condition = sprintf("input['%s'] == 'png'", paste0("plotly_format_", plot_id)),
          ns = ns,
          numericInput(ns(paste0("plotly_width_", plot_id)), "Width (px):",
                      value = 1200, min = 400, max = 3000, step = 100),
          numericInput(ns(paste0("plotly_height_", plot_id)), "Height (px):",
                      value = 800, min = 300, max = 2000, step = 100)
        ),

        hr(),

        downloadButton(ns(paste0("plotly_download_", plot_id)), "Download",
                      class = "btn-primary btn-block"),

        br(), br(),

        div(
          style = "font-size: 11px; color: #666;",
          icon("info-circle"),
          " Tip: Use the camera icon in the plot toolbar for quick PNG export."
        )
      ),

      footer = modalButton("Close")
    ))
  })

  # Download handler
  output[[paste0("plotly_download_", plot_id)]] <- downloadHandler(
    filename = function() {
      format <- input[[paste0("plotly_format_", plot_id)]]
      ext <- switch(format, "png" = "png", "html" = "html", "csv" = "csv", "png")
      paste0(filename_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      format <- input[[paste0("plotly_format_", plot_id)]]

      if (format == "csv" && !is.null(data_fn)) {
        data <- tryCatch(data_fn(), error = function(e) NULL)
        if (!is.null(data)) {
          write.csv(data, file, row.names = FALSE)
        }
      } else if (format == "html") {
        p <- tryCatch(plotly_fn(), error = function(e) NULL)
        if (!is.null(p)) {
          htmlwidgets::saveWidget(as_widget(p), file, selfcontained = TRUE)
        }
      } else {
        # For PNG, we'd need kaleido or orca - show notification
        p <- tryCatch(plotly_fn(), error = function(e) NULL)
        if (!is.null(p)) {
          # Try kaleido first
          if (requireNamespace("kaleido", quietly = TRUE)) {
            kaleido::save_image(p, file,
                               width = input[[paste0("plotly_width_", plot_id)]],
                               height = input[[paste0("plotly_height_", plot_id)]])
          } else {
            # Fallback to HTML
            htmlwidgets::saveWidget(as_widget(p), sub("\\.png$", ".html", file), selfcontained = TRUE)
            showNotification("Saved as HTML. Install 'kaleido' for PNG export.", type = "warning")
          }
        }
      }
    }
  )
}
