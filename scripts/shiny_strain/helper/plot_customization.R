# Load required libraries
library(ggplot2)
library(paletteer)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(bsicons) # Added for UI icons
library(shiny) # Added for modalDialog etc in error handler

# 1. Theme and Palette Definitions
# -----------------------------------

# List of available ggplot themes
theme_options <- list(
  "Classic" = theme_classic(),
  "Black & White" = theme_bw(),
  "Gray" = theme_gray(),
  "Dark" = theme_dark(),
  "Minimal" = theme_minimal(),
  "Void" = theme_void(),
  "Cowplot" = cowplot::theme_cowplot(),
  "Wall Street Journal" = ggthemes::theme_wsj(),
  "Economist" = ggthemes::theme_economist()
)

# List of available color palettes (categorical and gradient)
available_palettes <- list(
  categorical = list(
    "ggsci_npg" = "ggsci::nrc_npg",
    "ggsci_aaas" = "ggsci::default_aaas",
    "ggsci_nejm" = "ggsci::default_nejm",
    "ggsci_lancet" = "ggsci::lanonc_lancet",
    "ggsci_jama" = "ggsci::default_jama",
    "RColorBrewer_Spectral" = "RColorBrewer::Spectral",
    "RColorBrewer_Set3" = "RColorBrewer::Set3",
    "RColorBrewer_Paired" = "RColorBrewer::Paired",
    "ggsci_category10" = "ggsci::category10_d3",
    "ggthemes_Classic10" = "ggthemes::Classic_10_Light"
  ),
  gradient = list(
    "viridis" = viridis::viridis(256),
    "plasma" = viridis::plasma(256),
    "inferno" = viridis::inferno(256),
    "magma" = viridis::magma(256),
    "cividis" = viridis::cividis(256),
    "blues" = RColorBrewer::brewer.pal(9, "Blues"),
    "reds" = RColorBrewer::brewer.pal(9, "Reds"),
    "greens" = RColorBrewer::brewer.pal(9, "Greens")
  )
)

# 2. Core Functions (Apply plot style and get palette)
# --------------------------------------------------

#' Function to return color vector for palette name
#' @param palette_name Palette name
#' @param palette_type "categorical" or "gradient"
#' @return Color vector
get_palette <- function(palette_name, palette_type = "categorical") {
  if (palette_type == "categorical") {
    palette_path <- available_palettes$categorical[[palette_name]]
    if (is.null(palette_path)) {
      palette_path <- available_palettes$categorical[["ggsci_npg"]]
    }
    return(paletteer::paletteer_d(palette_path))
  } else if (palette_type == "gradient") {
    gradient_colors <- available_palettes$gradient[[palette_name]]
    if (is.null(gradient_colors)) {
      gradient_colors <- available_palettes$gradient[["viridis"]]
    }
    return(gradient_colors)
  }
}

#' Function to apply custom style to ggplot object
#' @param p ggplot object
#' @param input Shiny input object
#' @param plot_number (optional) Number to identify specific plot when there are multiple plots
#' @return ggplot object with style applied
apply_plot_style <- function(p, input, plot_number = NULL) {
  # Apply selected theme
  if (!is.null(input$plot_theme)) {
    p <- p + theme_options[[input$plot_theme]]
  }
  
  # Common style customization
  p <- p + theme(
    # Axis settings
    axis.text.x = element_text(
      size = input$xtext_size,
      angle = input$xtext_angle,
      hjust = ifelse(input$xtext_angle == 0, 0.5, 1),
      family = input$font_family
    ),
    axis.text.y = element_text(size = input$ytext_size, family = input$font_family),
    axis.title.x = element_text(size = input$xtitle_size, family = input$font_family),
    axis.title.y = element_text(size = input$ytitle_size, family = input$font_family),
    
    # Legend settings
    legend.title = element_text(size = input$legend_title_size, family = input$font_family),
    legend.text = element_text(size = input$legend_text_size, family = input$font_family),
    legend.position = input$legend_position,
    
    # Plot margin
    plot.margin = unit(rep(input$plot_margin, 4), "cm")
  )

  # Facet plot strip text (title)
  if (!is.null(input$strip_text)) {
    p <- p + theme(
      strip.text = element_text(size = input$strip_text, family = input$font_family)
    )
  }
  
  # Whether to show grid lines
  if (!input$show_grid) {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }
  
  # Flip coordinates
  if (!is.null(input$coord_flip) && input$coord_flip) {
    p <- p + coord_flip()
  }

  # Axis title (per plot or common)
  if (!is.null(plot_number) && !is.null(input[[paste0("x_axis_title", plot_number)]])) {
    p <- p + labs(
      x = input[[paste0("x_axis_title", plot_number)]],
      y = input[[paste0("y_axis_title", plot_number)]]
    )
  } else { 
    p <- p + labs(
      x = input$x_axis_title,
      y = input$y_axis_title
    )
  }

  # Legend title
  if (!is.null(input$legend_title) && input$legend_title != "") {
    p <- p + guides(fill = guide_legend(title = input$legend_title),
                   color = guide_legend(title = input$legend_title))
  }
  
  return(p)
}

#' Function to apply style to Plotly object
apply_plotly_style <- function(p, input) {
  p %>% layout(
    margin = get_plotly_margin(input$plot_margin),
    legend = get_legend_position(input$legend_position),
    width = input$plot_width,
    height = input$plot_height
  )
}


# 3. UI Creation Functions
# -----------------------------------

create_text_controls <- function(ns) {
  accordion_panel(
    "Text & Axis Settings",
    icon = bsicons::bs_icon("fonts"),
    textInput(inputId = ns("x_axis_title"), label = "X-axis Title", value = ""),
    textInput(inputId = ns("y_axis_title"), label = "Y-axis Title", value = ""),
    sliderInput(inputId = ns("xtext_size"), label = "X-axis Text Size", min = 6, max = 20, value = 10),
    sliderInput(inputId = ns("xtext_angle"), label = "X-axis Text Angle", min = 0, max = 90, value = 45),
    sliderInput(inputId = ns("ytext_size"), label = "Y-axis Text Size", min = 6, max = 20, value = 10),
    sliderInput(inputId = ns("xtitle_size"), label = "X-axis Title Size", min = 6, max = 20, value = 12),
    sliderInput(inputId = ns("ytitle_size"), label = "Y-axis Title Size", min = 6, max = 20, value = 12)
  )
}

create_legend_controls <- function(ns) {
  accordion_panel(
    "Legend Settings",
    icon = bsicons::bs_icon("list-ul"),
    textInput(inputId = ns("legend_title"), label = "Legend Title", value = ""),
    sliderInput(inputId = ns("legend_title_size"), label = "Legend Title Size", min = 6, max = 30, value = 10.5),
    sliderInput(inputId = ns("legend_text_size"), label = "Legend Text Size", min = 6, max = 30, value = 8.5),
    selectInput(inputId = ns("legend_position"), label = "Legend Position", choices = c("top", "bottom", "left", "right", "none"), selected = "bottom")
  )
}

create_plot_size_controls <- function(ns) {
  accordion_panel(
    "Plot Dimensions",
    icon = bsicons::bs_icon("aspect-ratio"),
    sliderInput(ns("plot_width"), "Width (px)", value = 800, min = 400, max = 1600, step = 50),
    sliderInput(ns("plot_height"), "Height (px)", value = 600, min = 300, max = 1200, step = 50),
    sliderInput(ns("plot_margin"), "Plot Margin (cm)", min = 0, max = 2, value = 0.5, step = 0.1)
  )
}

create_theme_controls <- function(ns) {
  accordion_panel(
    "Theme & Style",
    icon = bsicons::bs_icon("palette"),
    selectInput(inputId = ns("plot_theme"), label = "Plot Theme", choices = names(theme_options), selected = "Black & White"),
    selectInput(inputId = ns("color_palette"), label = "Color Palette", choices = names(available_palettes$categorical), selected = "ggsci_npg"),
    checkboxInput(inputId = ns("show_grid"), label = "Show Grid Lines", value = TRUE),
    checkboxInput(inputId = ns("coord_flip"), label = "Flip Coordinates", value = FALSE),
    selectInput(inputId = ns("font_family"), label = "Font Family", choices = c("Arial", "Times New Roman", "Courier", "Helvetica", "Georgia"), selected = "Arial")
  )
}

create_download_controls <- function(ns) {
  accordion_panel(
    "Download Options",
    icon = bsicons::bs_icon("download"),
    selectInput(ns("file_format"), "File Format", choices = c("PNG", "JPEG", "TIFF", "PDF"), selected = "PNG"),
    sliderInput(ns("download_width"), "Download Width (inches)", min = 4, max = 20, value = 10, step = 0.5),
    sliderInput(ns("download_height"), "Download Height (inches)", min = 3, max = 15, value = 8, step = 0.5),
    sliderInput(ns("download_dpi"), "Resolution (DPI)", min = 72, max = 600, value = 300, step = 24),
    downloadButton(ns("download_plot"), "Download Plot", class = "btn-primary btn-block")
  )
}

#' Create a download button for individual plots
#' @param ns Namespace function
#' @param plot_id ID of the plot
#' @param btn_class Optional CSS class for the button
#' @return Action button for opening download modal
create_plot_download_btn <- function(ns, plot_id, btn_class = "btn-sm btn-outline-primary") {
  actionButton(
    ns(paste0("download_btn_", plot_id)),
    label = tagList(bsicons::bs_icon("download"), "Download"),
    class = btn_class,
    style = "margin-left: 10px;"
  )
}

#' Create download modal dialog content
#' @param ns Namespace function
#' @param plot_id ID of the plot
#' @return Modal dialog UI
create_download_modal <- function(ns, plot_id) {
  modalDialog(
    title = tagList(bsicons::bs_icon("download"), " Download Plot"),
    size = "l",
    easyClose = TRUE,

    fluidRow(
      column(4,
        wellPanel(
          style = "background: #f8f9fa; padding: 15px;",
          h5("Download Settings", style = "margin-top: 0;"),

          selectInput(
            ns(paste0("modal_format_", plot_id)),
            "File Format:",
            choices = c("PNG" = "png", "JPEG" = "jpeg", "PDF" = "pdf", "SVG" = "svg"),
            selected = "png"
          ),

          numericInput(
            ns(paste0("modal_width_", plot_id)),
            "Width (inches):",
            value = 10,
            min = 4,
            max = 20,
            step = 0.5
          ),

          numericInput(
            ns(paste0("modal_height_", plot_id)),
            "Height (inches):",
            value = 8,
            min = 3,
            max = 15,
            step = 0.5
          ),

          numericInput(
            ns(paste0("modal_dpi_", plot_id)),
            "Resolution (DPI):",
            value = 300,
            min = 72,
            max = 600,
            step = 24
          ),

          hr(),

          actionButton(
            ns(paste0("preview_btn_", plot_id)),
            label = tagList(bsicons::bs_icon("eye"), " Update Preview"),
            class = "btn-info btn-block",
            style = "margin-bottom: 10px;"
          ),

          downloadButton(
            ns(paste0("modal_download_", plot_id)),
            label = tagList(bsicons::bs_icon("download"), " Download"),
            class = "btn-primary btn-block"
          )
        )
      ),

      column(8,
        div(
          style = "border: 1px solid #ddd; border-radius: 4px; padding: 10px; background: white; min-height: 400px;",
          h5("Preview", style = "margin-top: 0; color: #666;"),
          div(
            style = "overflow: auto; max-height: 500px;",
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

#' Server-side handler for download modal
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param plot_id ID of the plot
#' @param plot_reactive Reactive expression returning the ggplot object
#' @param filename_prefix Prefix for the downloaded file name
setup_download_modal_handler <- function(input, output, session, plot_id, plot_reactive, filename_prefix = "plot") {
  ns <- session$ns

  # Open modal when download button clicked
  observeEvent(input[[paste0("download_btn_", plot_id)]], {
    showModal(create_download_modal(ns, plot_id))

    # Trigger initial preview
    shinyjs::delay(100, {
      shinyjs::click(paste0("preview_btn_", plot_id))
    })
  })

  # Render preview
  output[[paste0("modal_preview_", plot_id)]] <- renderPlot({
    # React to preview button
    input[[paste0("preview_btn_", plot_id)]]

    p <- plot_reactive()
    req(p)

    p
  }, width = function() {
    w <- input[[paste0("modal_width_", plot_id)]]
    if (is.null(w)) w <- 10
    w * 50  # Convert inches to approximate pixels for preview
  }, height = function() {
    h <- input[[paste0("modal_height_", plot_id)]]
    if (is.null(h)) h <- 8
    h * 50
  })

  # Download handler
  output[[paste0("modal_download_", plot_id)]] <- downloadHandler(
    filename = function() {
      format <- input[[paste0("modal_format_", plot_id)]]
      if (is.null(format)) format <- "png"
      paste0(filename_prefix, "_", Sys.Date(), ".", format)
    },
    content = function(file) {
      p <- plot_reactive()
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
}


# 4. Helper Functions (Plotly, Error Handling, State Storage)
# --------------------------------------------------

# Convert Plotly legend position
get_legend_position <- function(position) {
  switch(position,
    "top" = list(orientation = "h", x = 0.5, y = 1.1, xanchor = "center", yanchor = "bottom"),
    "bottom" = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center", yanchor = "top"),
    "left" = list(orientation = "v", x = -0.2, y = 0.5, xanchor = "right", yanchor = "middle"),
    "right" = list(orientation = "v", x = 1.1, y = 0.5, xanchor = "left", yanchor = "middle"),
    "none" = list(visible = FALSE),
    list()
  )
}

# Convert Plotly margin
get_plotly_margin <- function(margin_cm) {
  margin_px <- margin_cm * 37.8  # Approximate cm to px conversion
  list(
    l = margin_px,
    r = margin_px,
    t = margin_px + 40,  # Extra space for title
    b = margin_px + 40,  # Extra space for axis labels
    pad = 4
  )
}

#' Global error handler
global_error_handler <- function(e, session, output, code_snippet = NULL, input_list = NULL) {

  # Display simple error message only
  error_message <- paste0("Error: ", e$message)
  
  showModal(modalDialog(
    title = "An Error Occurred",
    size = "m",
    easyClose = TRUE,
    div(
      style = "color: red; font-weight: bold; margin-bottom: 15px;",
      error_message
    ),
    div(
      style = "color: #666; font-size: 0.9em;",
      "Please check your data and try again."
    ),
    footer = modalButton("Close")
  ))
  
  message("Error handled: ", e$message)
}


#' Function to return current plot style settings as a list
get_plot_style_state <- function(input) {
  list(
    plot_theme = input$plot_theme,
    show_grid = input$show_grid,
    coord_flip = input$coord_flip,
    font_family = input$font_family,
    color_palette = input$color_palette,
    xtext_size = input$xtext_size,
    ytext_size = input$ytext_size,
    xtitle_size = input$xtitle_size,
    ytitle_size = input$ytitle_size,
    xtext_angle = input$xtext_angle,
    legend_position = input$legend_position,
    legend_title_size = input$legend_title_size,
    legend_text_size = input$legend_text_size,
    plot_width = input$plot_width,
    plot_height = input$plot_height,
    plot_margin = input$plot_margin,
    x_axis_title = input$x_axis_title,
    y_axis_title = input$y_axis_title,
    legend_title = input$legend_title
  )
}