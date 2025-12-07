# =============================================================================
# Error Handler Functions - helper/error_handlers.R
# =============================================================================

# Main error handler for strain analysis modules
strain_error_handler <- function(error, session, output, context = "", additional_info = list()) {
  error_msg <- paste("Error in", context, ":", error$message)
  
  # Log to console
  cat("\n", rep("=", 80), "\n")
  cat("STRAIN ANALYSIS ERROR:\n")
  cat("Context:", context, "\n")
  cat("Error:", error$message, "\n")
  cat("Additional Info:", paste(names(additional_info), additional_info, sep = "=", collapse = ", "), "\n")
  cat(rep("=", 80), "\n\n")
  
  # Show notification to user
  showNotification(
    paste("Error:", error$message), 
    type = "error",
    duration = 10
  )
  
  # Try to show error in plots if output exists
  tryCatch({
    if (!is.null(output) && exists("plotly_empty")) {
      # Return empty plotly if it's a plot error
      return(plotly_empty())
    }
  }, error = function(e) {
    # Ignore if plotly_empty doesn't exist
  })
  
  return(NULL)
}

# Empty plotly plot for error cases
plotly_empty <- function(message = "Error generating plot") {
  plot_ly() %>%
    add_text(
      text = message,
      x = 0.5, 
      y = 0.5,
      textfont = list(size = 16, color = "red")
    ) %>%
    layout(
      showlegend = FALSE,
      xaxis = list(showgrid = FALSE, showticklabels = FALSE, title = ""),
      yaxis = list(showgrid = FALSE, showticklabels = FALSE, title = "")
    )
}

# NOTE: Plot styling and UI control functions are defined in helper/plot_customization.R
# This file contains only error handling functions to avoid function name collisions
