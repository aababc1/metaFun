# modules/nucleotide_diversity_stats_helpers.R

#' @title Statistical Helper Functions for Nucleotide Diversity Analysis
#' @description
#' This script contains a set of helper functions for performing statistical
#' tests commonly used in microbiome analysis, specifically for nucleotide
#' diversity. These functions are designed to be used within a Shiny application,
#' providing support for group comparisons and significance testing.
#'
#' The functions included are:
#' - apply_bh_correction: Conditionally applies Benjamini-Hochberg correction.
#' - add_significance_notation: Adds significance asterisks based on p-values.
#' - perform_wilcoxon_test: Performs a two-group comparison (Wilcoxon rank-sum).
#' - perform_kruskal_wallis_test: Performs a multi-group comparison (Kruskal-Wallis).
#' - perform_dunn_posthoc: Performs Dunn's post-hoc test for pairwise comparisons.
#' - perform_spearman_correlation: Performs Spearman correlation for numerical variables.
#'
#' Dependencies: dplyr, tibble, FSA. These packages must be installed.

#' Apply Benjamini-Hochberg correction if multiple comparisons are made.
#'
#' This function takes a vector of p-values and applies the Benjamini-Hochberg
#' (BH) correction for multiple testing. The correction is only applied if more
#' than one p-value is provided, aligning with the principle that correction is
#' for multiple comparisons.
#'
#' @param p_values A numeric vector of p-values.
#' @return A numeric vector of p-values. If the input vector has a length
#'   greater than 1, the returned p-values are adjusted using the BH method.
#'   Otherwise, the original p-values are returned.
apply_bh_correction <- function(p_values) {
  if (!is.numeric(p_values)) {
    stop("Input must be a numeric vector of p-values.")
  }
  # Only apply correction if there are multiple p-values to correct
  if (length(p_values) > 1) {
    return(stats::p.adjust(p_values, method = "BH"))
  } else {
    return(p_values)
  }
}

#' Add significance notation to a results data frame.
#'
#' This function adds a 'significance' column to a data frame containing
#' p-values. The notation is based on standard p-value thresholds.
#'
#' @param data A data frame.
#' @param p_col A string specifying the name of the column that contains the
#'   p-values (or adjusted p-values). Defaults to "p.adj".
#' @return The input data frame with an additional 'significance' column
#'   containing notations ('***', '**', '*', 'ns'). NA p-values result in NA.
add_significance_notation <- function(data, p_col = "p.adj") {
  if (!p_col %in% names(data)) {
    stop(paste0("Column '", p_col, "' not found in the data frame."))
  }

  p_values <- data[[p_col]]

  significance_levels <- cut(
    p_values,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns"),
    right = TRUE
  )

  data$significance <- as.character(significance_levels)

  return(data)
}

#' Perform Wilcoxon Rank-Sum test for two-group comparisons.
#'
#' This function performs a Wilcoxon rank-sum (Mann-Whitney U) test and
#' calculates summary statistics (n, mean, median) for each group. It does not
#' compute effect sizes.
#'
#' @param data A data frame containing the data for the test.
#' @param formula A formula of the form `y ~ x`, where `y` is a numeric
#'   variable and `x` is a factor with exactly two levels.
#' @return A tibble with one row containing the test method, statistic,
#'   p-value, and a nested tibble with summary statistics for each group.
perform_wilcoxon_test <- function(data, formula) {
  vars <- all.vars(formula)
  outcome_var <- vars[1]
  group_var <- vars[2]

  if (!all(c(outcome_var, group_var) %in% names(data))) {
    stop("Variables specified in the formula are not present in the data frame.")
  }

  # Ensure the grouping variable is a factor
  data[[group_var]] <- as.factor(data[[group_var]])

  if (nlevels(data[[group_var]]) != 2) {
    stop("Wilcoxon test requires the grouping variable to have exactly two levels.")
  }

  # Perform the test with error handling for insufficient data
  test_result <- tryCatch({
    stats::wilcox.test(formula, data = data)
  }, error = function(e) {
    warning("Wilcoxon test failed. This can happen with insufficient data. Error: ", e$message)
    list(method = "Wilcoxon rank sum test", statistic = NA_real_, p.value = NA_real_)
  })

  # Calculate summary statistics
  summary_stats <- tryCatch({
    data %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean = mean(.data[[outcome_var]], na.rm = TRUE),
        median = median(.data[[outcome_var]], na.rm = TRUE),
        .groups = 'drop'
      )
  }, error = function(e) {
    warning("Could not compute summary statistics. Error: ", e$message)
    NULL
  })

  # Wilcoxon test: p.adj = p.value (no correction for single comparison)
  result_df <- tibble::tibble(
    method = test_result$method,
    statistic.W = test_result$statistic,
    p.value = test_result$p.value,
    p.adj = test_result$p.value,  # Raw p-value, no BH correction
    summary = list(summary_stats)
  )

  return(result_df)
}

#' Perform Kruskal-Wallis test for multi-group comparisons.
#'
#' This function performs a Kruskal-Wallis rank sum test for comparing three
#' or more groups. It also calculates summary statistics (n, mean, median) for
#' each group.
#'
#' @param data A data frame containing the data.
#' @param formula A formula of the form `y ~ x`, where `y` is numeric and `x`
#'   is a factor with three or more levels.
#' @return A tibble with the test method, chi-squared statistic, p-value, and
#'   a nested tibble with summary statistics.
perform_kruskal_wallis_test <- function(data, formula) {
  vars <- all.vars(formula)
  outcome_var <- vars[1]
  group_var <- vars[2]

  if (!all(c(outcome_var, group_var) %in% names(data))) {
    stop("Variables specified in the formula are not present in the data frame.")
  }

  # Ensure the grouping variable is a factor
  data[[group_var]] <- as.factor(data[[group_var]])

  # Perform the test with error handling
  test_result <- tryCatch({
    stats::kruskal.test(formula, data = data)
  }, error = function(e) {
    warning("Kruskal-Wallis test failed. Error: ", e$message)
    list(method = "Kruskal-Wallis rank sum test", statistic = NA_real_, p.value = NA_real_)
  })

  # Calculate summary statistics
  summary_stats <- tryCatch({
    data %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean = mean(.data[[outcome_var]], na.rm = TRUE),
        median = median(.data[[outcome_var]], na.rm = TRUE),
        .groups = 'drop'
      )
  }, error = function(e) {
    warning("Could not compute summary statistics. Error: ", e$message)
    NULL
  })

  # Kruskal-Wallis test: p.adj = p.value (BH correction handled by Dunn's post-hoc)
  result_df <- tibble::tibble(
    method = test_result$method,
    statistic.chi.squared = test_result$statistic,
    p.value = test_result$p.value,
    p.adj = test_result$p.value,  # Raw p-value, BH handled by Dunn's test
    summary = list(summary_stats)
  )

  return(result_df)
}


#' Perform Dunn's post-hoc test for pairwise comparisons.
#'
#' This function is intended to be called after a significant Kruskal-Wallis
#' test. It uses the `dunnTest` function from the `FSA` package to perform
#' pairwise comparisons among groups.
#'
#' @param data A data frame containing the data.
#' @param formula A formula of the form `y ~ x`.
#' @param p.adjust.method The method for p-value adjustment. Defaults to "fdr",
#'   which is equivalent to Benjamini-Hochberg.
#' @return A tidy tibble with pairwise comparison results, including the
#'   Z-statistic, unadjusted p-value, and adjusted p-value.
perform_dunn_posthoc <- function(data, formula, p.adjust.method = "fdr") {
  # dunnTest requires the FSA package
  if (!requireNamespace("FSA", quietly = TRUE)) {
    stop("Package 'FSA' is required. Please install it with install.packages('FSA')")
  }

  # Perform the test with error handling
  dunn_result <- tryCatch({
    FSA::dunnTest(formula, data = data, method = p.adjust.method)
  }, error = function(e) {
    warning("Dunn's test failed. This can occur with groups that lack variance. Error: ", e$message)
    return(NULL)
  })

  if (is.null(dunn_result)) {
    # Return a structured empty tibble on failure
    return(tibble::tibble(
      comparison = NA_character_,
      Z = NA_real_,
      p.unadj = NA_real_,
      p.adj = NA_real_
    ))
  }

  # Extract the results data frame and format it as a tibble
  result_df <- tibble::as_tibble(dunn_result$res)

  # Rename columns for consistency and clarity
  result_df <- result_df %>%
    dplyr::rename(
      comparison = .data$Comparison,
      Z = .data$Z,
      p.unadj = .data$P.unadj,
      p.adj = .data$P.adj
    )

  return(result_df)
}

#' Perform Spearman Correlation test for numerical variables.
#'
#' This function performs a Spearman rank correlation test to assess the
#' relationship between a numerical variable and a continuous outcome (e.g.,
#' nucleotide diversity). Unlike Wilcoxon which tests group differences,
#' Spearman tests for monotonic association.
#'
#' @param data A data frame containing the data for the test.
#' @param formula A formula of the form `y ~ x`, where `y` is the outcome
#'   variable (e.g., nucleotide diversity) and `x` is a numerical predictor.
#' @return A tibble with one row containing the test method, correlation
#'   coefficient (rho), p-value, and sample size.
perform_spearman_correlation <- function(data, formula) {
  vars <- all.vars(formula)
  outcome_var <- vars[1]
  predictor_var <- vars[2]

  if (!all(c(outcome_var, predictor_var) %in% names(data))) {
    stop("Variables specified in the formula are not present in the data frame.")
  }

  # Extract vectors and remove NAs
  x <- data[[predictor_var]]
  y <- data[[outcome_var]]

  # Remove pairs with missing values
  complete_cases <- complete.cases(x, y)
  x <- x[complete_cases]
  y <- y[complete_cases]

  n <- length(x)

  # Ensure both variables are numeric
  x <- as.numeric(as.character(x))
  y <- as.numeric(y)

  if (any(is.na(x)) || any(is.na(y))) {
    warning("Predictor variable contains non-numeric values. Spearman correlation requires numeric data.")
    return(tibble::tibble(
      method = "Spearman's rank correlation",
      statistic.rho = NA_real_,
      p.value = NA_real_,
      p.adj = NA_real_,
      n = 0
    ))
  }

  if (n < 3) {
    warning("Insufficient data for Spearman correlation (n < 3)")
    return(tibble::tibble(
      method = "Spearman's rank correlation",
      statistic.rho = NA_real_,
      p.value = NA_real_,
      p.adj = NA_real_,
      n = n
    ))
  }

  # Perform Spearman correlation test
  test_result <- tryCatch({
    cor.test(x, y, method = "spearman", exact = FALSE)
  }, error = function(e) {
    warning("Spearman correlation test failed. Error: ", e$message)
    list(method = "Spearman's rank correlation rho", estimate = NA_real_, p.value = NA_real_)
  })

  # Spearman correlation: p.adj = p.value (no correction for single test)
  result_df <- tibble::tibble(
    method = "Spearman's rank correlation",
    statistic.rho = as.numeric(test_result$estimate),
    p.value = test_result$p.value,
    p.adj = test_result$p.value,  # Raw p-value, no BH correction for single correlation
    n = n
  )

  return(result_df)
}