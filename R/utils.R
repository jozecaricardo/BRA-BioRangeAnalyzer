#' Utility Functions for BioGeoBEARS Shiny App
#'
#' Collection of helper functions for data processing and validation

#' Calculate Likelihood Ratio Test
#'
#' Performs likelihood ratio test between two nested models
#'
#' @param res_simple Results from simpler (nested) model
#' @param res_complex Results from more complex model
#'
#' @return List with LRT statistics
#'
#' @keywords internal
calc_lrt <- function(res_simple, res_complex) {
  lnl_simple <- res_simple$total_loglikelihood
  lnl_complex <- res_complex$total_loglikelihood

  # Calculate degrees of freedom
  df_simple <- length(res_simple$inputs$BioGeoBEARS_model_object@params_table[
    res_simple$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"])
  df_complex <- length(res_complex$inputs$BioGeoBEARS_model_object@params_table[
    res_complex$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"])

  df <- df_complex - df_simple

  # Test statistic
  chisq <- 2 * (lnl_complex - lnl_simple)

  # P-value
  pval <- stats::pchisq(chisq, df = df, lower.tail = FALSE)

  list(
    model_simple = res_simple$inputs$description,
    model_complex = res_complex$inputs$description,
    lnl_simple = lnl_simple,
    lnl_complex = lnl_complex,
    df = df,
    chisq = chisq,
    pval = pval,
    significant = pval < 0.05
  )
}

#' Calculate AICc
#'
#' Calculate Akaike Information Criterion corrected for small samples
#'
#' @param lnl Log-likelihood value
#' @param k Number of parameters
#' @param n Sample size
#'
#' @return AICc value
#'
#' @keywords internal
calc_aicc <- function(lnl, k, n) {
  aic <- 2 * k - 2 * lnl
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  aicc
}

#' Validate Occurrence Data
#'
#' Check if occurrence data has required columns and valid format
#'
#' @param data Data frame with occurrence data
#'
#' @return Logical, TRUE if valid
#'
#' @keywords internal
validate_occurrence_data <- function(data) {
  required_cols <- c("species", "longitude", "latitude")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: species, longitude, latitude")
  }

  if (nrow(data) == 0) {
    stop("Data is empty")
  }

  if (any(is.na(data$longitude)) || any(is.na(data$latitude))) {
    stop("Coordinates contain NA values")
  }

  TRUE
}

#' Validate Tree File
#'
#' Check if tree file is valid Newick format
#'
#' @param tree_file Path to tree file
#'
#' @return Logical, TRUE if valid
#'
#' @keywords internal
validate_tree_file <- function(tree_file) {
  tryCatch({
    tree <- ape::read.tree(tree_file)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Format Model Results
#'
#' Format BioGeoBEARS results for display
#'
#' @param results List of BioGeoBEARS results
#'
#' @return Data frame with formatted results
#'
#' @keywords internal
format_results <- function(results) {
  if (is.null(results) || length(results) == 0) {
    return(data.frame())
  }

  # Extract key statistics for each model
  result_list <- lapply(results, function(res) {
    data.frame(
      Model = res$inputs$description,
      LnL = round(res$total_loglikelihood, 3),
      d = round(res$inputs$BioGeoBEARS_model_object@params_table["d", "est"], 4),
      e = round(res$inputs$BioGeoBEARS_model_object@params_table["e", "est"], 4),
      j = if ("j" %in% rownames(res$inputs$BioGeoBEARS_model_object@params_table)) {
        round(res$inputs$BioGeoBEARS_model_object@params_table["j", "est"], 4)
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  })

  results_df <- do.call(rbind, result_list)
  rownames(results_df) <- NULL
  results_df
}
