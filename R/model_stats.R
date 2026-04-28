#' Calculate AIC and AICc
#'
#' Calculate Akaike Information Criterion and corrected AIC
#'
#' @param lnl Log-likelihood value
#' @param k Number of parameters
#' @param n Sample size (number of areas)
#'
#' @return List with aic and aicc values
#'
#' @export
calc_aic_aicc <- function(lnl, k, n) {
  aic <- 2 * k - 2 * lnl
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  
  list(
    aic = aic,
    aicc = aicc
  )
}

#' Calculate Delta AIC and Delta AICc
#'
#' Calculate differences in AIC and AICc relative to the best model
#'
#' @param aic Vector of AIC values
#' @param aicc Vector of AICc values
#'
#' @return List with delta_aic and delta_aicc vectors
#'
#' @export
calc_delta_aic <- function(aic, aicc) {
  delta_aic <- aic - min(aic, na.rm = TRUE)
  delta_aicc <- aicc - min(aicc, na.rm = TRUE)
  
  list(
    delta_aic = delta_aic,
    delta_aicc = delta_aicc
  )
}

#' Perform Likelihood Ratio Test
#'
#' Test if a complex model is significantly better than a simple model
#'
#' @param lnl_simple Log-likelihood of simple model
#' @param lnl_complex Log-likelihood of complex model
#' @param df Degrees of freedom (difference in number of parameters)
#'
#' @return List with test statistic, p-value, and significance
#'
#' @export
perform_lrt <- function(lnl_simple, lnl_complex, df) {
  # LRT statistic: 2 * (lnl_complex - lnl_simple)
  # Follows chi-square distribution with df degrees of freedom
  
  if (lnl_complex < lnl_simple) {
    warning("Complex model has lower likelihood than simple model")
  }
  
  test_stat <- 2 * (lnl_complex - lnl_simple)
  pval <- stats::pchisq(test_stat, df = df, lower.tail = FALSE)
  
  list(
    test_statistic = test_stat,
    df = df,
    p_value = pval,
    significant = pval < 0.05
  )
}

#' Create Model Comparison Table
#'
#' Create a comprehensive comparison table for multiple models
#'
#' @param models Character vector of model names
#' @param lnl Numeric vector of log-likelihood values
#' @param params List of parameter values for each model
#' @param n_areas Number of areas (for AICc calculation)
#'
#' @return Data frame with model comparison
#'
#' @export
create_model_comparison <- function(models, lnl, params, n_areas) {
  
  # Calculate number of parameters for each model
  n_params <- sapply(params, function(x) {
    n <- 2  # d and e
    if (!is.na(x$j) && x$j > 0) n <- n + 1  # Add j if present
    n
  })
  
  # Calculate AIC and AICc
  aic <- sapply(seq_along(models), function(i) {
    calc_aic_aicc(lnl[i], n_params[i], n_areas)$aic
  })
  
  aicc <- sapply(seq_along(models), function(i) {
    calc_aic_aicc(lnl[i], n_params[i], n_areas)$aicc
  })
  
  # Calculate delta AIC and delta AICc
  delta_aic <- aic - min(aic, na.rm = TRUE)
  delta_aicc <- aicc - min(aicc, na.rm = TRUE)
  
  # Create comparison table
  comparison <- data.frame(
    Model = models,
    LnL = round(lnl, 3),
    K = n_params,
    AIC = round(aic, 2),
    DeltaAIC = round(delta_aic, 2),
    AICc = round(aicc, 2),
    DeltaAICc = round(delta_aicc, 2),
    stringsAsFactors = FALSE
  )
  
  # Add parameter values
  comparison$d <- round(sapply(params, function(x) x$d), 4)
  comparison$e <- round(sapply(params, function(x) x$e), 4)
  comparison$j <- sapply(params, function(x) {
    if (is.na(x$j)) "-" else round(x$j, 4)
  })
  
  # Sort by AICc
  comparison <- comparison[order(comparison$AICc), ]
  rownames(comparison) <- NULL
  
  comparison
}

#' Calculate Akaike Weights
#'
#' Calculate Akaike weights for model comparison
#'
#' @param aicc Vector of AICc values
#'
#' @return Vector of Akaike weights
#'
#' @export
calc_akaike_weights <- function(aicc) {
  delta_aicc <- aicc - min(aicc, na.rm = TRUE)
  L <- exp(-0.5 * delta_aicc)
  w <- L / sum(L, na.rm = TRUE)
  w
}

#' Generate Mock BioGeoBEARS Results
#'
#' Generate realistic mock results for demonstration
#'
#' @param models Character vector of model names
#' @param n_areas Number of areas
#' @param d_bounds Numeric vector with min and max for d parameter
#' @param e_bounds Numeric vector with min and max for e parameter
#' @param j_bounds Numeric vector with min and max for j parameter
#'
#' @return List with lnl and params
#'
#' @export
generate_mock_results <- function(models, n_areas, d_bounds, e_bounds, j_bounds) {
  
  lnl <- numeric(length(models))
  params <- list()
  
  for (i in seq_along(models)) {
    # Generate realistic log-likelihood values
    # Better models have higher (less negative) lnl
    base_lnl <- -50 - (i - 1) * 5
    lnl[i] <- base_lnl + rnorm(1, 0, 2)
    
    # Generate parameter values
    d <- runif(1, d_bounds[1], d_bounds[2])
    e <- runif(1, e_bounds[1], e_bounds[2])
    j <- NA
    
    if (grepl("\\+J", models[i])) {
      j <- runif(1, j_bounds[1], j_bounds[2])
    }
    
    params[[i]] <- list(d = d, e = e, j = j)
  }
  
  list(
    lnl = lnl,
    params = params
  )
}
