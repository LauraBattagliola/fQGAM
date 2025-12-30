#' @title Utility Functions for Functional Quantile GAM
#' 
#' @description 
#' Helper functions for data preparation, prediction, bias adjustment, 
#' and confidence interval construction in functional quantile regression.
#' 
#' @name fqgam-utilities
NULL


#' Prepare Data for Functional Quantile GAM
#'
#' Helper function to prepare longitudinal functional data for use with
#' \code{\link[qgam]{qgam}}.
#'
#' @details
#' This function creates the required matrices for the tensor product smooth 
#' representation of the functional coefficient. The model uses:
#' \itemize{
#'   \item \code{sMat}: Matrix of functional domain coordinates, where each row 
#'     corresponds to an observation and contains the functional grid points
#'   \item \code{tMat}: Matrix of longitudinal time coordinates, where each row 
#'     contains the same time value repeated S times
#' }
#' 
#' The functional covariate enters the model formula as:
#' \code{te(sMat, tMat, by = Xobs, bs = c("cr", "cr"), k = c(k_s, k_t))}
#' 
#' where the \code{by} argument allows for the integration approximation.
#'
#' @param y Numeric vector. Response variable.
#' @param id Factor or character vector. Subject identifiers.
#' @param time Numeric vector. Longitudinal time stamps.
#' @param Xobs Matrix. Functional covariate observations where rows are 
#'   observations and columns are functional grid points.
#' @param s_range Numeric vector of length 2. Range of the functional domain.
#'   Default is \code{c(0, 1)}.
#'
#' @return A data.frame suitable for use with \code{qgam}, containing:
#'   \describe{
#'     \item{y}{Response variable}
#'     \item{id}{Subject identifier (as factor)}
#'     \item{time}{Longitudinal time}
#'     \item{Xobs}{Functional covariate matrix}
#'     \item{sMat}{Matrix of functional domain coordinates}
#'     \item{tMat}{Matrix of longitudinal time coordinates}
#'   }
#'
#' @examples
#' \dontrun{
#' library(refund)
#' data(DTI)
#' myData <- DTI[complete.cases(DTI$cca), ]
#' myData <- subset(myData, Nscans > 1)
#' 
#' myDataFit <- prep_fqgam_data(
#'   y = myData$pasat,
#'   id = myData$ID,
#'   time = myData$visit.time,
#'   Xobs = myData$cca
#' )
#' }
#'
#' @export
prep_fqgam_data <- function(y, id, time, Xobs, s_range = c(0, 1)) {
  
  if (!is.matrix(Xobs)) {
    Xobs <- as.matrix(Xobs)
  }
  
  N <- nrow(Xobs)
  S <- ncol(Xobs)
  
  if (length(y) != N) {
    stop("Length of 'y' must match number of rows in 'Xobs'")
  }
  if (length(id) != N) {
    stop("Length of 'id' must match number of rows in 'Xobs'")
  }
  if (length(time) != N) {
    stop("Length of 'time' must match number of rows in 'Xobs'")
  }
  
  sVec <- seq(s_range[1], s_range[2], length = S)
  sGrid <- matrix(sVec, N, S, byrow = TRUE)
  tGrid <- matrix(time, N, S, byrow = FALSE)
  
  d <- data.frame(
    y = y,
    id = factor(id),
    time = time
  )
  d$Xobs <- Xobs
  d$sMat <- sGrid
  d$tMat <- tGrid
  
  return(d)
}


#' Compute Predictions from Functional Quantile GAM
#'
#' Computes predicted quantiles at specified reference curves and 
#' longitudinal time points from a fitted \code{qgam} model.
#'
#' @details
#' The predicted quantile at reference curve \eqn{X(s)} and time \eqn{t} is:
#' \deqn{\hat{Q}_{X,0}^\tau(t) = \hat{\alpha}^\tau(t) + \int_S X(s) \hat{\beta}^\tau(s,t) ds}
#' 
#' The integral is approximated as a Riemann sum. Note that random effects are 
#' excluded from predictions (setting \eqn{u_i = 0}), so these represent the 
#' quantile for a "typical" subject.
#'
#' @param model A fitted \code{qgam} model object.
#' @param Xref Matrix or list of numeric vectors. Reference curves at which 
#'   to compute predictions. If a list, each element should be a numeric 
#'   vector of the same length as the functional grid.
#' @param t_range Integer vector. Longitudinal time points at which to predict.
#'   Default is \code{NULL}, which uses \code{1:max(time)} from the model data.
#' @param exclude Character. Terms to exclude from predictions. 
#'   Default is \code{"s(id)"} to exclude random effects.
#' @param s_range Numeric vector of length 2. Range of the functional domain.
#'   Default is \code{c(0, 1)}.
#'
#' @return A list of class \code{"fqgam_pred"} containing:
#'   \describe{
#'     \item{predictions}{Named list of prediction vectors, one per reference curve}
#'     \item{time}{Vector of time points}
#'     \item{Xref}{The reference curves used}
#'   }
#'
#' @export
predict_fqgam <- function(model, Xref, t_range = NULL, exclude = "s(id)",
                          s_range = c(0, 1)) {
  
  # Handle matrix or list input
  if (is.matrix(Xref)) {
    Xref_list <- lapply(1:nrow(Xref), function(i) Xref[i, ])
    names(Xref_list) <- paste0("ref", 1:nrow(Xref))
  } else if (is.list(Xref)) {
    Xref_list <- Xref
    if (is.null(names(Xref_list))) {
      names(Xref_list) <- paste0("ref", seq_along(Xref_list))
    }
  } else {
    # Single vector
    Xref_list <- list(ref1 = Xref)
  }
  
  S <- length(Xref_list[[1]])
  sVec <- seq(s_range[1], s_range[2], length = S)
  
  # Determine time range
  if (is.null(t_range)) {
    model_data <- model$model
    if ("time" %in% names(model_data)) {
      maxT <- max(model_data$time)
      t_range <- 1:maxT
    } else {
      stop("Cannot determine time range. Please provide 't_range'.")
    }
  }
  
  # Compute predictions for each reference curve
  predictions <- lapply(Xref_list, function(X) {
    pred <- rep(NA, length(t_range))
    
    for (i in seq_along(t_range)) {
      t <- t_range[i]
      nd <- data.frame(
        time = rep(t, S),
        tMat = rep(t, S),
        sMat = sVec,
        Xobs = X * S  # Multiply by S for integration approximation
      )
      pred[i] <- mean(stats::predict(model, exclude = exclude, 
                                     newdata = nd, newdata.guaranteed = TRUE))
    }
    
    return(pred)
  })
  
  result <- list(
    predictions = predictions,
    time = t_range,
    Xref = Xref_list
  )
  
  class(result) <- c("fqgam_pred", "list")
  
  return(result)
}


#' Compute Bias Adjustment from Wild Bootstrap
#'
#' Computes bias-adjusted estimates from wild bootstrap results following 
#' the methodology of Battagliola et al. (2022, 2024).
#'
#' @details
#' The bias is estimated as:
#' \deqn{bias_{boot}(\hat{\theta}) = \frac{1}{B} \sum_{b=1}^{B} (\check{\theta}_b - \hat{\theta})}
#' 
#' where \eqn{\check{\theta}_b} are bootstrap estimates and \eqn{\hat{\theta}} is 
#' the original estimate. The bias-adjusted estimator is:
#' \deqn{\hat{\theta}_{adj} = \hat{\theta} - bias_{boot}(\hat{\theta}) = 2\hat{\theta} - \bar{\theta}^*}
#'
#' @param wild_boot A wild bootstrap result object from \code{\link{boot_pred_wild}}.
#' @param original_pred Numeric vector. Original (non-bootstrap) predictions
#'   for the quantity of interest (typically predDiff = pred20 - pred80).
#'
#' @return A list of class \code{"fqgam_bias_adj"} containing:
#'   \describe{
#'     \item{adjusted}{Bias-adjusted estimates}
#'     \item{original}{Original estimates}
#'     \item{bias}{Estimated bias at each time point}
#'     \item{boot_mean}{Mean of bootstrap estimates}
#'   }
#'
#' @references
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2022).
#' A bias-adjusted estimator in quantile regression for clustered data.
#' Econometrics and Statistics, 23, 165-186.
#'
#' @export
compute_bias_adjustment <- function(wild_boot, original_pred) {
  
  if (!inherits(wild_boot, "fqgam_boot")) {
    stop("'wild_boot' must be an object of class 'fqgam_boot'")
  }
  
  if (!is.null(wild_boot$bootstrap_type) && wild_boot$bootstrap_type != "wild") {
    warning("Bias adjustment is typically computed from wild bootstrap results. ",
            "Block bootstrap was provided.")
  }
  
  maxT <- ncol(wild_boot$predDiff)
  
  if (length(original_pred) != maxT) {
    stop("Length of 'original_pred' must match number of time points in bootstrap")
  }
  
  # Compute difference: bootstrap estimate - original estimate
  diff_adj <- sweep(wild_boot$predDiff, 2, original_pred, "-")
  
  # Bias = mean(bootstrap - original)
  bias <- colMeans(diff_adj)
  
  # Bias-adjusted estimate = original - bias = 2*original - mean(bootstrap)
  adjusted <- original_pred - bias
  
  result <- list(
    adjusted = adjusted,
    original = original_pred,
    bias = bias,
    boot_mean = colMeans(wild_boot$predDiff)
  )
  
  class(result) <- c("fqgam_bias_adj", "list")
  
  return(result)
}


#' Compute Bootstrap Confidence Intervals
#'
#' Computes confidence intervals from bootstrap results, combining wild bootstrap 
#' for bias adjustment and block bootstrap for variance estimation as recommended 
#' in Battagliola et al. (2024).
#'
#' @details
#' \strong{Recommended Approach:}
#' The recommended approach combines wild bootstrap for bias adjustment and
#' block bootstrap for standard error estimation. The approximate (1-alpha) 
#' confidence interval is: theta - bias +/- z * sd_boot, where bias is estimated 
#' from wild bootstrap and sd_boot from block bootstrap.
#' 
#' \strong{Alternative Methods:}
#' Three methods are available: "normal" (uses normal approximation with 
#' bootstrap SE), "percentile" (uses quantiles of bootstrap distribution 
#' directly), and "basic" (uses basic bootstrap intervals with reflected 
#' percentiles).
#' 
#' Model-based standard errors from qgam tend to underestimate the true 
#' sampling variation due to penalization of random effects.
#'
#' @param block_boot A block bootstrap result object from \code{\link{boot_pred_block}}.
#'   Used for variance/SE estimation.
#' @param wild_boot Optional. A wild bootstrap result object from \code{\link{boot_pred_wild}}.
#'   Used for bias adjustment. If \code{NULL}, no bias adjustment is performed.
#' @param original_pred Numeric vector. Original (non-bootstrap) predictions.
#' @param model_se Optional. Numeric vector of model-based standard errors.
#' @param alpha Numeric. Significance level for confidence intervals.
#'   Default is 0.05.
#' @param method Character. One of \code{"normal"}, \code{"percentile"}, or \code{"basic"}.
#'
#' @return A list of class \code{"fqgam_ci"} containing:
#'   \describe{
#'     \item{estimate}{Original or bias-adjusted point estimate}
#'     \item{bias}{Estimated bias}
#'     \item{se_boot}{Bootstrap standard error}
#'     \item{se_model}{Model-based standard error (if provided)}
#'     \item{ci_lower}{Lower confidence bound}
#'     \item{ci_upper}{Upper confidence bound}
#'     \item{alpha}{Significance level}
#'     \item{method}{Method used}
#'     \item{bias_adjusted}{Logical indicating if bias adjustment was applied}
#'   }
#'
#' @references
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2024).
#' Quantile regression for longitudinal functional data with application to
#' feed intake of lactating sows. Journal of Agricultural, Biological, and 
#' Environmental Statistics, 30(1), 211-230.
#'
#' @export
boot_ci <- function(block_boot, wild_boot = NULL, original_pred, 
                    model_se = NULL, alpha = 0.05, 
                    method = c("normal", "percentile", "basic")) {
  
  method <- match.arg(method)
  
  if (!inherits(block_boot, "fqgam_boot")) {
    stop("'block_boot' must be an object of class 'fqgam_boot'")
  }
  
  B <- nrow(block_boot$predDiff)
  maxT <- ncol(block_boot$predDiff)
  
  if (length(original_pred) != maxT) {
    stop("Length of 'original_pred' must match number of time points in bootstrap")
  }
  
  # Bootstrap standard error from block bootstrap
  se_boot <- apply(block_boot$predDiff, 2, stats::sd)
  
  # Bias from wild bootstrap (if provided)
  if (!is.null(wild_boot)) {
    if (!inherits(wild_boot, "fqgam_boot")) {
      stop("'wild_boot' must be an object of class 'fqgam_boot'")
    }
    bias_adj <- compute_bias_adjustment(wild_boot, original_pred)
    bias <- bias_adj$bias
    estimate <- bias_adj$adjusted
    bias_adjusted <- TRUE
  } else {
    bias <- rep(0, maxT)
    estimate <- original_pred
    bias_adjusted <- FALSE
  }
  
  # Construct confidence intervals
  z <- stats::qnorm(1 - alpha / 2)
  
  if (method == "normal") {
    ci_lower <- estimate - z * se_boot
    ci_upper <- estimate + z * se_boot
  } else if (method == "percentile") {
    ci_lower <- apply(block_boot$predDiff, 2, stats::quantile, probs = alpha / 2)
    ci_upper <- apply(block_boot$predDiff, 2, stats::quantile, probs = 1 - alpha / 2)
    estimate <- original_pred
    bias_adjusted <- FALSE
  } else if (method == "basic") {
    q_lower <- apply(block_boot$predDiff, 2, stats::quantile, probs = alpha / 2)
    q_upper <- apply(block_boot$predDiff, 2, stats::quantile, probs = 1 - alpha / 2)
    ci_lower <- 2 * original_pred - q_upper
    ci_upper <- 2 * original_pred - q_lower
    estimate <- original_pred
    bias_adjusted <- FALSE
  }
  
  result <- list(
    estimate = estimate,
    bias = bias,
    se_boot = se_boot,
    se_model = model_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    alpha = alpha,
    method = method,
    bias_adjusted = bias_adjusted
  )
  
  class(result) <- c("fqgam_ci", "list")
  
  return(result)
}


#' Compute Model-Based Standard Errors
#'
#' Computes parametric standard errors for prediction differences using 
#' the delta method.
#'
#' @details
#' Model-based SEs typically underestimate the true sampling variation. 
#' Use \code{\link{boot_ci}} with block bootstrap for reliable inference.
#'
#' @param model A fitted \code{qgam} model object.
#' @param X20 Numeric vector. Lower reference curve.
#' @param X80 Numeric vector. Upper reference curve.
#' @param t_range Integer vector. Time points.
#' @param s_range Numeric vector of length 2. Functional domain range.
#'
#' @return A list with predictions and model-based standard errors.
#'
#' @export
compute_model_se <- function(model, X20, X80, t_range = NULL, s_range = c(0, 1)) {
  
  S <- length(X20)
  sVec <- seq(s_range[1], s_range[2], length = S)
  
  if (is.null(t_range)) {
    model_data <- model$model
    if ("time" %in% names(model_data)) {
      maxT <- max(model_data$time)
      t_range <- 1:maxT
    } else {
      stop("Cannot determine time range. Please provide 't_range'.")
    }
  }
  
  Nid <- length(unique(model$model$id))
  Nfixed <- length(model$coefficients) - Nid
  cov_mat <- model$Vp[1:Nfixed, 1:Nfixed]
  
  pred20 <- pred80 <- predDiff <- se_model <- rep(NA, length(t_range))
  
  for (i in seq_along(t_range)) {
    t <- t_range[i]
    
    nd20 <- data.frame(time = rep(t, S), tMat = rep(t, S), 
                       sMat = sVec, Xobs = X20 * S)
    nd80 <- data.frame(time = rep(t, S), tMat = rep(t, S), 
                       sMat = sVec, Xobs = X80 * S)
    
    pred20[i] <- mean(stats::predict(model, exclude = 's(id)', 
                                     newdata = nd20, newdata.guaranteed = TRUE))
    pred80[i] <- mean(stats::predict(model, exclude = 's(id)', 
                                     newdata = nd80, newdata.guaranteed = TRUE))
    predDiff[i] <- pred20[i] - pred80[i]
    
    pred20_mat <- stats::predict(model, exclude = 's(id)', type = 'lpmatrix',
                                 newdata = nd20, newdata.guaranteed = TRUE)[, 1:Nfixed]
    pred80_mat <- stats::predict(model, exclude = 's(id)', type = 'lpmatrix',
                                 newdata = nd80, newdata.guaranteed = TRUE)[, 1:Nfixed]
    
    A <- as.matrix(colMeans(pred20_mat) - colMeans(pred80_mat))
    se_model[i] <- sqrt(t(A) %*% cov_mat %*% A)
  }
  
  list(pred20 = pred20, pred80 = pred80, predDiff = predDiff, 
       se_model = se_model, time = t_range)
}


#' Print Method for fqgam_boot Objects
#' @param x An object of class \code{fqgam_boot}.
#' @param ... Additional arguments (ignored).
#' @export
print.fqgam_boot <- function(x, ...) {
  cat("Functional Quantile GAM Bootstrap Results\n")
  cat("-----------------------------------------\n")
  cat("Bootstrap type:", x$bootstrap_type, "\n")
  cat("Quantile level (tau):", x$tau, "\n")
  cat("Number of bootstrap samples:", nrow(x$pred20), "\n")
  cat("Number of time points:", ncol(x$pred20), "\n")
  cat("Random seed:", x$seed, "\n")
  cat("Mean running time per iteration:", 
      round(mean(x$running_time), 2), "seconds\n")
  invisible(x)
}


#' Print Method for fqgam_ci Objects
#' @param x An object of class \code{fqgam_ci}.
#' @param ... Additional arguments (ignored).
#' @export
print.fqgam_ci <- function(x, ...) {
  cat("Bootstrap Confidence Intervals\n")
  cat("------------------------------\n")
  cat("Method:", x$method, "\n")
  cat("Confidence level:", (1 - x$alpha) * 100, "%\n")
  cat("Bias adjustment:", ifelse(x$bias_adjusted, "Yes", "No"), "\n")
  cat("Number of time points:", length(x$estimate), "\n")
  cat("Mean bootstrap SE:", round(mean(x$se_boot), 4), "\n")
  invisible(x)
}


#' Print Method for fqgam_bias_adj Objects
#' @param x An object of class \code{fqgam_bias_adj}.
#' @param ... Additional arguments (ignored).
#' @export
print.fqgam_bias_adj <- function(x, ...) {
  cat("Bias Adjustment Results\n")
  cat("-----------------------\n")
  cat("Number of time points:", length(x$adjusted), "\n")
  cat("Mean absolute bias:", round(mean(abs(x$bias)), 4), "\n")
  invisible(x)
}
