#' @title Bootstrap Methods for Functional Quantile GAM
#'
#' @description
#' This file implements bootstrap procedures for inference in functional quantile 
#' regression models for longitudinal data. Two bootstrap schemes are provided:
#'
#' 1. **Block Bootstrap** (`boot_pred_block`): Resamples complete subject 
#'    trajectories with replacement. This preserves the within-subject dependence 
#'    structure and is used primarily for **variance estimation** and construction 
#'    of standard errors.
#'    
#' 2. **Wild Bootstrap** (`boot_pred_wild`): Combines standard resampling of 
#'    estimated random effects with wild bootstrap of residuals. This method is 
#'    used primarily for **bias adjustment** of the estimators.
#'
#' @section Theoretical Background:
#'
#' The functional quantile regression model is:
#' \deqn{Q_{Y_{ij}|X_{ij}, u_i}^\tau(t_{ij}) = \alpha^\tau(t_{ij}) + \int_S \beta^\tau(s, t_{ij}) X_{ij}(s) ds + u_i}
#'
#' where:
#' \itemize{
#'   \item \eqn{\tau \in (0,1)} is the quantile level
#'   \item \eqn{\alpha^\tau(t)} is the time-varying intercept
#'   \item \eqn{\beta^\tau(s,t)} is the bivariate functional coefficient
#'   \item \eqn{X_{ij}(s)} are functional covariates observed on domain \eqn{S}
#'   \item \eqn{u_i} are subject-specific random intercepts with \eqn{E[u_i] = 0}
#' }
#'
#' @section Why Two Bootstrap Methods:
#'
#' As documented in Battagliola et al. (2022), estimators in quantile regression 
#' for longitudinal data can be biased when clusters are small. The bias arises 
#' from a combination of:
#' \itemize{
#'   \item The incidental parameter problem (number of parameters grows with N)
#'   \item Nonlinearity of quantiles
#'   \item Penalization of subject-specific intercepts
#' }
#'
#' Block resampling cannot be used for bias adjustment because the target parameter 
#' is not computable under the bootstrap distribution. Instead, the wild bootstrap 
#' generates data under a distribution where the true parameter equals the estimate 
#' from the observed data, allowing bias estimation.
#'
#' The recommended approach (Battagliola et al. 2024) combines:
#' \itemize{
#'   \item Wild bootstrap for bias adjustment
#'   \item Block bootstrap for variance/standard error estimation
#' }
#'
#' @section References:
#'
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2024).
#' Quantile regression for longitudinal functional data with application to
#' feed intake of lactating sows. Journal of Agricultural, Biological, and 
#' Environmental Statistics, 30(1), 211-230.
#'
#'
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2022).
#' A bias-adjusted estimator in quantile regression for clustered data.
#' Econometrics and Statistics, 23, 165-186.
#' 
#'
#' @name bootstrap-methods
NULL


#' Create Block Bootstrapped Dataset
#'
#' Internal function to create a block bootstrapped dataset by resampling
#' complete subjects with replacement. This preserves the within-subject 
#' dependence structure inherent in longitudinal data.
#'
#' @details
#' The block bootstrap procedure:
#' \enumerate{
#'   \item Sample N subject IDs with replacement from \code{unique(mydata$id)}
#'   \item For each sampled ID, extract all observations for that subject
#'   \item If a subject is sampled multiple times, create copies with modified IDs
#'   \item Stack all extracted observations into the bootstrap dataset
#' }
#'
#' This approach maintains the correlation structure within subjects because 
#' complete subject trajectories are resampled together. It is equivalent to 
#' the "cross-sectional resampling" described in Galvao and Montes-Rojas (2015).
#'
#' @param mydata A data.frame containing the original data with columns:
#'   \describe{
#'     \item{id}{Subject identifier (factor or character)}
#'     \item{y}{Response variable}
#'     \item{time}{Longitudinal time stamps}
#'     \item{Xobs}{Matrix of functional covariate observations}
#'     \item{sMat}{Matrix of functional domain coordinates}
#'     \item{tMat}{Matrix of longitudinal time coordinates}
#'   }
#'
#' @return A data.frame with the same structure as the input, containing
#'   resampled data. Subject IDs are modified to be unique (e.g., "5.1", "5.2" 
#'   if subject 5 is sampled twice).
#'
#' @references
#' Galvao, A., & Montes-Rojas, G. (2015). On bootstrap inference for quantile 
#' regression panel data: A Monte Carlo study. Econometrics, 3(3), 654-666.
#'
#' @keywords internal
new_dataset_boot <- function(mydata) {
  
  S <- dim(mydata$Xobs)[2]
  N <- length(unique(mydata$id))
  
  new_ind <- sample(unique(mydata$id), replace = TRUE)
  freq <- data.frame(table(new_ind))
  
  new_data <- NULL
  new_Xobs <- NULL
  new_sMat <- NULL
  new_tMat <- NULL
  
  for (i in 1:N) {
    if (freq$Freq[i] != 0) {
      ind <- which(mydata$id == freq$new_ind[i])
      for (j in 1:freq$Freq[i]) {
        new_data <- rbind.data.frame(
          new_data, 
          cbind.data.frame(
            id = rep(paste(freq$new_ind[i], ".", j, sep = ""), length(ind)),
            y = mydata$y[ind],
            time = mydata$time[ind]
          )
        )
        new_Xobs <- rbind(new_Xobs, mydata$Xobs[ind, ])
        new_sMat <- rbind(new_sMat, mydata$sMat[ind, ])
        new_tMat <- rbind(new_tMat, mydata$tMat[ind, ])
      }
    }
  }
  
  new_data$id <- factor(new_data$id)
  new_data$Xobs <- new_Xobs
  new_data$sMat <- new_sMat
  new_data$tMat <- new_tMat
  
  return(new_data)
}


#' Block Bootstrap for Functional Quantile GAM
#'
#' @description
#' Performs block bootstrap inference for functional quantile regression models
#' fitted with \code{\link[qgam]{qgam}}. The bootstrap resamples complete subject 
#' trajectories with replacement, preserving the longitudinal dependence structure.
#'
#' @details
#' \strong{Purpose:}
#' Block bootstrap is primarily used for variance estimation and construction 
#' of standard errors. As noted in Battagliola et al. (2024), model-based standard 
#' errors from \code{qgam} may underestimate the true sampling variation due to 
#' penalization of random effects.
#'
#' \strong{Algorithm:}
#' For each bootstrap iteration b = 1, ..., B:
#' (1) Create a bootstrap dataset by sampling N subjects with replacement.
#' (2) Fit the quantile regression model to the bootstrap data.
#' (3) Compute predictions at reference curves X20 and X80 for each time point.
#' (4) Store the difference in predictions.
#'
#' The bootstrap standard error is: sd_boot = sd(theta_1, ..., theta_B).
#'
#' \strong{Predictions:}
#' For a reference functional covariate X(s) at longitudinal time t, 
#' the predicted quantile (excluding random effects) represents the 
#' tau-quantile for a typical subject with u_i = 0.
#'
#' \strong{Important Note on Bias:}
#' Block bootstrap estimates are centered around the original estimate, so this 
#' method is not suitable for bias adjustment. For bias correction, use 
#' \code{\link{boot_pred_wild}} instead.
#'
#' @param d A data.frame used to fit the model. Must contain columns:
#'   \describe{
#'     \item{id}{Subject identifier (factor). Observations from the same subject 
#'       must share the same id.}
#'     \item{y}{Response variable (numeric)}
#'     \item{time}{Longitudinal time stamps (integer or numeric)}
#'     \item{Xobs}{Matrix of functional covariate observations (n x S matrix, 
#'       where S is the number of functional grid points)}
#'     \item{sMat}{Matrix of functional domain coordinates (n x S matrix)}
#'     \item{tMat}{Matrix of longitudinal time coordinates (n x S matrix)}
#'   }
#' @param B Integer. Number of bootstrap iterations. Battagliola et al. (2024) 
#'   use B = 100 in their application.
#' @param seed Integer. Random seed for reproducibility.
#' @param tau Numeric. Quantile level in (0,1). Should match the quantile 
#'   level used in the fitted model.
#' @param X20 Numeric vector of length S. Point-wise lower reference curve 
#'   (e.g., 0.2 quantile of functional covariates).
#' @param X80 Numeric vector of length S. Point-wise upper reference curve 
#'   (e.g., 0.8 quantile of functional covariates).
#' @param model A fitted \code{qgam} model object.
#' @param multicore Logical. If \code{TRUE}, uses parallel computation in 
#'   \code{qgam}. Default is \code{TRUE}.
#' @param ncores Integer. Number of cores for parallel computation. 
#'   Default is 3.
#' @param verbose Logical. If \code{TRUE}, prints progress messages. 
#'   Default is \code{TRUE}.
#'
#' @return A list of class \code{"fqgam_boot"} with components:
#'   \describe{
#'     \item{pred20}{Matrix (B x maxT) of predicted quantiles at each 
#'       longitudinal time for X20}
#'     \item{pred80}{Matrix (B x maxT) of predicted quantiles at each 
#'       longitudinal time for X80}
#'     \item{predDiff}{Matrix (B x maxT) of differences \eqn{\hat{Q}_{20}^\tau(t) - \hat{Q}_{80}^\tau(t)}}
#'     \item{seed}{The random seed used}
#'     \item{d}{The input data.frame}
#'     \item{tau}{The quantile level}
#'     \item{running_time}{Vector of execution times (in seconds) for each iteration}
#'     \item{bootstrap_type}{"block"}
#'   }
#'
#' @seealso 
#' \code{\link{boot_pred_wild}} for wild bootstrap (bias adjustment),
#' \code{\link{boot_ci}} for computing confidence intervals,
#' \code{\link[qgam]{qgam}} for model fitting
#'
#' @references
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2024).
#' Quantile regression for longitudinal functional data with application to
#' feed intake of lactating sows. Journal of Agricultural, Biological, and 
#' Environmental Statistics, 30(1), 211-230.
#'
#' @examples
#' \dontrun{
#' # Prepare data and reference curves
#' # Run block bootstrap for variance estimation
#' bb <- boot_pred_block(myDataFit, B = 100, seed = 1234, tau = 0.1,
#'                       X20 = cca20, X80 = cca80, model = fit)
#' }
#'
#' @export
boot_pred_block <- function(d, B, seed, tau, X20, X80, model,
                            multicore = TRUE, ncores = 3, verbose = TRUE) {
  
  S <- dim(d$Xobs)[2]
  N <- length(unique(d$id))
  sVec <- seq(0, 1, length = S)
  maxT <- max(d$time)
  myformula <- stats::formula(model)
  
  set.seed(seed)
  mySeed <- sample(1:10^7, size = B, replace = FALSE)
  
  pred20 <- matrix(NA, B, maxT)
  pred80 <- matrix(NA, B, maxT)
  predDiff <- matrix(NA, B, maxT)
  running_time <- rep(NA, B)
  
  for (b in 1:B) {
    start_time <- Sys.time()
    set.seed(mySeed[b])
    
    if (verbose) {
      message(paste("Block bootstrap iteration", b, "of", B))
    }
    
    new_d <- new_dataset_boot(d)
    fit <- qgam::qgam(myformula, qu = tau, data = new_d, 
                      multicore = multicore, ncores = ncores)
    
    for (t in 1:maxT) {
      nd20 <- data.frame(time = rep(t, S), tMat = rep(t, S), sMat = sVec, Xobs = X20 * S)
      nd80 <- data.frame(time = rep(t, S), tMat = rep(t, S), sMat = sVec, Xobs = X80 * S)
      
      pred20[b, t] <- mean(stats::predict(fit, exclude = 's(id)', newdata = nd20, newdata.guaranteed = TRUE))
      pred80[b, t] <- mean(stats::predict(fit, exclude = 's(id)', newdata = nd80, newdata.guaranteed = TRUE))
      predDiff[b, t] <- pred20[b, t] - pred80[b, t]
    }
    
    end_time <- Sys.time()
    running_time[b] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  res <- list(pred20 = pred20, pred80 = pred80, predDiff = predDiff, seed = seed, 
              d = d, tau = tau, running_time = running_time, bootstrap_type = "block")
  class(res) <- c("fqgam_boot", "list")
  return(res)
}


#' Wild Bootstrap for Functional Quantile GAM
#'
#' @description
#' Performs wild bootstrap inference for functional quantile regression models
#' fitted with \code{\link[qgam]{qgam}}. This method combines standard resampling 
#' of estimated random effects with wild bootstrap of residuals, and is primarily 
#' used for bias adjustment.
#'
#' @details
#' \strong{Purpose:}
#' Wild bootstrap is used for bias adjustment of the quantile regression 
#' estimators. As documented in Battagliola et al. (2022), estimators can be 
#' severely biased when clusters are small.
#'
#' \strong{Algorithm:}
#' For each bootstrap iteration b = 1, ..., B:
#' (1) Extract estimated terms and residuals from the fitted model.
#' (2) Resample random effects with replacement from the estimated values.
#' (3) Generate wild bootstrap weights w independently from a two-point 
#'     distribution: w = 2(1-tau) with probability (1-tau), and w = -2*tau 
#'     with probability tau. This distribution has tau-quantile equal to 0.
#' (4) Create bootstrap responses: Y* = fitted_fixed + u* + w*|residual|.
#' (5) Fit the model to (Y*, X, t) and compute predictions.
#'
#' @inheritParams boot_pred_block
#'
#' @return A list of class \code{"fqgam_boot"} with components:
#'   \describe{
#'     \item{pred20}{Matrix (B x maxT) of predicted quantiles at each 
#'       longitudinal time for X20}
#'     \item{pred80}{Matrix (B x maxT) of predicted quantiles at each 
#'       longitudinal time for X80}
#'     \item{predDiff}{Matrix (B x maxT) of differences \eqn{\hat{Q}_{20}^\tau(t) - \hat{Q}_{80}^\tau(t)}}
#'     \item{seed}{The random seed used}
#'     \item{d}{The input data.frame}
#'     \item{tau}{The quantile level}
#'     \item{running_time}{Vector of execution times (in seconds) for each iteration}
#'     \item{bootstrap_type}{"wild"}
#'   }
#'
#' @seealso 
#' \code{\link{boot_pred_block}} for block bootstrap (variance estimation),
#' \code{\link{boot_ci}} for computing confidence intervals
#'
#' @references
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2022).
#' A bias-adjusted estimator in quantile regression for clustered data.
#' Econometrics and Statistics, 23, 165-186.
#'
#' @export
boot_pred_wild <- function(d, B, seed, tau, X20, X80, model,
                           multicore = TRUE, ncores = 3, verbose = TRUE) {
  
  S <- dim(d$Xobs)[2]
  N <- length(unique(d$id))
  M <- length(d$id)
  sVec <- seq(0, 1, length = S)
  maxT <- max(d$time)
  
  terms <- stats::predict(model, type = 'terms')
  eps <- d$y - stats::fitted(model)
  ind <- which(colnames(terms) == "s(id)")
  u <- unique(terms[, ind])
  myformula <- stats::formula(model)
  
  set.seed(seed)
  mySeed <- sample(1:10^7, size = B, replace = FALSE)
  pred20 <- matrix(NA, B, maxT)
  pred80 <- matrix(NA, B, maxT)
  predDiff <- matrix(NA, B, maxT)
  running_time <- rep(NA, B)
  orig.id <- unique(d$id)
  
  for (b in 1:B) {
    start_time <- Sys.time()
    set.seed(mySeed[b])
    
    if (verbose) {
      message(paste("Wild bootstrap iteration", b, "of", B))
    }
    
    bootData <- d
    w <- sample(c(2 * (1 - tau), -2 * tau), size = M, replace = TRUE, prob = c(1 - tau, tau))
    u.hat <- sample(u, replace = TRUE)
    uBoot <- NULL
    for (i in 1:N) {
      len <- length(which(d$id == orig.id[i]))
      uBoot <- c(uBoot, rep(u.hat[i], len))
    }
    
    if (!is.null(dim(terms[, -ind]))) {
      bootData$y <- model$coefficients[1] + rowSums(terms[, -ind]) + uBoot + w * abs(eps)
    } else {
      bootData$y <- model$coefficients[1] + terms[, -ind] + uBoot + w * abs(eps)
    }
    
    fit <- qgam::qgam(myformula, qu = tau, data = bootData, 
                      multicore = multicore, ncores = ncores)
    
    for (t in 1:maxT) {
      nd20 <- data.frame(time = rep(t, S), tMat = rep(t, S), sMat = sVec, Xobs = X20 * S)
      nd80 <- data.frame(time = rep(t, S), tMat = rep(t, S), sMat = sVec, Xobs = X80 * S)
      pred20[b, t] <- mean(stats::predict(fit, exclude = 's(id)', newdata = nd20, newdata.guaranteed = TRUE))
      pred80[b, t] <- mean(stats::predict(fit, exclude = 's(id)', newdata = nd80, newdata.guaranteed = TRUE))
      predDiff[b, t] <- pred20[b, t] - pred80[b, t]
    }
    
    end_time <- Sys.time()
    running_time[b] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  res <- list(pred20 = pred20, pred80 = pred80, predDiff = predDiff, seed = seed, 
              d = d, tau = tau, running_time = running_time, bootstrap_type = "wild")
  class(res) <- c("fqgam_boot", "list")
  return(res)
}