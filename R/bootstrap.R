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
#' \deqn{Q_{Y_{ij}|X_{ij}, u_i}^\tau(t_{ij}) = \alpha^\tau(t_{ij}) + 
#'       \int_S \beta^\tau(s, t_{ij}) X_{ij}(s) ds + u_i}
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
#' \doi{10.1007/s13253-024-00601-5}
#' 
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2022).
#' A bias-adjusted estimator in quantile regression for clustered data.
#' Econometrics and Statistics, 23, 165-186.
#' \doi{10.1016/j.ecosta.2021.07.003}
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
  
  ## Length of functional grid
  S <- dim(mydata$Xobs)[2]
  ## Number of subjects
  N <- length(unique(mydata$id))
  
  ## Sample the subjects id with replacement 
  new_ind <- sample(unique(mydata$id), replace = TRUE)
  
  freq <- data.frame(table(new_ind))
  
  new_data <- NULL
  new_Xobs <- NULL
  new_sMat <- NULL
  new_tMat <- NULL
  
  ## For every subject id, if it is included in new_ind, 
  ## extract the data relative to it and stack it in new_data
  for (i in 1:N) {
    
    if (freq$Freq[i] != 0) {
      
      ind <- which(mydata$id == freq$new_ind[i])
      
      for (j in 1:freq$Freq[i]) {
        # Create unique ID for repeated subjects (e.g., "5.1", "5.2")
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
#'   (e.g., 20\% quantile of functional covariates).
#' @param X80 Numeric vector of length S. Point-wise upper reference curve 
#'   (e.g., 80\% quantile of functional covariates).
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
#' Galvao, A., & Montes-Rojas, G. (2015). On bootstrap inference for quantile 
#' regression panel data: A Monte Carlo study. Econometrics, 3(3), 654-666.
#'
#' @examples
#' \dontrun{
#' library(refund)
#' library(qgam)
#' 
#' # Prepare DTI data
#' data(DTI)
#' myData <- DTI[complete.cases(DTI$cca), ]
#' myData <- subset(myData, Nscans > 1)
#' 
#' # Create data.frame for qgam
#' N <- nrow(myData)
#' S <- ncol(myData$cca)
#' sVec <- seq(0, 1, length = S)
#' 
#' myDataFit <- data.frame(
#'   y = myData$pasat,
#'   id = factor(myData$ID),
#'   time = myData$visit.time
#' )
#' myDataFit$Xobs <- myData$cca
#' myDataFit$sMat <- matrix(sVec, N, S, byrow = TRUE)
#' myDataFit$tMat <- matrix(myData$visit.time, N, S, byrow = FALSE)
#' 
#' # Fit model at tau = 0.1
#' tau <- 0.1
#' formula <- y ~ s(time, bs = "cr", k = 10) + 
#'               te(sMat, tMat, by = Xobs, bs = c("cr", "cr"), k = c(10, 10)) + 
#'               s(id, bs = "re")
#' fit <- qgam(formula, qu = tau, data = myDataFit)
#' 
#' # Get reference curves (20% and 80% pointwise quantiles)
#' cca20 <- apply(myData$cca, 2, quantile, probs = 0.2)
#' cca80 <- apply(myData$cca, 2, quantile, probs = 0.8)
#' 
#' # Run block bootstrap for variance estimation
#' bb <- boot_pred_block(myDataFit, B = 100, seed = 1234, tau = tau,
#'                       X20 = cca20, X80 = cca80, model = fit)
#' 
#' # Compute bootstrap standard errors
#' sd_block <- apply(bb$predDiff, 2, sd)
#' }
#'
#' @export
boot_pred_block <- function(d, B, seed, tau, X20, X80, model,
                            multicore = TRUE, ncores = 3, verbose = TRUE) {
  
  ## Length of functional grid
  S <- dim(d$Xobs)[2]
  ## Number of subjects
  N <- length(unique(d$id))
  
  ## Functional grid, taken as [0,1]
  sVec <- seq(0, 1, length = S)
  
  ## Maximum time stamp
  maxT <- max(d$time)
  
  ## Formula of the model extracted
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
    
    ## Create a block bootstrapped dataset
    new_d <- new_dataset_boot(d)
    
    ## Fit the model to bootstrap data
    fit <- qgam::qgam(myformula, qu = tau, data = new_d, 
                      multicore = multicore, ncores = ncores)
    
    ## For every longitudinal time stamp, predict the quantile 
    ## for both X20 and X80, as well as their difference
    for (t in 1:maxT) {
      # Create prediction data frames
      # Note: Xobs is multiplied by S for proper integration approximation
      nd20 <- data.frame(
        time = rep(t, S), 
        tMat = rep(t, S), 
        sMat = sVec, 
        Xobs = X20 * S
      )
      nd80 <- data.frame(
        time = rep(t, S), 
        tMat = rep(t, S), 
        sMat = sVec, 
        Xobs = X80 * S
      )
      
      # Predict excluding random effects (for "typical" subject)
      pred20[b, t] <- mean(stats::predict(
        fit, exclude = 's(id)', newdata = nd20, newdata.guaranteed = TRUE
      ))
      
      pred80[b, t] <- mean(stats::predict(
        fit, exclude = 's(id)', newdata = nd80, newdata.guaranteed = TRUE
      ))
      
      predDiff[b, t] <- pred20[b, t] - pred80[b, t]
    }
    
    end_time <- Sys.time()
    running_time[b] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  res <- list(
    pred20 = pred20,
    pred80 = pred80,
    predDiff = predDiff,
    seed = seed,
    d = d,
    tau = tau,
    running_time = running_time,
    bootstrap_type = "block"
  )
  
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
#' severely biased when clusters are small due to the incidental parameter 
#' problem, nonlinearity of quantiles, and penalization.
#' 
#' \strong{Why Wild Bootstrap for Bias:}
#' Block resampling cannot be used for bias adjustment because the target 
#' parameter is not computable under the bootstrap distribution (the resampled 
#' data do not have the same underlying parameters). Wild bootstrap generates 
#' bootstrap data under a distribution where the true parameter equals the 
#' estimate from observed data, allowing bias to be measured as the deviation 
#' between bootstrap estimates and the original estimate.
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
#' \strong{Bias Estimation:}
#' The bias is estimated as: bias = mean(theta_boot - theta_original).
#' The bias-adjusted estimator is: theta_adj = theta - bias = 2*theta - mean(theta_boot).
#' 
#' \strong{Wild Bootstrap Properties:}
#' The wild bootstrap with asymmetric weights was proposed by Feng et al. (2011) 
#' for quantile regression. As noted by Wang et al. (2018), wild bootstrap captures 
#' asymmetry and heteroscedasticity better than ordinary resampling of residuals.
#' 
#' \strong{Important Notes:}
#' Covariates X and time points t are kept fixed (same as in the original data).
#' This method may underestimate variance due to shrinkage of predicted 
#' random effects; use \code{\link{boot_pred_block}} for variance estimation.
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
#' \code{\link{boot_ci}} for computing confidence intervals,
#' \code{\link{compute_bias_adjustment}} for computing bias-adjusted estimates
#'
#' @references
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2024).
#' Quantile regression for longitudinal functional data with application to
#' feed intake of lactating sows. Journal of Agricultural, Biological, and 
#' Environmental Statistics, 30(1), 211-230.
#' 
#' Battagliola, M.L., Sorensen, H., Tolver, A., & Staicu, A.-M. (2022).
#' A bias-adjusted estimator in quantile regression for clustered data.
#' Econometrics and Statistics, 23, 165-186.
#' 
#' Feng, X., He, X., & Hu, J. (2011). Wild bootstrap for quantile regression.
#' Biometrika, 98(4), 995-999.
#' 
#' Wang, L., Van Keilegom, I., & Maidman, A. (2018). Wild residual bootstrap 
#' inference for penalized quantile regression with heteroscedastic errors.
#' Biometrika, 105(4), 859-872.
#'
#' @examples
#' \dontrun{
#' library(refund)
#' library(qgam)
#' 
#' # (Data preparation as in boot_pred_block example)
#' 
#' # Run wild bootstrap for bias adjustment
#' wb <- boot_pred_wild(myDataFit, B = 100, seed = 4321, tau = tau,
#'                      X20 = cca20, X80 = cca80, model = fit)
#' 
#' # Compute bias-adjusted estimates
#' # Original prediction difference
#' predDiff_original <- ... # computed from original fit
#' 
#' # Bias estimation
#' diff_adj <- sweep(wb$predDiff, 2, predDiff_original)
#' bias <- colMeans(diff_adj)
#' 
#' # Bias-adjusted prediction
#' predDiff_adjusted <- predDiff_original - bias
#' }
#'
#' @export
boot_pred_wild <- function(d, B, seed, tau, X20, X80, model,
                           multicore = TRUE, ncores = 3, verbose = TRUE) {
  
  ## Length of functional grid
  S <- dim(d$Xobs)[2]
  ## Number of subjects
  N <- length(unique(d$id))
  ## Total number of observations in the dataset
  M <- length(d$id)
  
  ## Functional grid, taken as [0,1]
  sVec <- seq(0, 1, length = S)
  
  ## Maximum time stamp
  maxT <- max(d$time)
  
  ## ============================================================
  ## STEP 1: Extract components from the fitted model
  ## ============================================================
  
  ## Single estimated terms from the model (intercept, smooths, etc.)
  ## Each column corresponds to one term in the model
  terms <- stats::predict(model, type = 'terms')
  
  ## Residuals: Y - fitted values (not the working residuals from gamObject)
  ## These are: eps_ij = Y_ij - alpha(t) - integral(beta*X) - u_i
  eps <- d$y - stats::fitted(model)
  
  ## Estimated subject-specific intercepts
  ## Find the column corresponding to s(id) random effect
  ind <- which(colnames(terms) == "s(id)")
  
  ## Extract unique random effect values (one per subject)
  ## Note: terms[,ind] preserves the original ordering of d$id
  u <- unique(terms[, ind])
  
  ## Formula of the model extracted
  myformula <- stats::formula(model)
  
  ## ============================================================
  ## STEP 2: Set up bootstrap iterations
  ## ============================================================
  
  set.seed(seed)
  mySeed <- sample(1:10^7, size = B, replace = FALSE)
  
  pred20 <- matrix(NA, B, maxT)
  pred80 <- matrix(NA, B, maxT)
  predDiff <- matrix(NA, B, maxT)
  running_time <- rep(NA, B)
  
  orig.id <- unique(d$id)
  
  ## ============================================================
  ## STEP 3: Bootstrap iterations
  ## ============================================================
  
  for (b in 1:B) {
    start_time <- Sys.time()
    
    set.seed(mySeed[b])
    
    if (verbose) {
      message(paste("Wild bootstrap iteration", b, "of", B))
    }
    
    ## Initialize bootstrapped data (keep X and t fixed)
    bootData <- d
    
    ## --------------------------------------------------------
    ## Generate wild bootstrap weights
    ## Distribution: w = 2(1-tau) with prob (1-tau), -2*tau with prob tau
    ## This ensures the tau-quantile of w is 0
    ## --------------------------------------------------------
    w <- sample(c(2 * (1 - tau), -2 * tau), size = M, replace = TRUE, 
                prob = c(1 - tau, tau))
    
    ## --------------------------------------------------------
    ## Resample the estimated subject-specific intercepts
    ## Sample N values with replacement from the N estimated u_i
    ## --------------------------------------------------------
    u.hat <- sample(u, replace = TRUE)
    
    ## Map resampled u.hat to all observations
    ## Each u.hat[i] should be repeated for all observations of subject i
    uBoot <- NULL
    
    for (i in 1:N) {
      len <- length(which(d$id == orig.id[i]))
      uBoot <- c(uBoot, rep(u.hat[i], len))
    }
    
    ## --------------------------------------------------------
    ## Construct bootstrap response Y*
    ## Y* = intercept + fixed_effects + resampled_random_effect + w*|residual|
    ## --------------------------------------------------------
    
    # Check if there are multiple fixed effect terms (matrix) or just one (vector)
    if (!is.null(dim(terms[, -ind]))) {
      # Multiple terms: sum them rowwise
      bootData$y <- model$coefficients[1] + rowSums(terms[, -ind]) + 
                    uBoot + w * abs(eps)
    } else {
      # Single term (vector)
      bootData$y <- model$coefficients[1] + terms[, -ind] + 
                    uBoot + w * abs(eps)
    }
    
    ## --------------------------------------------------------
    ## Fit the model with bootstrapped data
    ## --------------------------------------------------------
    fit <- qgam::qgam(myformula, qu = tau, data = bootData, 
                      multicore = multicore, ncores = ncores)
    
    ## --------------------------------------------------------
    ## Compute predictions for each longitudinal time stamp
    ## --------------------------------------------------------
    for (t in 1:maxT) {
      nd20 <- data.frame(
        time = rep(t, S), 
        tMat = rep(t, S), 
        sMat = sVec, 
        Xobs = X20 * S
      )
      nd80 <- data.frame(
        time = rep(t, S), 
        tMat = rep(t, S), 
        sMat = sVec, 
        Xobs = X80 * S
      )
      
      pred20[b, t] <- mean(stats::predict(
        fit, exclude = 's(id)', newdata = nd20, newdata.guaranteed = TRUE
      ))
      
      pred80[b, t] <- mean(stats::predict(
        fit, exclude = 's(id)', newdata = nd80, newdata.guaranteed = TRUE
      ))
      
      predDiff[b, t] <- pred20[b, t] - pred80[b, t]
    }
    
    end_time <- Sys.time()
    running_time[b] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  res <- list(
    pred20 = pred20,
    pred80 = pred80,
    predDiff = predDiff,
    seed = seed,
    d = d,
    tau = tau,
    running_time = running_time,
    bootstrap_type = "wild"
  )
  
  class(res) <- c("fqgam_boot", "list")
  
  return(res)
}
