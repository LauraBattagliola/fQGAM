#' fQGAM: Functional Quantile Regression via Generalized Additive Models
#'
#' The fQGAM package provides tools for fitting functional quantile regression
#' models for longitudinal functional data using generalized additive models.
#'
#' @section Model:
#' The functional quantile regression model is specified as:
#' \deqn{Q_{Y_{ij}| X_{ij}, u_i}^\tau(t_{ij}) = \alpha^\tau(t_{ij}) + 
#'       \int_S \beta^\tau(s,t_{ij}) X_{ij}(s)ds + u_i}
#' where \eqn{\tau \in (0,1)} is the quantile level, \eqn{X_{ij}(\cdot)} are 
#' functional covariates, and \eqn{u_i} are subject-specific random effects.
#' 
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{prep_fqgam_data}}}{Prepare data for model fitting}
#'   \item{\code{\link{boot_pred_block}}}{Block bootstrap for variance estimation}
#'   \item{\code{\link{boot_pred_wild}}}{Wild bootstrap for bias adjustment}
#'   \item{\code{\link{compute_bias_adjustment}}}{Compute bias-adjusted estimates}
#'   \item{\code{\link{boot_ci}}}{Compute bootstrap confidence intervals}
#'   \item{\code{\link{predict_fqgam}}}{Predict quantiles at reference curves}
#'   \item{\code{\link{compute_model_se}}}{Compute model-based standard errors}
#' }
#' 
#' @section Bootstrap Methods:
#' Two bootstrap schemes are implemented:
#' \itemize{
#'   \item **Block Bootstrap**: Resamples complete subject trajectories with 
#'     replacement. Used for variance estimation.
#'   \item **Wild Bootstrap**: Combines resampling of random effects with wild 
#'     bootstrap of residuals. Used for bias adjustment.
#' }
#' 
#' The recommended approach combines both methods:
#' \itemize{
#'   \item Use wild bootstrap to estimate and correct for bias
#'   \item Use block bootstrap to estimate standard errors
#'   \item Construct confidence intervals as: estimate - bias +/- z * SE
#' }
#'
#' @section Dependencies:
#' This package relies on the \code{qgam} package for quantile regression
#' and \code{mgcv} for the underlying GAM machinery.
#'
#' @references
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
#' @docType package
#' @name fQGAM-package
#' @aliases fQGAM
#' @keywords package
"_PACKAGE"
