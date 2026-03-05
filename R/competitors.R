#' Precision-weighted estimator
#'
#' Combines RCT and observational estimates using inverse-variance weights.
#'
#' @param ed Data frame with columns `rctEst`, `obsEst`, `rctVar`, `obsVar`.
#' @param k Integer. Number of strata (kept for interface consistency).
#' @param trueVar Unused. Kept for interface consistency.
#'
#' @return Numeric vector of precision-weighted estimates, one per stratum.
#' @export
precistion.wtd <- function(ed, k, trueVar = NULL) {
  ed$rctVar / (ed$rctVar + ed$obsVar) * ed$obsEst +
    ed$obsVar / (ed$rctVar + ed$obsVar) * ed$rctEst
}


#' Kappa-1 estimator (Rosenman et al. 2023)
#'
#' Scalar shrinkage toward the observational estimate using a single pooled
#' shrinkage factor.
#'
#' @inheritParams precistion.wtd
#' @return Numeric vector of kappa-1 estimates, one per stratum.
#' @export
kappa.1 <- function(ed, k) {
  lambda.1 <- min(sum(ed$rctVar) / sum((ed$rctEst - ed$obsEst)^2), 1)
  ed$obsEst + (1 - lambda.1) * (ed$rctEst - ed$obsEst)
}


#' Kappa-2 estimator (Rosenman et al. 2023)
#'
#' Diagonal shrinkage toward the observational estimate using stratum-specific
#' shrinkage factors.
#'
#' @inheritParams precistion.wtd
#' @return Numeric vector of kappa-2 estimates, one per stratum.
#' @export
kappa.2 <- function(ed, k, trueVar = NULL) {
  lambda.2 <- pmin(pmax(sum(ed$rctVar)^2 /
                          sum(ed$rctVar^2 * (ed$rctEst - ed$obsEst)^2) * diag(ed$rctVar), 0), 1)
  diag(lambda.2) %*% ed$obsEst + diag(1 - lambda.2) %*% ed$rctEst
}


#' Delta-1 estimator (Green et al. 2005)
#'
#' James-Stein-type scalar shrinkage toward the observational estimate.
#'
#' @inheritParams precistion.wtd
#' @return Numeric vector of delta-1 estimates, one per stratum.
#' @export
delta.1 <- function(ed, k) {
  lambda.1 <- max(1 - (k - 2) / sum(1 / ed$rctVar * (ed$rctEst - ed$obsEst)^2), 0)
  ed$obsEst + lambda.1 * (ed$rctEst - ed$obsEst)
}


#' Delta-2 estimator (Green et al. 2005)
#'
#' James-Stein-type diagonal shrinkage toward the observational estimate.
#'
#' @inheritParams precistion.wtd
#' @return Numeric vector of delta-2 estimates, one per stratum.
#' @export
delta.2 <- function(ed, k) {
  lambda.2 <- pmax(diag(k) - (k - 2) / sum(1 / ed$rctVar^2 * (ed$rctEst - ed$obsEst)^2) *
                     diag(1 / ed$rctVar), 0)
  ed$obsEst + lambda.2 %*% (ed$rctEst - ed$obsEst)
}
