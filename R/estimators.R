#' Core double shrinker estimator
#'
#' Computes the double shrinkage estimator given hyperparameters and per-stratum
#' estimates from both an RCT and an observational study.
#'
#' @param gammaSq Non-negative scalar. Variance of the bias (discrepancy between
#'   RCT and observational treatment effects).
#' @param etaSq Non-negative scalar. Variance of the true treatment effects
#'   (signal variance).
#' @param ed Data frame with columns `rctEst`, `obsEst`, `rctVar`, `obsVar` —
#'   per-stratum point estimates and variances from the RCT and observational study.
#'
#' @return Numeric vector of double-shrunken treatment effect estimates, one per stratum.
#' @export
shrinker <- function(gammaSq, etaSq, ed) {

  x <- ed$rctEst
  y <- ed$obsEst
  var.r <- ed$rctVar
  var.o <- ed$obsVar

  lambda <- (gammaSq + var.o) / (gammaSq + var.o + var.r)
  a <- etaSq * (gammaSq + var.r + var.o) /
    (var.r * (gammaSq + var.o) + etaSq * (gammaSq + var.r + var.o))

  return(a * (lambda * x + (1 - lambda) * y))
}


#' Double shrinker: method of moments estimator 1 (eb.mm1)
#'
#' Estimates hyperparameters `etaSq` and `gammaSq` via method of moments using
#' cross-terms of RCT and observational estimates, then applies the double shrinker.
#'
#' @param ed Data frame with columns `rctEst`, `obsEst`, `rctVar`, `obsVar`.
#' @param k Integer. Number of strata (kept for interface consistency).
#' @param lowerBound Logical. If `TRUE`, applies a lower bound to the hyperparameter
#'   estimates to prevent near-zero shrinkage. Default `FALSE`.
#' @param returnCIComps Logical. If `TRUE`, returns a list with the shrunken estimates
#'   and intermediate quantities needed for confidence interval construction.
#'   Default `FALSE`.
#'
#' @return If `returnCIComps = FALSE`, a numeric vector of shrunken estimates.
#'   If `returnCIComps = TRUE`, a list with elements `shrinker`, `lambda`, `a`, `c`,
#'   `etaSqLwr`, and `gammaSqLwr`.
#' @export
eb.mm1 <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {

  etaSq.pre <- max(mean(ed$rctEst^2 - ed$rctVar), 0)
  gammaSq.pre <- max(mean((ed$rctEst - ed$obsEst)^2 - ed$rctVar - ed$obsVar), 0)

  if (lowerBound) {
    etaSq <- max(etaSq.pre, 2 * sum(ed$rctVar^2) / nrow(ed) / sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2 * sum((ed$obsVar + ed$rctVar)^2) /
                     nrow(ed) / sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }

  if (!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar) / (gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq * (gammaSq + ed$obsVar + ed$rctVar) /
                  (ed$rctVar * (gammaSq + ed$obsVar) + etaSq * (gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar * (gammaSq * (etaSq + 2 * ed$obsVar) + gammaSq^2 + ed$obsVar^2) /
                  (etaSq * ((gammaSq + ed$obsVar)^2 + ed$obsVar * ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}


#' Double shrinker: method of moments estimator 2 (eb.mm2)
#'
#' Alternative method of moments estimator using a different set of moment
#' conditions for `gammaSq` estimation.
#'
#' @inheritParams eb.mm1
#' @inherit eb.mm1 return
#' @export
eb.mm2 <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {

  etaSq.pre <- max(mean(ed$rctEst^2 - ed$rctVar), 0)
  gammaSq.pre <- max(mean(-(ed$rctEst^2) + (ed$obsEst)^2 + ed$rctVar - ed$obsVar), 0)

  if (lowerBound) {
    etaSq <- max(etaSq.pre, 2 * sum(ed$rctVar^2) / nrow(ed) / sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2 * sum((ed$obsVar + ed$rctVar)^2) /
                     nrow(ed) / sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }

  if (!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar) / (gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq * (gammaSq + ed$obsVar + ed$rctVar) /
                  (ed$rctVar * (gammaSq + ed$obsVar) + etaSq * (gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar * (gammaSq * (etaSq + 2 * ed$obsVar) + gammaSq^2 + ed$obsVar^2) /
                  (etaSq * ((gammaSq + ed$obsVar)^2 + ed$obsVar * ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}


#' Double shrinker: maximum likelihood estimator (eb.mle)
#'
#' Estimates hyperparameters by maximizing the empirical Bayes marginal likelihood
#' using `constrOptim` and `L-BFGS-B`, then applies the double shrinker.
#'
#' @inheritParams eb.mm1
#' @inherit eb.mm1 return
#' @importFrom stats constrOptim optim
#' @export
eb.mle <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {

  hyperParams.constrOptim <- tryCatch(
    {constrOptim(c(1e-12, 1e-12), ebLik, grad = ebLikGrad,
                 ui = diag(1, 2), ci = c(0, 0), ed = ed)$par},
    error = function(e) {c(0, 0)})
  hyperParams.constrOptim <- pmax(hyperParams.constrOptim, 0)

  hyperParams.optim <- optim(c(0, 0), ebLik, gr = ebLikGrad, lower = c(0, 0),
                             method = 'L-BFGS-B', ed = ed)$par
  hyperParams.optim <- pmax(hyperParams.optim, 0)

  if (ebLik(hyperParams.constrOptim, ed) < ebLik(hyperParams.optim, ed)) {
    gammaSq.pre <- hyperParams.constrOptim[1]
    etaSq.pre <- hyperParams.constrOptim[2]
  } else {
    gammaSq.pre <- hyperParams.optim[1]
    etaSq.pre <- hyperParams.optim[2]
  }

  if (lowerBound) {
    etaSq <- max(etaSq.pre, 2 * sum(ed$rctVar^2) / nrow(ed) / sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2 * sum((ed$obsVar + ed$rctVar)^2) /
                     nrow(ed) / sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }

  if (!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar) / (gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq * (gammaSq + ed$obsVar + ed$rctVar) /
                  (ed$rctVar * (gammaSq + ed$obsVar) + etaSq * (gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar * (gammaSq * (etaSq + 2 * ed$obsVar) + gammaSq^2 + ed$obsVar^2) /
                  (etaSq * ((gammaSq + ed$obsVar)^2 + ed$obsVar * ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}


#' Double shrinker: unbiased risk estimation (eb.ure)
#'
#' Estimates hyperparameters by minimizing an unbiased risk estimate (URE),
#' then applies the double shrinker.
#'
#' @inheritParams eb.mm1
#' @inherit eb.mm1 return
#' @importFrom stats constrOptim optim
#' @export
eb.ure <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {

  hyperParams.constrOptim <- tryCatch(
    {constrOptim(c(1e-6, 1e-6), URE, grad = UREgrad,
                 ui = diag(1, 2), ci = c(0, 0), ed = ed)$par},
    error = function(e) {c(0, 0)})
  hyperParams.optim <- optim(c(0, 0), URE, lower = c(0, 0), method = 'L-BFGS-B',
                             ed = ed)$par

  if (URE(hyperParams.constrOptim, ed) < URE(hyperParams.optim, ed)) {
    gammaSq.pre <- hyperParams.constrOptim[1]
    etaSq.pre <- hyperParams.constrOptim[2]
  } else {
    gammaSq.pre <- hyperParams.optim[1]
    etaSq.pre <- hyperParams.optim[2]
  }

  if (lowerBound) {
    etaSq <- max(etaSq.pre, 2 * sum(ed$rctVar^2) / nrow(ed) / sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2 * sum((ed$obsVar + ed$rctVar)^2) /
                     nrow(ed) / sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }

  if (!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar) / (gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq * (gammaSq + ed$obsVar + ed$rctVar) /
                  (ed$rctVar * (gammaSq + ed$obsVar) + etaSq * (gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar * (gammaSq * (etaSq + 2 * ed$obsVar) + gammaSq^2 + ed$obsVar^2) /
                  (etaSq * ((gammaSq + ed$obsVar)^2 + ed$obsVar * ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}
