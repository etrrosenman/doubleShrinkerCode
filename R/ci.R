#' Empirical Bayes confidence interval half-widths
#'
#' Constructs confidence interval half-widths for the double shrinker estimator
#' using the critical value adjustment from the `ebci` package.
#'
#' @param ebData List returned by an `eb.*` estimator when called with
#'   `returnCIComps = TRUE`. Must contain elements `c`, `a`, and `lambda`.
#' @param ed Data frame with columns `rctVar` and `obsVar`.
#' @param alpha Numeric. Significance level. Default `0.05`.
#'
#' @return Numeric vector of half-widths for each stratum's confidence interval.
#' @importFrom ebci cva
#' @importFrom stats qnorm
#' @export
ebCiFunc <- function(ebData, ed, alpha = 0.05) {
  sapply(ebData$c, FUN = function(x) {cva(x, alpha = alpha)$cv}) *
    ebData$a * sqrt(ebData$lambda^2 * ed$rctVar +
                      (1 - ebData$lambda)^2 * ed$obsVar)
}


#' Confidence interval coverage check
#'
#' Tests whether true values fall within symmetric confidence intervals centered
#' at point estimates.
#'
#' @param ptEst Numeric vector of point estimates.
#' @param interval Numeric vector of half-widths.
#' @param trueVals Numeric vector of true parameter values.
#'
#' @return Logical vector; `TRUE` where the interval covers the true value.
#' @export
coverage <- function(ptEst, interval, trueVals) {
  ptEst - interval < trueVals &
    trueVals < ptEst + interval
}
