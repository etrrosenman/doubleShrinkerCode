


#' Estimate per-stratum causal effects from RCT and observational data
#'
#' Computes stratum-level treatment effect estimates and variances for both
#' an RCT sample and an observational study sample. The observational estimate
#' uses inverse probability weighting when `propScoreAdjust = TRUE`.
#'
#' @param rctSample Data frame of RCT data. Must have columns `subgroup`, `Test`,
#'   and the outcome column named by `outcome`.
#' @param obsSample Data frame of observational data. Must have columns `subgroup`,
#'   `Test`, `propScore`, and the outcome column named by `outcome`.
#' @param outcome Character string. Name of the binary outcome column.
#' @param propScoreAdjust Logical. If `TRUE` (default), uses inverse probability
#'   weighting for the observational estimate.
#'
#' @return Data frame with one row per stratum and columns `stratum`, `rctEst`,
#'   `rctVar`, `obsEst`, `obsVar`.
#' @importFrom plyr ddply .
#' @export
getParamEstimates <- function(rctSample, obsSample, outcome, propScoreAdjust = TRUE) {
  
  # build the data frame stratum by stratum 
  strata <- sort(unique(rctSample$subgroup))
  
  # rct 
  rctStrataResults <- ddply(rctSample, .(subgroup), .fun = function(d) {
    
    # rct
    n1 <- sum(d$Test)
    n0 <- sum(1 - d$Test)
    p1 <- ifelse((mu1 <- mean(d[[outcome]][d$Test == 1])) == 0, 
                 0.5/n1, ifelse(mu1 == 1, (n1 - 0.5)/n1, mu1))
    p0 <- ifelse((mu0 <- mean(d[[outcome]][d$Test == 0])) == 0, 
                 0.5/n0, ifelse(mu0 == 1, (n0 - 0.5)/n0, mu0))
    rctEst <- ifelse(mu1 == 0 && mu0 == 0, 0, p1 - p0)
    rctVar <- p1*(1-p1)/ifelse(n1 > 1, n1 - 1, n1) + p0*(1-p0)/ifelse(n0 > 1, n0 - 1, n0)
    
    c(rctEst = rctEst, rctVar = rctVar)
  })
  
  # obs
  obsStrataResults <- ddply(obsSample, .(subgroup), .fun = function(d) {
    
    # obs
    if(propScoreAdjust) {
      p1 <- sum(d[[outcome]][d$Test == 1]/d$propScore[d$Test == 1])/
        sum(1/d$propScore[d$Test == 1])
      p0 <- sum(d[[outcome]][d$Test == 0]/(1-d$propScore[d$Test == 0]))/
        sum(1/(1 - d$propScore[d$Test == 0]))
      
      obsEst <- p1 - p0 
      
      wts.trt <- 1/d$propScore[d$Test == 1]/sum(1/d$propScore[d$Test == 1])
      wts.ctl <-  1/(1-d$propScore[d$Test == 0])/sum(1/(1-d$propScore[d$Test == 0]))
      
      obsVar <- (sum(wts.trt^2*p1*(1-p1)) + sum(wts.ctl^2*p0*(1-p0)))
      
    } else {
      p1 <- mean(d[[outcome]][d$Test == 1]) 
      p0 <- mean(d[[outcome]][d$Test == 0]) 
      
      obsEst <- p1 - p0 
      obsVar <- p1*(1-p1)/sum(d$Test == 1) + p0*(1-p0)/sum(d$Test == 0)
    }
    
    c(obsEst = obsEst, obsVar = obsVar)
    
    
  })
  
  return(data.frame(stratum = sort(unique(rctSample$subgroup)), 
                    rctEst = rctStrataResults$rctEst, rctVar = rctStrataResults$rctVar, 
                    obsEst = obsStrataResults$obsEst, obsVar = obsStrataResults$obsVar))
}


#' Define subgroups and compute true treatment effects
#'
#' Tags units in both datasets with a subgroup label based on one or more
#' variables, and computes the "true" per-stratum treatment effects from the
#' RCT population.
#'
#' @param subgroupVars Character vector of column names to use for subgrouping.
#'   If multiple names are provided, subgroups are formed by crossing all variables.
#' @param data.os Data frame of observational study data.
#' @param data.rct Data frame of RCT data.
#' @param outcome Character string. Name of the binary outcome column used to
#'   compute true treatment effects.
#'
#' @return A list with elements:
#'   - `data.os`: observational data with an added `subgroup` column.
#'   - `data.rct`: RCT data with an added `subgroup` column.
#'   - `subgroupTaus.true`: named numeric vector of true per-stratum treatment effects.
#' @export
subgroupDefs_noSplit <- function(subgroupVars, data.os, data.rct, outcome ) {
  
  # define the subgroups in the observational study and RCT
  if(length(subgroupVars) > 1) {
    data.os$subgroup <- factor(apply(data.os[,subgroupVars], 1, FUN = function(x) {
      paste(x, collapse = " ")
    }))
  } else {
    data.os$subgroup <- factor(data.os[,subgroupVars])
  }
  
  if(length(subgroupVars) > 1) {
    data.rct$subgroup <- factor(apply(data.rct[,subgroupVars], 1, FUN = function(x) {
      paste(x, collapse = " ")
    }))
  } else {
    data.rct$subgroup <- factor(data.rct[,subgroupVars])
  }
  
  # estimate "true" causal effects
  subgroupTaus.gold <- sapply(sort(unique(data.rct$subgroup)), FUN = function(i) {
    mean(data.rct[data.rct$subgroup == i & data.rct$Test == 1,][[outcome]]) - 
      mean(data.rct[data.rct$subgroup == i & data.rct$Test == 0,][[outcome]])
  })
  names(subgroupTaus.gold) <- sort(unique(data.rct$subgroup))
  
  # return the true causal effects
  return(list(data.os = data.os, 
              data.rct = data.rct,
              subgroupTaus.true = subgroupTaus.gold))
}





#' Stratified bootstrap sample
#'
#' Returns a vector of row indices for a stratified bootstrap sample, sampling
#' treated and control units separately within each stratum.
#'
#' @param sampleIds List of lists, one per stratum, each with elements `ids.0`
#'   (control row indices) and `ids.1` (treated row indices).
#' @param delta Numeric. Scaling factor for the bootstrap sample size within
#'   each stratum.
#'
#' @return Integer vector of sampled row indices.
#' @importFrom dqrng dqsample
stratifiedBootstrapSample <- function(sampleIds, delta) {
  
  ids <- unlist(l <- lapply(1:length(sampleIds), FUN = function(k) {
    s <- sampleIds[[k]]
    
    c(dqsample(s$ids.0, ceiling(delta*length(s$ids.0)), replace = TRUE), 
      dqsample(s$ids.1, ceiling(delta*length(s$ids.1)), replace = TRUE))
  }))
  
  return(ids)
}

