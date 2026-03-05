rm(list = ls())

library(mvtnorm)
library(plyr)
library(dqrng)
library(rootSolve)
library(here)
library(parallel)
set.seed(2023)

# load the helper functions 
setwd(here::here())
library( doubleShrinker )

# load the synthetic data
load('data/synthetic_data.os.Rdata')
load('data/synthetic_data.rct.Rdata')


#################################################
####                helper functions                ####
#################################################

# function to compute loss
loss <- function(tauTrue, tauHat) {
  mean((tauTrue - tauHat)^2)
}


#################################################
####                constants                ####
#################################################

# control variables
outcome <- "Outcome.CHD"
numRctUnits <- 1000
numSamples <- 1000

# subgroups
subgroupVarLists <- list("CVD", "AGER", "LANGLEYSCAT", c("AGER", "CVD"),
                         c("CVD", "LANGLEYSCAT"), c("AGER", "LANGLEYSCAT"), 
                         c("CVD", "AGER", "LANGLEYSCAT"))

#################################################
####      iterate through the subgroups      ####
#################################################

results <- lapply(subgroupVarLists, FUN = function(subgroupVars) {
  
  ###############################################
  ##      stratum tagging & preliminaries      ##
  ###############################################
  
  # obtain the strata indicators for each unit in the obs & rct data 
  print(subgroupVars)
  allData <- subgroupDefs_noSplit(subgroupVars, data.os, data.rct, 
                                  outcome=outcome) 
  obsData <- allData$data.os
  rctData <- allData$data.rct
  k <- length(allData$subgroupTaus.true)
  
  pe <- getParamEstimates(rctData, obsData, outcome)

  # pre-store IDs for treated & control units within each dataset
  obsSampleIds <- lapply(s <- sort(unique(obsData$subgroup)), FUN = function(k) {
    list(ids.0 = which(obsData$subgroup == k & obsData$Test == 0), 
         ids.1 = which(obsData$subgroup == k & obsData$Test == 1))
  })
  names(obsSampleIds) <- s
  
  rctSampleIds <- lapply(s <- sort(unique(rctData$subgroup)), FUN = function(k) {
    list(ids.0 = which(rctData$subgroup == k & rctData$Test == 0), 
         ids.1 = which(rctData$subgroup == k & rctData$Test == 1))
  })
  names(rctSampleIds) <- s
  
  ###############################################
  ##          bootstrap from the data          ##
  ###############################################
  estimatorValues <- mclapply(1:numSamples, FUN = function(iter) {
    
    # report progress
    if(iter %% 500 == 0)
      system(sprintf('echo "%s"', iter))
    
    # extract stratified bootstrap samples from obs and rct data
    obsSample <- obsData[stratifiedBootstrapSample(obsSampleIds, 1),]
    rctSample <- rctData[stratifiedBootstrapSample(rctSampleIds, 
                                                   numRctUnits/nrow(rctData)),]
    
    # estimate all the necessary parameters
    estimatorData <- getParamEstimates(rctSample, obsSample, outcome)
    
    # compute precision-weighted estimators
    estimatorData$precisionWtd.est <- precistion.wtd(estimatorData, k)
      
    # compute kappa estimators (from Rosenman et al. 2023)
    estimatorData$kappa1.plus.est <- kappa.1(estimatorData, k)
    estimatorData$kappa2.plus.est <- kappa.2(estimatorData, k)
    
    # compute delta estimators (from Green et al. 2005)
    estimatorData$delta.1.est <- delta.1(estimatorData, k)
    estimatorData$delta.2.est <- delta.2(estimatorData, k)
    
    # compute double shrinkage estimators (from manuscript)
    estimatorData$eb.mm1 <- eb.mm1(estimatorData, k, lowerBound = TRUE)
    estimatorData$eb.mm2 <- eb.mm2(estimatorData, k, lowerBound = TRUE)
    estimatorData$eb.mle <- eb.mle(estimatorData, k, lowerBound = TRUE)
    estimatorData$eb.ure <- eb.ure(estimatorData, k, lowerBound = TRUE)
    
    list(estimators = estimatorData[, -c(1, 3, 5)])
  })
  
  # reshape and return the results
  estimatorList <- lapply(estimatorValues, 
                          FUN = function(res) {res$estimators})
  
  # compute the MSE, bias squared, and variance for each estimator
  mseVals <- rowMeans(sapply(estimatorList, FUN = function(res) {
    apply(res, 2, FUN = function(x) {loss(x, allData$subgroupTaus.true)})
  }), na.rm = TRUE)
  biasSqVals <- apply(Reduce('+', estimatorList)/numSamples, 2, 
                      FUN = function(x) {loss(x, allData$subgroupTaus.true)})
  varVals <- (numSamples - 1)/numSamples*
    colMeans(matrix(apply(sapply(estimatorList, unlist), 1, var), nrow = k))
  names(varVals) <- names(biasSqVals)
  
  # return everything 
  list(numStrata = k, mse = mseVals, 
       bias = biasSqVals, var = varVals)
  
})

#################################################
####          summarize results              ####
#################################################

# mse table
rawMse <- data.frame(t(sapply(results, FUN = function(x) {c(numStrata = x$numStrata, x$mse)})))
mseTable <- tibble(
  data.frame(numStrata = rawMse$numStrata, apply(rawMse[,3:ncol(rawMse)], 2, 
                                                           FUN = function(x) {x/rawMse$rctEst}))
)
print(mseTable)
