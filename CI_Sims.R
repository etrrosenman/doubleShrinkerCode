rm(list = ls())

library(mvtnorm)
library(plyr)
library(dqrng)
library(rootSolve)
library(ebci)
set.seed(2023)

# load the helper functions 
setwd(here())
source("Helper Functions.R")

# load the synthetic data
load('synthetic_data.os.Rdata')
load('synthetic_data.rct.Rdata')

#################################################
####                constants                ####
#################################################

# control variables
outcome <- "Outcome.CHD"
numRctUnits <- 1000
numSamples <- 1000
alpha <- 0.05 

# subgroups
subgroupVarLists <- list("CVD", "AGER", "LANGLEYSCAT", c("AGER", "CVD"),
                         c("CVD", "LANGLEYSCAT"), c("AGER", "LANGLEYSCAT"), 
                         c("CVD", "AGER", "LANGLEYSCAT"))

#################################################
####                constants                ####
#################################################

# resample and compute MSE
results <- lapply(subgroupVarLists, FUN = function(subgroupVars) {
  
  ###############################################
  ##      stratum tagging & preliminaries      ##
  ###############################################
  
  # obtain the strata indicators for each unit in the obs & rct data 
  print(subgroupVars)
  allData <- subgroupDefs_noSplit(subgroupVars, data.os, data.rct) 
  obsData <- allData$data.os
  rctData <- allData$data.rct
  k <- length(allData$subgroupTaus.true)
  
  pe <- getParamEstimates(rctData, obsData, outcome)
  
  # pre-store id data
  obsSampleIds <- lapply(s <- sort(unique(obsData$subgroup)), FUN = function(k) {
    list(ids.0 = which(obsData$subgroup == k & obsData$Test == 0), ids.1 = which(obsData$subgroup == k & obsData$Test == 1))
  })
  names(obsSampleIds) <- s
  
  rctSampleIds <- lapply(s <- sort(unique(rctData$subgroup)), FUN = function(k) {
    list(ids.0 = which(rctData$subgroup == k & rctData$Test == 0), ids.1 = which(rctData$subgroup == k & rctData$Test == 1))
  })
  names(rctSampleIds) <- s
  
  ###############################################
  ##          bootstrap from the data          ##
  ###############################################

  coverageValues <- mclapply(1:numSamples, FUN = function(iter) {
    
    if(iter %% 200 == 0)
      system(sprintf('echo "%s"', iter))
    
    # extract bootstrap samples from observational and rct data
    obsSample <- obsData[stratifiedBootstrapSample(obsSampleIds, 1), 
                         which(names(obsData) %in% c(outcome, "propScore", "Test", "subgroup"))]
    rctSample <- rctData[stratifiedBootstrapSample(rctSampleIds, numRctUnits/nrow(rctData)), 
                         which(names(rctData) %in% c(outcome, "Test", "subgroup"))]
    
    # estimate all the necessary parameters
    estimatorData <- getParamEstimates(rctSample, obsSample, outcome)
    
    # extract the needed parameters for the double shrinker EBCI construction
    eb.mm1.params <- eb.mm1(estimatorData, k, lowerBound = TRUE, returnCIComps = TRUE)
    eb.mm2.params <- eb.mm2(estimatorData, k, lowerBound = TRUE, returnCIComps = TRUE)
    eb.mle.params <- eb.mle(estimatorData, k, lowerBound = TRUE, returnCIComps = TRUE)
    eb.ure.params <- eb.ure(estimatorData, k, lowerBound = TRUE, returnCIComps = TRUE)
    
    # construct the intervals
    intervals.mm1 <- ebCiFunc(eb.mm1.params, estimatorData)
    intervals.mm2 <- ebCiFunc(eb.mm2.params, estimatorData)
    intervals.mle <- ebCiFunc(eb.mle.params, estimatorData)
    intervals.ure <- ebCiFunc(eb.ure.params, estimatorData)

    # get the coverage
    z_alpha <- qnorm(1 - alpha/2)
    true <- allData$subgroupTaus.true
    
    list(baseCoverage = coverage(estimatorData$rctEst, z_alpha*sqrt(estimatorData$rctVar), true),
         mm1Coverage = coverage(eb.mm1.params$shrinker, intervals.mm1, true),
         mm2Coverage = coverage(eb.mm2.params$shrinker, intervals.mm2, true),
         mleCoverage = coverage(eb.mle.params$shrinker, intervals.mle, true),
         ureCoverage = coverage(eb.ure.params$shrinker, intervals.ure, true),
         intervals.base = z_alpha*sqrt(estimatorData$rctVar),
         intervals.mmOne = intervals.mm1,
         intervals.mmTwo = intervals.mm2,
         intervals.mle = intervals.mle,
         intervals.ure = intervals.ure)
  })
  
  # compute average and minum coverage rates 
  coverage <- rowMeans(sapply(coverageValues, FUN = function(x) {
    unlist(x[grepl("Coverage", names(x))])
  }))
  avgCoverage <- colMeans(m <- matrix(coverage, nrow = k))
  minCoverage <- apply(m, 2, min)
  names(avgCoverage) <- unique(gsub('\\..+', '', names(coverage)))
  names(minCoverage) <- names(avgCoverage)
  
  # compute interval lengths
  intervalLengths <- rowMeans(sapply(coverageValues, FUN = function(x) {
    unlist(x[grepl("intervals", names(x))])
  }))
  avgLength <- colMeans(matrix(intervalLengths, nrow = k))
  names(avgLength) <- unique(gsub('intervals\\.|[0-9]+', '', names(intervalLengths)))
  
  # return the average coverage statistics
  list(avgCoverage = avgCoverage, 
       minCoverage = minCoverage, 
       avgLength = avgLength)
})


#################################################
####          summarize results              ####
#################################################

# average coverage
avgCoverage <- tibble(data.frame(
    t(sapply(results, FUN = function(x) {x[[1]]}))
  ))
print(avgCoverage)

# min coverage
minCoverage <- tibble(data.frame(
    t(sapply(results, FUN = function(x) {x[[2]]}))
  ))
print(minCoverage) 

# average length
intervalLengths <- tibble(data.frame(
    t(sapply(results, FUN = function(x) {x[[3]]}))
  ))
print(intervalLengths)
