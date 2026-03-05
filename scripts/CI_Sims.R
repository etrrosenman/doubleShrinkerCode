library(doubleShrinker)
library(mvtnorm)
library(parallel)
library(tibble)
library(here)
set.seed(2023)

# load the synthetic population data
load(here("data", "synthetic_data.os.RData"))
load(here("data", "synthetic_data.rct.RData"))

#################################################
####                constants                ####
#################################################

outcome     <- "Outcome.CHD"
numRctUnits <- 1000
numSamples  <- 1000
alpha       <- 0.05

subgroupVarLists <- list(
  "CVD", "AGER", "LANGLEYSCAT",
  c("AGER", "CVD"),
  c("CVD", "LANGLEYSCAT"),
  c("AGER", "LANGLEYSCAT"),
  c("CVD", "AGER", "LANGLEYSCAT")
)

#################################################
####      iterate through the subgroups      ####
#################################################

results <- lapply(subgroupVarLists, FUN = function(subgroupVars) {

  print(subgroupVars)
  allData <- subgroupDefs_noSplit(subgroupVars, data.os, data.rct, outcome)
  obsData <- allData$data.os
  rctData <- allData$data.rct
  k       <- length(allData$subgroupTaus.true)

  obsSampleIds <- lapply(s <- sort(unique(obsData$subgroup)), FUN = function(g) {
    list(ids.0 = which(obsData$subgroup == g & obsData$Test == 0),
         ids.1 = which(obsData$subgroup == g & obsData$Test == 1))
  })
  names(obsSampleIds) <- s

  rctSampleIds <- lapply(s <- sort(unique(rctData$subgroup)), FUN = function(g) {
    list(ids.0 = which(rctData$subgroup == g & rctData$Test == 0),
         ids.1 = which(rctData$subgroup == g & rctData$Test == 1))
  })
  names(rctSampleIds) <- s

  ###############################################
  ##          bootstrap from the data          ##
  ###############################################

  coverageValues <- mclapply(1:numSamples, FUN = function(iter) {

    if (iter %% 200 == 0) system(sprintf('echo "%s"', iter))

    obsSample <- obsData[stratifiedBootstrapSample(obsSampleIds, 1),
                         which(names(obsData) %in% c(outcome, "propScore", "Test", "subgroup"))]
    rctSample <- rctData[stratifiedBootstrapSample(rctSampleIds, numRctUnits / nrow(rctData)),
                         which(names(rctData) %in% c(outcome, "Test", "subgroup"))]

    ed <- getParamEstimates(rctSample, obsSample, outcome)

    mm1.params <- eb.mm1(ed, k, lowerBound = TRUE, returnCIComps = TRUE)
    mm2.params <- eb.mm2(ed, k, lowerBound = TRUE, returnCIComps = TRUE)
    mle.params <- eb.mle(ed, k, lowerBound = TRUE, returnCIComps = TRUE)
    ure.params <- eb.ure(ed, k, lowerBound = TRUE, returnCIComps = TRUE)

    intervals.mm1 <- ebCiFunc(mm1.params, ed)
    intervals.mm2 <- ebCiFunc(mm2.params, ed)
    intervals.mle <- ebCiFunc(mle.params, ed)
    intervals.ure <- ebCiFunc(ure.params, ed)

    z_alpha <- qnorm(1 - alpha / 2)
    true    <- allData$subgroupTaus.true

    list(
      baseCoverage = coverage(ed$rctEst, z_alpha * sqrt(ed$rctVar), true),
      mm1Coverage  = coverage(mm1.params$shrinker, intervals.mm1, true),
      mm2Coverage  = coverage(mm2.params$shrinker, intervals.mm2, true),
      mleCoverage  = coverage(mle.params$shrinker, intervals.mle, true),
      ureCoverage  = coverage(ure.params$shrinker, intervals.ure, true),
      intervals.base   = z_alpha * sqrt(ed$rctVar),
      intervals.mmOne  = intervals.mm1,
      intervals.mmTwo  = intervals.mm2,
      intervals.mle    = intervals.mle,
      intervals.ure    = intervals.ure
    )
  })

  coverage_rates <- rowMeans(sapply(coverageValues, FUN = function(x) {
    unlist(x[grepl("Coverage", names(x))])
  }))
  avgCoverage <- colMeans(m <- matrix(coverage_rates, nrow = k))
  minCoverage <- apply(m, 2, min)
  names(avgCoverage) <- unique(gsub("\\..+", "", names(coverage_rates)))
  names(minCoverage) <- names(avgCoverage)

  intervalLengths <- rowMeans(sapply(coverageValues, FUN = function(x) {
    unlist(x[grepl("intervals", names(x))])
  }))
  avgLength <- colMeans(matrix(intervalLengths, nrow = k))
  names(avgLength) <- unique(gsub("intervals\\.|[0-9]+", "", names(intervalLengths)))

  list(avgCoverage = avgCoverage, minCoverage = minCoverage, avgLength = avgLength)
})

#################################################
####          summarize results              ####
#################################################

avgCoverage <- tibble(data.frame(t(sapply(results, FUN = function(x) {x[[1]]}))))
print(avgCoverage)

minCoverage <- tibble(data.frame(t(sapply(results, FUN = function(x) {x[[2]]}))))
print(minCoverage)

intervalLengths <- tibble(data.frame(t(sapply(results, FUN = function(x) {x[[3]]}))))
print(intervalLengths)
