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

outcome       <- "Outcome.CHD"
numRctUnits   <- 1000
numSamples    <- 1000

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
  allData  <- subgroupDefs_noSplit(subgroupVars, data.os, data.rct, outcome)
  obsData  <- allData$data.os
  rctData  <- allData$data.rct
  k        <- length(allData$subgroupTaus.true)

  # pre-store IDs for treated & control units within each stratum
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
  estimatorValues <- mclapply(1:numSamples, FUN = function(iter) {

    if (iter %% 500 == 0) system(sprintf('echo "%s"', iter))

    obsSample <- obsData[stratifiedBootstrapSample(obsSampleIds, 1), ]
    rctSample <- rctData[stratifiedBootstrapSample(rctSampleIds,
                                                   numRctUnits / nrow(rctData)), ]

    ed <- getParamEstimates(rctSample, obsSample, outcome)

    ed$precisionWtd.est <- precistion.wtd(ed, k)
    ed$kappa1.plus.est  <- kappa.1(ed, k)
    ed$kappa2.plus.est  <- kappa.2(ed, k)
    ed$delta.1.est      <- delta.1(ed, k)
    ed$delta.2.est      <- delta.2(ed, k)
    ed$eb.mm1           <- eb.mm1(ed, k, lowerBound = TRUE)
    ed$eb.mm2           <- eb.mm2(ed, k, lowerBound = TRUE)
    ed$eb.mle           <- eb.mle(ed, k, lowerBound = TRUE)
    ed$eb.ure           <- eb.ure(ed, k, lowerBound = TRUE)

    list(estimators = ed[, -c(1, 3, 5)])
  })

  estimatorList <- lapply(estimatorValues, FUN = function(res) {res$estimators})

  mseVals   <- rowMeans(sapply(estimatorList, FUN = function(res) {
    apply(res, 2, FUN = function(x) {loss(x, allData$subgroupTaus.true)})
  }), na.rm = TRUE)
  biasSqVals <- apply(Reduce("+", estimatorList) / numSamples, 2,
                      FUN = function(x) {loss(x, allData$subgroupTaus.true)})
  varVals    <- (numSamples - 1) / numSamples *
    colMeans(matrix(apply(sapply(estimatorList, unlist), 1, var), nrow = k))
  names(varVals) <- names(biasSqVals)

  list(numStrata = k, mse = mseVals, bias = biasSqVals, var = varVals)
})

#################################################
####          summarize results              ####
#################################################

rawMse   <- data.frame(t(sapply(results, FUN = function(x) {c(numStrata = x$numStrata, x$mse)})))
mseTable <- tibble(
  data.frame(numStrata = rawMse$numStrata,
             apply(rawMse[, 3:ncol(rawMse)], 2,
                   FUN = function(x) {x / rawMse$rctEst}))
)
print(mseTable)
