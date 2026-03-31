############################################
##              MSE function              ##
############################################

# function to compute loss
loss <- function(tauTrue, tauHat) {
  mean((tauTrue - tauHat)^2)
}

#######################################################
##      function to get strata causal estimates      ##
#######################################################

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

###############################################
##       subgroup definition function        ##
###############################################

subgroupDefs_noSplit <- function(subgroupVars,
                                 data.os, data.rct,
                                 outcome,
                                 treatment = "Test") {
  
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
    mean(data.rct[data.rct$subgroup == i & data.rct[[treatment]] == 1,][[outcome]]) -
      mean(data.rct[data.rct$subgroup == i & data.rct[[treatment]] == 0,][[outcome]])
  })
  names(subgroupTaus.gold) <- sort(unique(data.rct$subgroup))
  
  # return the true causal effects
  return(list(data.os = data.os, 
              data.rct = data.rct,
              subgroupTaus.true = subgroupTaus.gold))
}

####################################################
##      stratified bootstrap sample function      ##
####################################################

# faster function to return RCT IDs under different sampling schemes
stratifiedBootstrapSample <- function(sampleIds, delta) {
  
  ids <- unlist(l <- lapply(1:length(sampleIds), FUN = function(k) {
    s <- sampleIds[[k]]
    
    c(dqsample(s$ids.0, ceiling(delta*length(s$ids.0)), replace = TRUE), 
      dqsample(s$ids.1, ceiling(delta*length(s$ids.1)), replace = TRUE))
  }))
  
  return(ids)
}

####################################################
##              CI coverage function              ##
####################################################

coverage <- function(ptEst, interval, trueVals) {
  ptEst - interval < trueVals &
    trueVals < ptEst + interval
}