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

subgroupDefs_noSplit <- function(subgroupVars, data.os, data.rct) {
  
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
##         competitor estimator functions         ##
####################################################

# precision-weighted estimator
precistion.wtd <- function(ed, k, trueVar = NULL) {
  ed$rctVar/(ed$rctVar + ed$obsVar)*ed$obsEst + 
    ed$obsVar/(ed$rctVar + ed$obsVar)*ed$rctEst
}

# kappa_1 estimator (from Rosenman et al. 2023)
kappa.1 <- function(ed, k) {
  
  lambda.1 <- min(sum(ed$rctVar)/sum((ed$rctEst - ed$obsEst)^2), 1)
  
  estimator <- ed$obsEst + (1 - lambda.1) * (ed$rctEst - ed$obsEst)
  return(estimator)
}

# kappa_2 estimator (from Rosenman et al. 2023)
kappa.2 <- function(ed, k, trueVar = NULL) {
  
  lambda.2 <- pmin(pmax(sum(ed$rctVar)^2/
                          sum(ed$rctVar^2*(ed$rctEst - ed$obsEst)^2)*diag(ed$rctVar), 0), 1)
  
  estimator <- diag(lambda.2) * ed$obsEst + diag(1 - lambda.2) * ed$rctEst
  return(estimator)
}

# delta_1 estimator (from Green et al. 2005)
delta.1 <- function(ed, k) {
  
  lambda.1 <- max(1 - (k-2)/sum(1/ed$rctVar*(ed$rctEst - ed$obsEst)^2), 0)
  
  estimator <- ed$obsEst + (lambda.1) * (ed$rctEst - ed$obsEst)
  return(estimator)
}

# delta_2 estimator (from Green et al. 2005)
delta.2 <- function(ed, k) {
  
  lambda.2 <- pmax(diag(k) - (k-2)/sum(1/ed$rctVar^2*(ed$rctEst - ed$obsEst)^2)*
                     diag(1/ed$rctVar), 0)
  
  estimator <- ed$obsEst + lambda.2 %*% (ed$rctEst - ed$obsEst)
  return(estimator)
}

#####################################################
##            double shrinker functions            ##
#####################################################

# generic double shrinker function
shrinker <- function(gammaSq, etaSq, ed) {
  
  # unpack the parameters
  x <- ed$rctEst
  y <- ed$obsEst
  var.r <- ed$rctVar
  var.o <- ed$obsVar
  
  # determine the weights
  lambda <- (gammaSq + var.o)/(gammaSq + var.o + var.r)
  a <- etaSq*(gammaSq + var.r + var.o)/(var.r*(gammaSq + var.o) + etaSq*(gammaSq + var.r + var.o))
  
  # return the shrinker
  return(a*(lambda*x + (1 - lambda)*y))
}

# eb method of moments functions
eb.mm1 <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {

  # kappa and gamma estimation
  etaSq.pre <- max(mean(ed$rctEst^2 - ed$rctVar), 0)
  gammaSq.pre <- max(mean((ed$rctEst - ed$obsEst)^2 - ed$rctVar - ed$obsVar), 0)

  # kappa and gamma lower bounding (if requested)
  if(lowerBound) {
    etaSq <- max(etaSq.pre, 2*sum(ed$rctVar^2)/nrow(ed)/sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2*sum((ed$obsVar + ed$rctVar)^2)/
                     nrow(ed)/sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }

  # return the shrinker
  if(!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar)/(gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq*(gammaSq + ed$obsVar + ed$rctVar)/
                  (ed$rctVar*(gammaSq + ed$obsVar) + etaSq*(gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar*(gammaSq*(etaSq + 2*ed$obsVar) + gammaSq^2 + ed$obsVar^2)/
                  (etaSq*((gammaSq + ed$obsVar)^2 + ed$obsVar*ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}

eb.mm2 <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {
  
  # kappa and gamma estimation
  etaSq.pre <- max(mean(ed$rctEst^2 - ed$rctVar), 0)
  gammaSq.pre <- max(mean(-(ed$rctEst^2) + (ed$obsEst)^2 + ed$rctVar - ed$obsVar), 0)
  
  # kappa and gamma lower bounding (if requested)
  if(lowerBound) {
    etaSq <- max(etaSq.pre, 2*sum(ed$rctVar^2)/nrow(ed)/sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2*sum((ed$obsVar + ed$rctVar)^2)/
                     nrow(ed)/sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }
  
  # return the shrinker
  if(!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar)/(gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq*(gammaSq + ed$obsVar + ed$rctVar)/
                  (ed$rctVar*(gammaSq + ed$obsVar) + etaSq*(gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar*(gammaSq*(etaSq + 2*ed$obsVar) + gammaSq^2 + ed$obsVar^2)/
                  (etaSq*((gammaSq + ed$obsVar)^2 + ed$obsVar*ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}

# eb MLE function
eb.mle <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {
  
  # solve for the optimizing parameters 
  hyperParams.constrOptim <- tryCatch({constrOptim(c(1e-12, 1e-12), ebLik, grad = ebLikGrad, 
                                                   ui = diag(1, 2), ci = c(0, 0), ed = ed)$par},
                                      error = function(e) {c(0, 0)})
  hyperParams.constrOptim <- pmax(hyperParams.constrOptim, 0)
  
  hyperParams.optim <- optim(c(0, 0), ebLik, gr = ebLikGrad, lower = c(0, 0), 
                             method = 'L-BFGS-B', ed = ed)$par
  hyperParams.optim <- pmax(hyperParams.optim, 0)
  
  # choose hyperparameters based on whichever yields lower function value
  if(ebLik(hyperParams.constrOptim, ed) < ebLik(hyperParams.optim, ed)) {
    gammaSq.pre <- hyperParams.constrOptim[1]
    etaSq.pre <- hyperParams.constrOptim[2]
  } else {
    gammaSq.pre <- hyperParams.optim[1]
    etaSq.pre <- hyperParams.optim[2]
  }
  
  # kappa and gamma lower bounding (if requested)
  if(lowerBound) {
    etaSq <- max(etaSq.pre, 2*sum(ed$rctVar^2)/nrow(ed)/sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2*sum((ed$obsVar + ed$rctVar)^2)/
                     nrow(ed)/sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }
  
  # return the shrinker
  if(!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar)/(gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq*(gammaSq + ed$obsVar + ed$rctVar)/
                  (ed$rctVar*(gammaSq + ed$obsVar) + etaSq*(gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar*(gammaSq*(etaSq + 2*ed$obsVar) + gammaSq^2 + ed$obsVar^2)/
                  (etaSq*((gammaSq + ed$obsVar)^2 + ed$obsVar*ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}

# eb URE minimization function
eb.ure <- function(ed, k, lowerBound = FALSE, returnCIComps = FALSE) {
  
  # solve for the optimizing parameters 
  hyperParams.constrOptim <- tryCatch({constrOptim(c(1e-6, 1e-6), URE, grad = UREgrad, 
                     ui = diag(1, 2), ci = c(0, 0), ed = ed)$par},
                     error = function(e) {c(0, 0)})
  hyperParams.optim <- optim(c(0, 0), URE, lower = c(0, 0), method = 'L-BFGS-B', 
                             ed = ed)$par
  
  # choose hyperparameters based on whichever yields lower function value
  if(URE(hyperParams.constrOptim, ed) < URE(hyperParams.optim, ed)) {
    gammaSq.pre <- hyperParams.constrOptim[1]
    etaSq.pre <- hyperParams.constrOptim[2]
  } else {
    gammaSq.pre <- hyperParams.optim[1]
    etaSq.pre <- hyperParams.optim[2]
  }
  
  # kappa and gamma lower bounding (if requested)
  if(lowerBound) {
    etaSq <- max(etaSq.pre, 2*sum(ed$rctVar^2)/nrow(ed)/sum(ed$rctVar))
    gammaSq <- max(gammaSq.pre, 2*sum((ed$obsVar + ed$rctVar)^2)/
                     nrow(ed)/sum(ed$obsVar + ed$rctVar))
  } else {
    etaSq <- etaSq.pre
    gammaSq <- gammaSq.pre
  }

  # return the value
  if(!returnCIComps)
    shrinker(gammaSq, etaSq, ed)
  else {
    return(list(shrinker = shrinker(gammaSq, etaSq, ed),
                lambda = (gammaSq + ed$obsVar)/(gammaSq + ed$obsVar + ed$rctVar),
                a = etaSq*(gammaSq + ed$obsVar + ed$rctVar)/
                  (ed$rctVar*(gammaSq + ed$obsVar) + etaSq*(gammaSq + ed$obsVar + ed$rctVar)),
                c = ed$rctVar*(gammaSq*(etaSq + 2*ed$obsVar) + gammaSq^2 + ed$obsVar^2)/
                  (etaSq*((gammaSq + ed$obsVar)^2 + ed$obsVar*ed$rctVar)),
                etaSqLwr = etaSq == etaSq.pre,
                gammaSqLwr = gammaSq == gammaSq.pre))
  }
}

#######################################################
##       auxiliary functions for eb estimators       ##
#######################################################

# URE function
URE <- function(hyperparams, ed) {
  
  # unpack the hyperparameters
  gammaSq <- hyperparams[1] + 1e-12 # fudge factors to avoid errors
  etaSq <- hyperparams[2] + 1e-12 # fudge factors to avoid errors
  
  # unpack the parameters
  var.r <- ed$rctVar
  var.o <- ed$obsVar
  
  # determine the weights
  lambda <- (gammaSq + var.o)/(gammaSq + var.o + var.r)
  a <- etaSq*(gammaSq + var.r + var.o)/(var.r*(gammaSq + var.o) + etaSq*(gammaSq + var.r + var.o))
  
  # return the URE
  biasTerm <- sum((shrinker(gammaSq, etaSq, ed) - ed$rctEst)^2)
  varTerm <- -2*sum(var.r*(1 - a*lambda))
  ureVal <- (sum(var.r) + biasTerm + varTerm) 
  
  return(ureVal)
  
}

UREgrad <- function(hyperparams, ed) {
  
  # unpack the hyperparameters
  gammaSq <- hyperparams[1]
  etaSq <- hyperparams[2]
  
  # unpack the parameters
  var.r <- ed$rctVar
  var.o <- ed$obsVar
  
  # bias grad
  biasGrad <- c(sum(-((2*etaSq*ed$rctVar^2* ((etaSq + ed$rctVar) *ed$obsEst - 
                  etaSq*ed$rctEst)*(etaSq*ed$obsEst - (etaSq + gammaSq + 
                  ed$obsVar)* ed$rctEst))/((gammaSq + ed$obsVar)*ed$rctVar + 
                  etaSq*(gammaSq + ed$obsVar + ed$rctVar))^3)), 
                sum(-((2*(gammaSq + ed$obsVar)*ed$rctVar^2* (ed$rctVar*ed$obsEst + 
                  (gammaSq + ed$obsVar)*ed$rctEst)*(-etaSq* ed$obsEst + (etaSq + gammaSq + 
                  ed$obsVar)*ed$rctEst))/((gammaSq + ed$obsVar)*ed$rctVar + 
                  etaSq*(gammaSq + ed$obsVar + ed$rctVar))^3)))
  
  # var grad
  varGrad <- c(sum((2*etaSq^2*ed$rctVar^2)/((gammaSq + ed$obsVar)*ed$rctVar + 
                etaSq*(gammaSq + ed$obsVar + ed$rctVar))^2), 
               sum((2*(gammaSq + ed$obsVar)^2*ed$rctVar^2)/((gammaSq + ed$obsVar)*
                ed$rctVar + etaSq*(gammaSq + ed$obsVar + ed$rctVar))^2))
  
  biasGrad + varGrad
  
  
}
  
# eb likelihood function 
ebLik <- function(hyperParams, ed) {
  
  # unpack the parameters
  gammaSq <- hyperParams[1]
  etaSq <- hyperParams[2]
  
  # return the likelihood
  lik <- sum(log(etaSq + ed$rctVar) +
    ed$rctEst^2/(etaSq + ed$rctVar)
  + log(gammaSq + etaSq + ed$obsVar + 1e-12) +  # fudge factor to avoid errors
    ed$obsEst^2/(gammaSq + etaSq + ed$obsVar + 1e-12)) # fudge factor to avoid errors
  
  return(lik)
}

ebLikGrad <- function(hyperParams, ed) {
  
  # unpack the parameters
  gammaSq <- hyperParams[1]
  etaSq <- hyperParams[2]
  
  # return the gradient
  -c(sum(ed$obsEst^2/(gammaSq + etaSq + ed$obsVar)^2 - 
          1/(gammaSq + etaSq + ed$obsVar)), 
    sum(ed$rctEst^2/(etaSq + ed$rctVar)^2 - 
               1/(etaSq + ed$rctVar) +
             ed$obsEst^2/(gammaSq + etaSq + ed$obsVar)^2 - 
               1/(gammaSq + etaSq + ed$obsVar)))
}


####################################################
##         confidence interval functions          ##
####################################################

ebCiFunc <- function(ebData, ed, alpha = 0.05) {
  
  sapply(ebData$c, FUN = function(x) {cva(x, alpha = alpha)$cv})*
    ebData$a*sqrt(ebData$lambda^2*ed$rctVar + 
                           (1 - ebData$lambda)^2*ed$obsVar)
}

coverage <- function(ptEst, interval, trueVals) {
  ptEst - interval < trueVals &
    trueVals < ptEst + interval
}