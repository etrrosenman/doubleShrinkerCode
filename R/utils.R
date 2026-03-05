############################################
##              MSE function              ##
############################################

# function to compute loss
loss <- function(tauTrue, tauHat) {
  mean((tauTrue - tauHat)^2)
}



#######################################################
#       auxiliary functions for eb estimators   ----
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



