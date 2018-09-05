SimulateRawCorrelation <- function(p, corrRange) {
  # Input
  # p          Number of variables
  # corrRange  list(min, max), the minimum and maximum amount of correlation between -1 and 1.
  R <- matrix(runif(n = p * p, min = corrRange$min, max = corrRange$max), ncol = p)       
  R <- R * lower.tri(R) + t(R * lower.tri(R))
  diag(R) <- 1
  return(R)  
}

NearestPosDef <- function(x, corr = T) {
  # Input
  # x	    numeric n * n approximately positive definite symmetric matrix, typically an approximation to a correlation            or covariance matrix.
  # corr  logical indicating if the matrix should be a correlation matrix.
  return(as.matrix(nearPD(x = x, corr = corr)$mat))
}

Cor2Cov <- function(cor.mat, sdRange) {
  # Input
  # cor.mat  correlation matrix
  # sdRange  The range for the standard deviation for each variable (diagonal of covariance mat)
  return(cor2cov(cor.mat, runif(n = ncol(cor.mat), 
                                min = sdRange$min, 
                                max = sdRange$max)))
}

SimulateCov <- function(p, corrRange, sdRange) {
  if(corrRange$min < -1.0 || corrRange$max > 1.0) {
    stop("Improper correlation range.")
  }
  
  if(sdRange$min < 0.0) {
    stop("Improper standard deviation range.")
  }
  
  if(p == 1) {
    return((runif(n = 1, min = sdRange$min, max = sdRange$max)) ^ 2)
  } else {
    rawCorr <- SimulateRawCorrelation(p = p, 
                                      corrRange = list(min = corrRange$min, 
                                                       max = corrRange$max)
    )
    corr <- NearestPosDef(x = rawCorr, corr = T)
    
    return(Cor2Cov(cor.mat = corr, sdRange = sdRange))  
  }
}

SimulateBlockCovariance <- function(totalVars, maxVarsInABlock, corrRange, sdRange) {

  numVars <- 1
  seedCovariance <- matrix(0, ncol = totalVars, nrow = totalVars)

  while(numVars <= totalVars) {
    
    p <- sample.int(n = maxVarsInABlock, size = 1)
    if((numVars + p) > totalVars) {
      seedCovariance[numVars:totalVars, numVars:totalVars] <- 
        SimulateCov(p = (totalVars - numVars + 1), 
                    corrRange = list(min = corrRange$min, max = corrRange$max), 
                    sdRange = list(min = sdRange$min, max = sdRange$max))
    } else {
      seedCovariance[numVars:(numVars + p - 1), numVars:(numVars + p - 1)] <- 
        SimulateCov(p = p, 
                    corrRange = list(min = corrRange$min, max = corrRange$max), 
                    sdRange = list(min = sdRange$min, max = sdRange$max))
    }
    numVars <- numVars + p
  }
  return(seedCovariance)  
}

ScaleDown <- function(mat, totalVars, scaleParam) {
  for(i in 1:nrow(mat)) {
    mat[i, ] <- mat[i, ] / (sum(abs(mat[i, ]) + scaleParam * totalVars))
  }
  return(mat)
}

SimulateRelationship <- function(totalVars, varRelationStructure, ...) {
  
  configs <- list(...)
  distr <- "uniform"
  mu <- 0
  sig <- 1
  r.min <- -0.5
  r.max <- 0.5
  scaleParam <- 5

  if(!is.null(configs$distr)) {
    distr <- configs$distr
  }
  if(!is.null(configs$mean)) {
    mu <- configs$mean
  }
  if(!is.null(configs$sd)) {
    sig <- configs$sd
  }
  if(!is.null(configs$r.min)) {
    r.min <- configs$r.min
  }
  if(!is.null(configs$sd)) {
    r.max <- configs$r.max
  }

  if(!is.null(configs$scaleParam)) {
    scaleParam <- configs$scaleParam
  }
  
  if(distr == "uniform") {
    out <- matrix(runif(n = totalVars * totalVars, 
                        min = r.min, 
                        max = r.max), 
                  ncol = totalVars, 
                  nrow = totalVars) * 
      varRelationStructure
  } else if (distr == "normal") {
    out <- matrix(rnorm(n = totalVars * totalVars, 
                        mean = mu, 
                        sd = sig), 
                  ncol = totalVars, 
                  nrow = totalVars) * 
      varRelationStructure
  }
  out <- ScaleDown(mat = out, totalVars = totalVars, scaleParam = scaleParam)
  return(out)
}

SimulateIntercept <- function(totalVars, allZero = T) {
  if(allZero) {
    out <- rep(0, totalVars)
  } else {
    out <- rnorm(n = totalVars, mean = 0, sd = 1)
  }
  return(VectorToColumn(out))
}

GenerateNoise <- function(covariance) {
  return(VectorToColumn(mvrnorm(n = 1, 
                                mu = rep(0, nrow(covariance)), 
                                Sigma = covariance)))
}

InitData <- function(totalVars, covariance, var.p = 1) {
  # Input
  # var.p  The order of Vector Auto Regression model.
  out <- matrix(0, nrow = totalVars, ncol = var.p)
  for(i in 1:var.p) {
    out[, i] <- mvrnorm(n = 1, mu = rep(0, totalVars), Sigma = covariance)
  }
  return(out)
}

SimulateMultipleRelationship <- function(var.p, covariance, gradualChange = T, drift = 0.1, scaleParam = 5) {
  # Input
  # var.p  The order of vector autoregressive regression model.
  # covariance  Covariance Matrix.
  # gradualChange  If we want the subsequent relationship matrix to be gradually drifting
  # drift  The drift parameter. The drift will be random N(0, drift).
  
  # Initialize
  totalVars <- nrow(covariance)
  varRelationStructure <- VariablesRelationshipStructure(covariance = covariance)
  relationships <- list(SimulateRelationship(totalVars = totalVars, 
                                             varRelationStructure = varRelationStructure,
                                             scaleParam = scaleParam))
  
  # Generate subsequent relationship matrices
  if(var.p > 1) {
    p <- 2
    while(p <= var.p) {
      
      if(gradualChange) {
        tmp <- relationships[[length(relationships)]] + 
          SimulateRelationship(totalVars = totalVars, 
                               varRelationStructure = varRelationStructure, 
                               mean = 0, 
                               sd = drift,
                               distr = "normal")
      } else {
        tmp <- SimulateRelationship(totalVars = totalVars, 
                                    varRelationStructure = varRelationStructure)
      }
      
      relationships[[length(relationships) + 1]] <- ScaleDown(mat = tmp, 
                                                              totalVars = totalVars, 
                                                              scaleParam = scaleParam)
      
      p <- p + 1
    }  
  }
  
  return(relationships)
}

SimVARseries <- function(numObservations, seedRelationships, seedCovariance) {
  totalVars <- nrow(seedCovariance)
  var.p <- length(seedRelationships)
  
  X <- matrix(0, nrow = totalVars, ncol = numObservations)
  seedIntercept <- SimulateIntercept(totalVars = totalVars,
                                     allZero = F)
  
  iter <- var.p + 1
  while(iter <= numObservations) {
    out <- seedIntercept
    for(i in 1:var.p) {
      out <- out + seedRelationships[[i]] %*% X[, iter - i]
    }
    X[, iter] <- out + GenerateNoise(covariance = seedCovariance)
    iter <- iter + 1
  }
  
  return(t(X[,-(1:var.p)]))
}
