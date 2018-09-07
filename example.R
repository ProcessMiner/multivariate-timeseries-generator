rm(list = ls())
source("methods/libraries.R")
source("methods/datasimulation.R")
source("methods/utilities.R")

## Generate covariance matrix
seedCovariance <- SimulateBlockCovariance(totalVars = 6, 
                                          maxVarsInABlock = 3, 
                                          corrRange = list(min = -0.5, max = 0.9), 
                                          sdRange = list(min = 0.5, max = 0.75))

seedCovariance
print(PlotMatrix(seedCovariance))

## Generate relationship matrices, Ai
### Independent Ai
seedRelationships <- SimulateMultipleRelationship(var.p = 3, 
                                                  covariance = seedCovariance,
                                                  gradualChange = F,
                                                  scaleParam = 0.1)

multiplot(plotlist = list(PlotMatrix(mat = seedRelationships[[1]], color = "blue"), 
                          PlotMatrix(mat = seedRelationships[[2]], color = "blue"), 
                          PlotMatrix(mat = seedRelationships[[3]], color = "blue")), 
          cols = 3)

### Co-dependent Ai
seedRelationships <- SimulateMultipleRelationship(var.p = 3, 
                                                  covariance = seedCovariance,
                                                  gradualChange = T,
                                                  drift = 0.1,
                                                  scaleParam = 0.1)

multiplot(plotlist = list(PlotMatrix(mat = seedRelationships[[1]], limits = c(-0.05,0.1), color = "blue"), 
                          PlotMatrix(mat = seedRelationships[[2]], limits = c(-0.05,0.1), color = "blue"), 
                          PlotMatrix(mat = seedRelationships[[3]], limits = c(-0.05,0.1), color = "blue")), 
          cols = 3)


## Simulate multivariate time series
X <- SimVARseries(numObservations = 200, 
                  seedRelationships = seedRelationships, 
                  seedCovariance = seedCovariance)

plot.ts(X) # Plot the series
