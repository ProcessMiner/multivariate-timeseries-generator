rm(list = ls())
source("shiny/functions/libraries.R")
source("shiny/functions/datasimulation.R")
source("shiny/functions/utilities.R")

seedCovariance <- SimulateBlockCovariance(totalVars = 6, 
                                          maxVarsInABlock = 3, 
                                          corrRange = list(min = -0.5, max = 0.9), 
                                          sdRange = list(min = 0.5, max = 0.75))

seedRelationships <- SimulateMultipleRelationship(var.p = 3, 
                                                  covariance = seedCovariance,
                                                  gradualChange = T,
                                                  drift = 0.1,
                                                  scaleParam = 0.1)

X <- SimVARseries(numObservations = 200, 
                  seedRelationships = seedRelationships, 
                  seedCovariance = seedCovariance)

plot.ts(X) # Plot the series

library(gridplot)
grow_layout(a, b, c)
multiplot(plotlist = list(PlotMatrix(mat = seedRelationships[[1]], limits = c(-0.1,0.1)), 
                          PlotMatrix(mat = seedRelationships[[2]], limits = c(-0.1,0.1)), 
                          PlotMatrix(mat = seedRelationships[[3]], limits = c(-0.1,0.1))), 
          cols = 3)

                                     