# GillespieSSA implementation of Cu Efflux Model B
library(GillespieSSA)

# initialize set model parameters and initial state
source('CuEfflux_Init.R')

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('CuEfflux_rxnDef.R')

# Final time
tf <- 1000
simName <- "Cu Efflux Model B"
