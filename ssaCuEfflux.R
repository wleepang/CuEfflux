# GillespieSSA implementation of Cu Efflux Model B
library(GillespieSSA)

# Define parameters
#parms <- c(c1=1.0, c2=0.002, c3=0.5, c4=0.04)

kGenomeCopy = 15
isKO = FALSE
isOE = FALSE
mu = 0 # cell growth rate 1/s

# Initial state vector
x0 <- c(Cu							= 0, 
				Cu.out					= 500000,
				D0700           = kGenomeCopy,
				D0702           = kGenomeCopy*2*(!isKO),
				Dgfp            = kGenomeCopy,
				M0700           = 0,
				M0702           = 0,
				Mgfp            = 0,
				OE0702          = kGenomeCopy*isOE,
				P0700           = 0,
				P0700.Cu        = 0,
				P0702           = 700*(!isKO),
				P0702.Cu        = 0,
				P0702.Cu.P0700  = 0,
				P0702.Cu.P1179  = 0,
				P0702.Cu.Q      = 0,
				P1179           = kGenomeCopy*2,
				P1179.Cu        = 0,
				P1179.Cu.D0700  = 0,
				P1179.Cu.D0702  = 0,
				P1179.Cu.Dgfp   = 0,
				Pgfp            = 0,
				Q               = 100,
				Q.Cu            = 0, 
				V               = 1000)

# this is for building the stoic matrix see below in reaction defs
x <- rep(0, length(x0))
names(x) = names(x0)

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('ssaCuEfflux_rxnDef.R')

# State-change matrix
# rows = species, cols = reactions
nu = NULL
for (r in rxn) nu = cbind(nu, r$nu)
colnames(nu) = names(rxn)

# Propensity vector
a = NULL
for (r in rxn) a = c(a, r$a)
names(a) = names(rxn)

# Final time
tf <- 1000
simName <- "Cu Efflux Model B"
