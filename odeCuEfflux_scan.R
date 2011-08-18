# library of MATLAB emulation functions - used to get tic() and toc()
library(matlab)
library(deSolve)
source('CuEfflux_Func.R')

# set global values
source('CuEfflux_GlobalParams.R')

# initialize set model parameters and initial state
source('CuEfflux_Init.R')

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('CuEfflux_rxnDef.R')

## scanning code
# scan chaperone abundance
# 1) remove chaperone dynamic expression / degradation
rxn[['D0702 P1179.Cu binding']] = NULL
rxn[['D0702 P1179.Cu dissociation']] = NULL
rxn[['P0702.Cu degradation']] = NULL
rxn[['P0702 degradation']] = NULL

# 2) scale chaperone level with growth rate
##################################
## MAINTAIN CHAPERONE ABUNDANCE
##################################
rxn[['P0702 maintenance']] = list(
	a = 'mu', 
	a.s = FALSE, 
	nu = setStoic(x, 'P0702', +1))
##################################

nu = set.nu(rxn)
a  = set.a(rxn)

# wt scan chap level
#Levels = 10^seq(0,5,length=100)
#scan.result.wt = model.scan('P0702', Levels, 'wt_scan')

Levels = 10^seq(3,6,length=50)

# wt scan cu.out level
scan.result.wt = model.scan('Cu.out', Levels, 'wt_scan_CuOut')

# ko scan cu.out level
source('CuEfflux_GlobalParams.R')
isKO = TRUE
source('CuEfflux_Init.R')
scan.result.ko = model.scan('Cu.out', Levels, 'ko_scan_CuOut')

# oe scan cu.out level
source('CuEfflux_GlobalParams.R')
isOE = TRUE
source('CuEfflux_Init.R')
scan.result.oe = model.scan('Cu.out', Levels, 'oe_scan_CuOut')