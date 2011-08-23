# library of MATLAB emulation functions - used to get tic() and toc()
library(matlab)
library(deSolve)
source('CuEfflux_Func.R')

################################################################################
################################################################################
# SCANS OF CHAPERONE ABUNDANCE

source('CuEfflux_GlobalParams.R') # set global values
source('CuEfflux_Init.R')         # initialize set model parameters and initial state
source('CuEfflux_rxnDef.R')       # define reactions

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
	a = paste('mu*(', paste(.T('P0702'), collapse='+'), ')', sep='', collapse=''), 
	a.s = FALSE, 
	nu = setStoic(x, 'P0702', +1))
##################################

nu = set.nu(rxn)
a  = set.a(rxn)

# wt scan chap level
Levels = 10^seq(0,5,length=50) # 1 - 100000 Chaperone molecules
scan.result.chap.wt = model.scan('P0702', Levels, 'wt_scan_chap', t.final=18000)

################################################################################
# Plot scan results
srchw = scan.result.chap.wt

###
mai.orig = par('mai')
mai = mai.orig
mai[c(1,3)] = c(0.05, 1)
par(mai=mai)
layout(matrix(1:2, nrow=2, byrow=T))
plot(srchw[,'Level'],.T('Cu', srchw)/srchw[,'V'], 
		 log='x', xaxt='n', 
		 type='l', lty=1, lwd=3, col='black', ann=F)
title(ylab='[Cu]')

mai[c(1,3)] = c(1, 0.05)
par(mai=mai)
plot(srchw[,'Level'], srchw[,'P0700']/srchw[,'V']/10, 
		 log='x', 
		 type='l', lty=1, lwd=3, col='black', ann=F, xaxp=c(1,5,1))
title(xlab='Chaperone Abundance', ylab='[P0700]')

layout(1)
par(mai.orig)
