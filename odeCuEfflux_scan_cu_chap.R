# library of MATLAB emulation functions - used to get tic() and toc()
library(matlab)
library(deSolve)
source('CuEfflux_Func.R')

################################################################################
################################################################################
# SCANS OF EXT. CU AND CHAPERONE ABUNDANCE

source('CuEfflux_GlobalParams.R')
source('CuEfflux_Init.R')
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
	a = paste('mu*(', paste(.T('P0702'), collapse='+'), ')', sep='', collapse=''), 
	a.s = FALSE, 
	nu = setStoic(x, 'P0702', +1))
##################################

nu = set.nu(rxn)
a  = set.a(rxn)

t.final = 18000
Levels.cu = 10^seq(2,7,length=25) # 1000 - 1000000 Cu.out molecules
Levels.ch = 10^seq(0,5,length=25) #    1 - 100000 chap molecules

scan.result.cuch.wt = model.scan.2d(c('Cu.out', 'P0702'), list(Cu.out = Levels.cu, P0702 = Levels.ch), 
																	 'wt_scan_CuOut_Chap', t.final=t.final)

################################################################################
# Plot scan results
srcuchw = scan.result.cuch.wt

###
# plot.cus = function(s, include=c('w', 'k', 'o', 'o2'), log='xy', ...) {
# 	X = NULL
# 	col = NULL
# 	if ('w' %in% include) {
# 		X = cbind(X, .T(s, srcuw))
# 		col = c(col, 'black')
# 	}
# 	if ('k' %in% include) {
# 		X = cbind(X, .T(s, srcuk))
# 		col = c(col, 'red')
# 	}
# 	if ('o' %in% include) {
# 		X = cbind(X, .T(s, srcuo))
# 		col = c(col, 'green')
# 	}
# 	if ('o2' %in% include) {
# 		X = cbind(X, .T(s, srcuo2))
# 		col = c(col, 'orange')
# 	}
# 	
# 	if (is.null(X)) stop('Nothing to plot!')
# 	
# 	X = X / srcuw[,'V']
# 	#X = cbind(.T(s, srqw), .T(s, srqo)) / srqw[,'V']
# 	matplot(srcuw[,'Level'], X, 
# 					log=log, 
# 					type='l', lty=1, lwd=3,
# 					col=col, 
# 					ann=F, ...)
# 	lines(c(5e5,5e5), c(1e-12, 1e12), col='#ffcc00', lty=2)
# 	title(xlab='Cu.out', ylab=s)
# }