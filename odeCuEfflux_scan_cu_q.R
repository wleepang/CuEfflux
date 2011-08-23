# library of MATLAB emulation functions - used to get tic() and toc()
library(matlab)
library(deSolve)
source('CuEfflux_Func.R')

################################################################################
################################################################################
# SCANS OF EXT. CU AND Q

t.final = 18000
Levels.cu = 10^seq(2,7,length=25) # 1000 - 1000000 Cu.out molecules
Levels.q  = 10^seq(0,4,length=25) # 1 - 10000 Q molecules

# wt scan cu.out level
source('CuEfflux_GlobalParams.R')
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cuq.wt = model.scan.2d(c('Cu.out', 'Q'), list(Cu.out = Levels.cu, Q = Levels.q), 
																	 'wt_scan_CuOut_Q', t.final=t.final)

# ko scan cu.out level
source('CuEfflux_GlobalParams.R')
isKO = TRUE
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cuq.ko = model.scan.2d(c('Cu.out', 'Q'), list(Cu.out = Levels.cu, Q = Levels.q), 
																	 'ko_scan_CuOut_Q', t.final=t.final)

# oe scan cu.out level
source('CuEfflux_GlobalParams.R')
isOE = TRUE
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cuq.oe = model.scan.2d(c('Cu.out', 'Q'), list(Cu.out = Levels.cu, Q = Levels.q), 
																	 'oe_scan_CuOut_Q', t.final=t.final)

# adjusted oe scan cu.out level
source('CuEfflux_GlobalParams.R')
isOE = TRUE
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')

rxn[['OE0702 transcription']] = list(
	a = '0.646571*OE0702/100', # in order to get OE to be like expt. divide by factor of 100
	a.s = FALSE, 
	nu = setStoic(x, 'M0702', +1))

nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cuq.oe2 = model.scan.2d(c('Cu.out', 'Q'), list(Cu.out = Levels.cu, Q = Levels.q), 
																	 'oe2_scan_CuOut_Q', t.final=t.final)

################################################################################
# Plot scan results
srcuqw = scan.result.cuq.wt
srcuqk = scan.result.cuq.ko
srcuqo = scan.result.cuq.oe
srcuqo2= scan.result.cuq.oe2

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