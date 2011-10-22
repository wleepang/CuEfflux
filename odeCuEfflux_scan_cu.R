# library of MATLAB emulation functions - used to get tic() and toc()
library(matlab)
library(deSolve)
source('CuEfflux_Func.R')

################################################################################
################################################################################
# SCANS OF EXT. CU

t.final = 18000
Levels = 10^seq(2,7,length=25) # 1000 - 1000000 Cu.out molecules

# wt scan cu.out level
source('CuEfflux_GlobalParams.R')
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cu.wt = model.scan('Cu.out', Levels, 'wt_scan_CuOut', t.final=t.final, t.step=100)

# ko scan cu.out level
source('CuEfflux_GlobalParams.R')
isKO0702 = TRUE
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cu.ko = model.scan('Cu.out', Levels, 'ko_scan_CuOut', t.final=t.final, t.step=100)

# oe scan cu.out level
source('CuEfflux_GlobalParams.R')
isOE0702 = TRUE
parms[['k.OE0702.transcription']] = 0.646571
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cu.oe = model.scan('Cu.out', Levels, 'oe_scan_CuOut', t.final=t.final, t.step=100)

# adjusted oe scan cu.out level
source('CuEfflux_GlobalParams.R')
isOE0702 = TRUE
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')

nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.cu.oe2 = model.scan('Cu.out', Levels, 'oe2_scan_CuOut', t.final=t.final, t.step=100)

################################################################################
# Plot scan results
srcuw = scan.result.cu.wt
srcuk = scan.result.cu.ko
srcuo = scan.result.cu.oe
srcuo2= scan.result.cu.oe2

###
plot.cus = function(s, include=c('w', 'k', 'o', 'o2'), log='xy', ...) {
	X = NULL
	col = NULL
	if ('w' %in% include) {
		X = cbind(X, .T(s, srcuw))
		col = c(col, 'black')
	}
	if ('k' %in% include) {
		X = cbind(X, .T(s, srcuk))
		col = c(col, 'red')
	}
	if ('o' %in% include) {
		X = cbind(X, .T(s, srcuo))
		col = c(col, 'green')
	}
	if ('o2' %in% include) {
		X = cbind(X, .T(s, srcuo2))
		col = c(col, 'orange')
	}
	
	if (is.null(X)) stop('Nothing to plot!')
	
	X = X / srcuw[,'V']
	#X = cbind(.T(s, srqw), .T(s, srqo)) / srqw[,'V']
	matplot(srcuw[,'Level'], X, 
					log=log, 
					type='l', lty=1, lwd=3,
					col=col, 
					ann=F, ...)
	lines(c(5e5,5e5), c(1e-12, 1e12), col='#ffcc00', lty=2)
	title(xlab='Cu.out', ylab=s)
}