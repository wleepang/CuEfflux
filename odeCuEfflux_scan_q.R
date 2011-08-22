# library of MATLAB emulation functions - used to get tic() and toc()
library(matlab)
library(deSolve)
source('CuEfflux_Func.R')

################################################################################
################################################################################
# SCANS OF Q CAPACITY

Levels = 10^seq(0,4,length=50) # 1 - 10000 Cu.out molecules
GENOME_COPY = 25

# wt scan q level
source('CuEfflux_GlobalParams.R')
kGenomeCopy = GENOME_COPY
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.q.wt = model.scan('Q', Levels, 'wt_scan_Q')

# ko scan q level
source('CuEfflux_GlobalParams.R')
kGenomeCopy = GENOME_COPY
isKO = T
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.q.ko = model.scan('Q', Levels, 'ko_scan_Q')

# oe scan q level
source('CuEfflux_GlobalParams.R')
kGenomeCopy = GENOME_COPY
isOE = T
source('CuEfflux_Init.R')
source('CuEfflux_rxnDef.R')
nu = set.nu(rxn)
a  = set.a(rxn)

scan.result.q.oe = model.scan('Q', Levels, 'oe_scan_Q')

################################################################################
# Plot scan results
srqw  = scan.result.q.wt
srqk  = scan.result.q.ko
srqo  = scan.result.q.oe

###
plot.qs = function(s, include=c('w', 'k', 'o')) {
	X = NULL
	col = NULL
	if ('w' %in% include) {
		X = cbind(X, .T(s, srqw))
		col = c(col, 'black')
	}
	if ('k' %in% include) {
		X = cbind(X, .T(s, srqk))
		col = c(col, 'red')
	}
	if ('o' %in% include) {
		X = cbind(X, .T(s, srqo))
		col = c(col, 'green')
	}
	
	if (is.null(X)) stop('Nothing to plot!')
	
	X = X / srqw[,'V']
	#X = cbind(.T(s, srqw), .T(s, srqo)) / srqw[,'V']
	matplot(srqw[,'Level'], X, 
					log='x', 
					type='l', lty=1, lwd=3,
					col=col, 
					ann=F)
	lines(c(600,600), c(1e-12, 1e12), col='#ffcc00', lty=2)
	lines(c(200,200), c(1e-12, 1e12), col='#ffcc00', lty=2)
	lines(c(2000,2000), c(1e-12, 1e12), col='#ffcc00', lty=2)
	title(xlab='Q', ylab=s)
}