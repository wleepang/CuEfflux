library(deSolve)
source('CuEfflux_Func.R')

# set global values
source('CuEfflux_GlobalParams.R')

# initialize set model parameters and initial state
source('CuEfflux_Init.R')

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('CuEfflux_rxnDef.R')

nu = set.nu(rxn)
a  = set.a(rxn)

out = ode(x0, seq(0,36000,by=1), dxdt, list(nu=nu, a=a), method='bdf', verbose=T)

## plot mRNA dynamics of default system
t.abs = out[,'time'] / 60 # sim time is in seconds, convert to minutes

rna = cbind(out[, 'M0700'], out[, 'M0702'], out[, 'Mgfp'])/out[,'V']
matplot(t.abs, rna, type='l', col=c('blue', 'red', 'green'), lwd=3, lty=1, xlim=c(0.006,600), ann=F)
legend('topright', c('0700', 'Chaperone', 'gfp'), cex=0.8, lwd=3, col=c('blue', 'red', 'green'), bty='n')
title(xlab='Time (min)', ylab='# mRNA per cell')

prot = cbind(.T('P0700', out), .T('P0702', out)/10, out[, 'Pgfp'])/out[,'V']
matplot(t.abs, prot, type='l', col=c('blue', 'red', 'green'), lwd=3, lty=1, xlim=c(0.006,600), ann=F)
legend('topright', c('0700', 'Chaperone (x10)', 'GFP'), cex=0.8, lwd=3, col=c('blue', 'red', 'green'), bty='n')
title(xlab='Time (min)', ylab='# Proteins per cell')
