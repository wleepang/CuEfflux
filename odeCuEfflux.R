library(deSolve)
library(matlab)
source('CuEfflux_Func.R')

# set global values
source('CuEfflux_GlobalParams.R')
isKO0702 = F
isKO2581 = F
isOE0702 = F
isOE2581 = F

# initialize set model parameters and initial state
source('CuEfflux_Init.R')

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('CuEfflux_rxnDef.R')

nu = set.nu(rxn)
a  = set.a(rxn)

tic()
out = ode(x0, seq(0,18000,by=10), dxdt, list(nu=nu, a=a), method='daspk')
toc()

## plot mRNA dynamics of default system
t.abs = out[,'time'] / 60 # sim time is in seconds, convert to minutes

plot(t.abs, .T('Cu', out)/out[,'V'], type='l')

rna = cbind(out[, 'M0700'], out[, 'M0702'], out[, 'M2581'], out[, 'Mgfp'])/out[,'V']
matplot(t.abs, rna, type='l', col=c('black','blue', 'red', 'green'), lwd=3, lty=1, ann=F)
legend('topright', c('0700', '0702', '2581', 'gfp'), cex=0.8, lwd=3, col=c('black', 'blue', 'red', 'green'), bty='n')
title(xlab='Time (min)', ylab='mRNA Conc')

prot = cbind(.T('P0700', out), .T('P0702', out)/10, .T('P2581', out)/10, out[, 'Pgfp'])/out[,'V']
matplot(t.abs, prot, type='l', col=c('black', 'blue', 'red', 'green'), lwd=3, lty=1, ann=F)
legend('topright', c('0700', '0702 (x10)', '2581 (x10)', 'GFP'), cex=0.8, lwd=3, col=c('black', 'blue', 'red', 'green'), bty='n')
title(xlab='Time (min)', ylab='Protein Conc')

#t.abs[which(rna[,1] == max(rna[,1]))]
print(as.vector(.T('Cu', out)/out[,'V'])[length(t.abs)])
print(as.vector(.T('P0700', out)/out[,'V'])[length(t.abs)])

#save(out, nu, a, x0, file='wt_ts.RData')