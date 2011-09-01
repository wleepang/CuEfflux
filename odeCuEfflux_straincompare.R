# compare strains
library(deSolve)
library(matlab)
source('CuEfflux_Func.R')

par(mfcol=c(3,7),mai=c(0.3,0.3,0.1,0.1))
result = NULL

Strain = data.frame(strain=c('wt', 'D0702', 'D2581', 'DD', 'O0702', 'O2581', 'OO'),
										KO0702=c(F, T, F, T, F, F, F),
										KO2581=c(F, F, T, T, F, F, F), 
										OE0702=c(F, F, F, T, F, F, T),
										OE2581=c(F, F, F, F, T, F, T))

for (i in 1:dim(Strain)[1]) {
	s = Strain[i,]
	
	source('CuEfflux_GlobalParams.R')
	isKO0702 = s$KO0702
	isKO2581 = s$KO2581
	isOE0702 = s$OE0702
	isOE2581 = s$OE2581
	
	source('CuEfflux_Init.R')
	source('CuEfflux_rxnDef.R')
	
	nu = set.nu(rxn)
	a  = set.a(rxn)
	
	tic()
	out = ode(x0, seq(0,18000,by=100), dxdt, list(nu=nu, a=a), method='daspk')
	toc()
	
	t.abs = out[,'time'] / 60 # sim time is in seconds, convert to minutes
	
	# Cu profile
	plot(t.abs, .T('Cu', out)/out[,'V'], type='l', ylim=c(0,500))
	
	# mRNA profile
	rna = cbind(out[, 'M0700'], out[, 'M0702'], out[, 'M2581'], out[, 'Mgfp'])/out[,'V']
	matplot(t.abs, rna, type='l', col=c('black','blue', 'red', 'green'), lwd=3, lty=1, ann=F)
	legend('topright', c('0700', '0702', '2581', 'gfp'), cex=0.8, lwd=3, col=c('black', 'blue', 'red', 'green'), bty='n')
	title(xlab='Time (min)', ylab='mRNA Conc')
	
	# prot profile
	prot = cbind(.T('P0700', out), .T('P0702', out)/10, .T('P2581', out)/10, out[, 'Pgfp'])/out[,'V']
	matplot(t.abs, prot, type='l', col=c('black', 'blue', 'red', 'green'), lwd=3, lty=1, ann=F)
	legend('topright', c('0700', '0702 (x10)', '2581 (x10)', 'GFP'), cex=0.8, lwd=3, col=c('black', 'blue', 'red', 'green'), bty='n')
	title(xlab='Time (min)', ylab='Protein Conc')
	
	cu.ss = as.vector(.T('Cu', out)/out[,'V'])[length(t.abs)]
	P0700.ss = as.vector(.T('P0700', out)/out[,'V'])[length(t.abs)]
	
	result = cbind(result, list(strain=s$strain, cu.ss=cu.ss, P0700.ss=P0700.ss))
}

# plot the result
par(mfcol=c(1,1), mai=c(0.8, 0.8, 0.2, 0.2))
