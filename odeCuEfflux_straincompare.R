# compare strains
library(deSolve)
library(matlab)
source('CuEfflux_Func.R')

result = matrix(0, nrow=7, ncol=2)
colnames(result) = c('cu.ss', 'P0700.ss')
rownames(result) = c('wt', 'D0702', 'D2581', 'DD', 'O0702', 'O2581', 'OO')

Strain = data.frame(strain=c('wt', 'D0702', 'D2581', 'DD', 'O0702', 'O2581', 'OO'),
										KO0702=c(F,    T,       F,       T,    F,       F,       F),
										KO2581=c(F,    F,       T,       T,    F,       F,       F), 
										OE0702=c(F,    F,       F,       F,    T,       F,       T),
										OE2581=c(F,    F,       F,       F,    F,       T,       T))

trends = vector('list', 7)
names(trends)=rownames(result)
															
for (i in 1:dim(Strain)[1]) {
	s = Strain[i,]
	print(s)
	source('CuEfflux_GlobalParams.R')
	isKO0702 = s$KO0702
	isKO2581 = s$KO2581
	isOE0702 = s$OE0702
	isOE2581 = s$OE2581
	modelName = 'full'
	
	source('CuEfflux_Init.R')
	source('CuEfflux_rxnDef.R')
	
	nu = set.nu(rxn)
	a  = set.a(rxn)
	
	tic()
	trends[[as.character(s$strain)]] = ode(x0, seq(0,18000,by=100), dxdt, c(list(nu=nu, a=a), parms), method='daspk')
	toc()
}

# plot the result
pdf(file=sprintf('odeCuEfflux_straincompare_%s.pdf', modelName), width=11, height=8.5, paper='USr')
par(mfcol=c(3,7),mar=c(4,4,0.1,0.1),oma=c(0,0,0,0))
for (n in names(trends)) {
	print(n)
	
	out = trends[[n]]
	t.abs = out[,'time'] / 60 # sim time is in seconds, convert to minutes
	
	# Cu profile
	#print(par('mfg'))
	plot(t.abs, .T('Cu', out)/out[,'V'], type='l', ylim=c(0,500), ann=F)
	
	if (n == 'wt') {
		title(main=n, ylab='[Cu]')
	} else {
		title(main=n)
	}
	
	# mRNA profile
	#print(par('mfg'))
	rna = cbind(out[, 'M0700'], out[, 'M0702'], out[, 'M2581'], out[, 'Mgfp'])/out[,'V']
	matplot(t.abs, rna, type='l', col=c('black','blue', 'red', 'green'), lwd=3, lty=1, ann=F)
	legend('topright', c('0700', '0702', '2581', 'gfp'), cex=0.8, lwd=3, col=c('black', 'blue', 'red', 'green'), bty='n')
	
	if (n == 'wt') {
		title(ylab='[mRNA]')
	}
	
	# prot profile
	#print(par('mfg'))
	prot = cbind(.T('P0700', out), .T('P0702', out)/10, .T('P2581', out)/10, out[, 'Pgfp'])/out[,'V']
	matplot(t.abs, prot, type='l', col=c('black', 'blue', 'red', 'green'), lwd=3, lty=1, ann=F)
	legend('topright', c('0700', '0702 (x10)', '2581 (x10)', 'GFP'), cex=0.8, lwd=3, col=c('black', 'blue', 'red', 'green'), bty='n')
	
	if (n == 'wt') {
		title(xlab='Time (min)', ylab='[Protein]')
	} else {
		title(xlab='Time (min)')
	}
	
	cu.ss = as.vector(.T('Cu', out)/out[,'V'])[length(t.abs)]
	P0700.ss = as.vector(.T('P0700', out)/out[,'V'])[length(t.abs)]
	
	result[n,] = c(cu.ss, P0700.ss)
}															

par(mfcol=c(1,1), oma=c(1.5, 2, 1, 1))
barplot(result[,'cu.ss'], names.arg=rownames(result), ylab='[Cu]')
barplot(result[,'P0700.ss'], names.arg=rownames(result), ylab='[P0700]', log='y')

dev.off()