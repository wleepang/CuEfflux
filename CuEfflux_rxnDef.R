# this function helps to set the stoic matrix with a little error checking
setStoic = function(tpl, species, values) {
	if (length(species) != length(values)) stop('species names and stoic vector length mismatch')
	
	bSpeciesNotFound = !species %in% names(tpl)
	if (any(bSpeciesNotFound)) {
		stop(
			cat('species not defined: ', 
					paste(species[bSpeciesNotFound], collapse=', ')))
	}
	
	tpl[species] = values
	return(tpl)
}

# this is for building the stoic matrix see below in reaction defs
x <- rep(0, length(x0))
names(x) = names(x0)

# time is in units of seconds

# reaction definitions
rxn = list()
#rxn[['']] = list(
#	a = '', 										# propensity equation
# a.s = TRUE, 								# is a time dependent propensity?
#	nu = setStoic(x, c(), c()))

rxn[['Cu export']] = list(
	a = '0.1*P0700.Cu', 
	nu = setStoic(x, c('P0700.Cu', 'P0700', 'Cu.out'), c(-1, +1, +1)))

k.Cu.Import = function(Cu.out, Cu.T, Cu.PartCoef = 1000, k.Cu.Import.BaseRate = 1e-5) {
	grad = Cu.out - Cu.T*Cu.PartCoef
	rate = 0
	if (grad >= 0) {
		rate = grad*k.Cu.Import.BaseRate
	}
	
	return(rate)
}
rxn[['Cu import']] = list(
	a = paste('k.Cu.Import(Cu.out, ',
						paste(names(x)[(regexpr('Cu', names(x)) > 0) & (names(x) != 'Cu.out')], collapse='+'),
						')*Cu.out', sep="", collapse = ""), 
	a.s = FALSE, 
	nu = setStoic(x, c('Cu.out', 'Cu'), c(-1, +1)))

rxn[['D0700 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'D0700', 'P1179.Cu.D0700'), c(-1, -1, +1)))

rxn[['D0700 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700'), c(-1, +1, +1)))

rxn[['D0702 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D0702', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'D0702', 'P1179.Cu.D0702'), c(-1, -1, +1)))

rxn[['D0702 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D0702', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702'), c(-1, +1, +1)))

rxn[['Dgfp P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*Dgfp', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'Dgfp', 'P1179.Cu.Dgfp'), c(-1, -1, +1)))

rxn[['Dgfp P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.Dgfp', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp'), c(-1, +1, +1)))

rxn[['M0700 degradation']] = list(
	a = '0.00180197*M0700', 
	a.s = FALSE, 
	nu = setStoic(x, 'M0700', -1))

rxn[['M0700 transcription']] = list(
	a = '0.20718*P1179.Cu.D0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700', 'M0700'), c(-1, +1, +1, +1)))

rxn[['M0702 degradation']] = list(
	a = '0.00136725*M0702', 
	a.s = FALSE, 
	nu = setStoic(x, 'M0702', -1))

rxn[['M0702 transcription']] = list(
	a = '0.646571*P1179.Cu.D0702', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702', 'M0702'), c(-1, +1, +1, +1)))

rxn[['Mgfp degradation']] = list(
	a = '0.001*Mgfp', 
	a.s = FALSE, 
	nu = setStoic(x, 'Mgfp', -1))

rxn[['Mgfp transcription']] = list(
	a = '0.138889*P1179.Cu.Dgfp', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp', 'Mgfp'), c(-1, +1, +1, +1)))

rxn[['OE0702 transcription']] = list(
	a = '0.646571*OE0702', 
	a.s = FALSE, 
	nu = setStoic(x, 'M0702', +1))

rxn[['P0700.Cu by chap +1']] = list(
	a = '0.01*P0702.Cu*P0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu', 'P0700', 'P0702.Cu.P0700'), c(-1, -1, +1)))

rxn[['P0700.Cu by chap -1']] = list(
	a = '0.1*P0702.Cu.P0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702.Cu', 'P0700'), c(-1, +1, +1)))

rxn[['P0700.Cu by chap +2']] = list(
	a = '1.0*P0702.Cu.P0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702', 'P0700.Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu dissociation']] = list(
	a = '0.0001*P0700.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0700.Cu', 'P0700', 'Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu non-specific']] = list(
	a = '0.001*P0700*Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0700', 'Cu', 'P0700.Cu'), c(-1, -1, +1)))

rxn[['P0700 degradation']] = list(
	a = 'mu*P0700', 
	a.s = FALSE, 
	nu = setStoic(x, 'P0700', -1))

rxn[['P0700 translation']] = list(
	a = '0.002331*M0700', 
	a.s = FALSE, 
	nu = setStoic(x, 'P0700', +1))

rxn[['P0702.Cu bind']] = list(
	a = '0.01*P0702*Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702', 'Cu', 'P0702.Cu'), c(-1, -1, +1)))

rxn[['P0702.Cu by P1179 +1']] = list(
	a = '0.01*P0702*P1179.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702', 'P1179.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)))

rxn[['P0702.Cu by Q +1']] = list(
	a = '0.01*P0702*Q.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702', 'Q.Cu', 'P0702.Cu.Q'), c(-1, -1, +1)))

rxn[['P0702.Cu debind']] = list(
	a = '0.001*P0702.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu', 'P0702', 'Cu'), c(-1, +1, +1)))

rxn[['P0702.Cu degradation']] = list(
	a = '0.00015*P0702.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu', 'Cu'), c(-1, +1)))

rxn[['P0702 degradation']] = list(
	a = 'mu*P0702', 
	a.s = FALSE, 
	nu = setStoic(x, 'P0702', -1))

rxn[['P0702 translation']] = list(
	a = '0.0289855*M0702', 
	a.s = FALSE, 
	nu = setStoic(x, 'P0702', +1))

rxn[['P1179.Cu by chap +1']] = list(
	a = '0.01*P1179*P0702.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179', 'P0702.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)))

rxn[['P1179.Cu by chap -1']] = list(
	a = '0.1*P0702.Cu.P1179', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P0702.Cu', 'P1179'), c(-1, +1, +1)))

rxn[['P1179.Cu by chap 2']] = list(
	a = '0.1*P0702.Cu.P1179', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P1179.Cu', 'P0702'), c(-1, +1, +1)))

rxn[['P1179.Cu dissociation']] = list(
	a = '0.0001*P1179.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'P1179', 'Cu'), c(-1, +1, +1)))

rxn[['P1179.Cu non-specific']] = list(
	a = '0.001*P1179*Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179', 'Cu', 'P1179.Cu'), c(-1, -1, +1)))

rxn[['Pgfp degradation']] = list(
	a = 'mu*Pgfp', 
	a.s = FALSE, 
	nu = setStoic(x, 'Pgfp', -1))

rxn[['Pgfp translation']] = list(
	a = '0.008333*Mgfp', 
	a.s = FALSE, 
	nu = setStoic(x, 'Pgfp', +1))

rxn[['Q.Cu by chap +1']] = list(
	a = '0.01*P0702.Cu*Q', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu', 'Q', 'P0702.Cu.Q'), c(-1, -1, +1)))

rxn[['Q.Cu by chap -1']] = list(
	a = '0.1*P0702.Cu.Q', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702.Cu', 'Q'), c(-1, +1, +1)))

rxn[['Q.Cu by chap +2']] = list(
	a = '1.0*P0702.Cu.Q', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702', 'Q.Cu'), c(-1, +1, +1)))

rxn[['Q.Cu non-specific']] = list(
	a = '0.001*Q*Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('Q', 'Cu', 'Q.Cu'), c(-1, -1, +1)))

rxn[['Q.Cu turnover']] = list(
	a = 'mu*Q.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('Q.Cu', 'Q', 'Cu'), c(-1, +1, +1)))

rxn[['V growth']] = list(
	a = 'mu*V',
	a.s = FALSE,
	nu = setStoic(x, 'V', +1))

# State-change matrix
# rows = species, cols = reactions
nu = NULL
for (r in rxn) nu = cbind(nu, r$nu)
colnames(nu) = names(rxn)

# Propensity vector
a = NULL
for (r in rxn) a = c(a, r$a)
names(a) = names(rxn)