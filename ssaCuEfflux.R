# GillespieSSA implementation of Cu Efflux Model B
library(GillespieSSA)

# Define parameters
#parms <- c(c1=1.0, c2=0.002, c3=0.5, c4=0.04)

kGenomeCopy = 15
isKO = FALSE
isOE = FALSE

# Initial state vector
x0 <- c(Cu							= 0, 
				Cu.out					= 500000,
				D0700           = kGenomeCopy,
				D0702           = kGenomeCopy*2*(!isKO),
				Dgfp            = kGenomeCopy,
				M0700           = 0,
				M0702           = 0,
				Mgfp            = 0,
				OE0702          = kGenomeCopy*isOE,
				P0700           = 0,
				P0700.Cu        = 0,
				P0702           = 700*(!isKO),
				P0702.Cu        = 0,
				P0702.Cu.P0700  = 0,
				P0702.Cu.P1179  = 0,
				P0702.Cu.Q      = 0,
				P1179           = kGenomeCopy*2,
				P1179.Cu        = 0,
				P1179.Cu.D0700  = 0,
				P1179.Cu.D0702  = 0,
				P1179.Cu.Dgfp   = 0,
				Pgfp            = 0,
				Q               = 100,
				Q.Cu            = 0)

# this is for building the stoic matrix see below in reaction defs
x <- rep(0, length(x0))
names(x) = names(x0)

setStoic = function(tpl, species, values) {
	tpl[species] = values
	return(tpl)
}

# reaction definitions
rxn = list()
#rxn[['']] = list(
#	a = '', 
#	nu = setStoic(x, c(), c()))

rxn[['Cu export']] = list(
	a = '0.1*P0700.Cu', 
	nu = setStoic(x, c('P0700.Cu', 'P0700', 'Cu.out'), c(-1, +1, +1)))

rxn[['Cu import']] = list(
	a = 'k.Cu.Import*Cu.out', 
	nu = setStoic(x, c('Cu.out', 'Cu'), c(-1, +1)))

rxn[['D0700 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D0700', 
	nu = setStoic(x, c('P1179.Cu', 'D0700', 'P1179.Cu.D0700'), c(-1, -1, +1)))

rxn[['D0700 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D0700', 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700'), c(-1, +1, +1)))

rxn[['D0702 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D0702', 
	nu = setStoic(x, c('P1179.Cu', 'D0702', 'P1179.Cu.D0702'), c(-1, -1, +1)))

rxn[['D0702 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D0702', 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702'), c(-1, +1, +1)))

rxn[['Dgfp P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*Dgfp', 
	nu = setStoic(x, c('P1179.Cu', 'Dgfp', 'P1179.Cu.Dgfp'), c(-1, -1, +1)))

rxn[['Dgfp P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.Dgfp', 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp'), c(-1, +1, +1)))

rxn[['M0700 degradation']] = list(
	a = '0.00180197*M0700', 
	nu = setStoic(x, 'M0700', -1))

rxn[['M0700 transcription']] = list(
	a = '0.20718*P1179.Cu.D0700', 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700', 'M0700'), c(-1, +1, +1, +1)))

rxn[['M0702 degradation']] = list(
	a = '0.00136725*M0702', 
	nu = setStoic(x, 'M0702', -1))

rxn[['M0702 transcription']] = list(
	a = '0.646571*P1179.Cu.D0702', 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702', 'M0702'), c(-1, +1, +1, +1)))

rxn[['Mgfp degradation']] = list(
	a = '0.001*Mgfp', 
	nu = setStoic(x, 'Mgfp', -1))

rxn[['Mgfp transcription']] = list(
	a = '0.138889*P1179.Cu.Dgfp', 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp', 'Mgfp'), c(-1, +1, +1, +1)))

rxn[['OE0702 transcription']] = list(
	a = '0.646571*OE0702', 
	nu = setStoic(x, 'M0702', +1))

rxn[['P0700.Cu by chap +1']] = list(
	a = '0.01*P0702.Cu*P0700', 
	nu = setStoic(x, c('P0702.Cu', 'P0700', 'P0702.Cu.P0700'), c(-1, -1, +1)))

rxn[['P0700.Cu by chap -1']] = list(
	a = '0.1*P0702.Cu.P0700', 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702.Cu', 'P0700'), c(-1, +1, +1)))

rxn[['P0700.Cu by chap +2']] = list(
	a = '1.0*P0702.Cu.P0700', 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702', 'P0700.Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu dissociation']] = list(
	a = '0.0001*P0700.Cu', 
	nu = setStoic(x, c('P0700.Cu', 'P0700', 'Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu non-specific']] = list(
	a = '0.001*P0700*Cu', 
	nu = setStoic(x, c('P0700', 'Cu', 'P0700.Cu'), c(-1, -1, +1)))

rxn[['P0700 degradation']] = list(
	a = 'mu*P0700', 
	nu = setStoic(x, 'P0700', -1))

rxn[['P0700 translation']] = list(
	a = '0.002331*M0700', 
	nu = setStoic(x, 'P0700', +1))

rxn[['P0702.Cu bind']] = list(
	a = '0.01*P0702*Cu', 
	nu = setStoic(x, c('P0702', 'Cu', 'P0702.Cu'), c(-1, -1, +1)))

rxn[['P0702.Cu by P1179 +1']] = list(
	a = '0.01*P0702*P1179.Cu', 
	nu = setStoic(x, c('P0702', 'P1179.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)))

rxn[['P0702.Cu by Q +1']] = list(
	a = '0.01*P0702*Q.Cu', 
	nu = setStoic(x, c('P0702', 'Q.Cu', 'P0702.Cu.Q'), c(-1, -1, +1)))

rxn[['P0702.Cu debind']] = list(
	a = '0.001*P0702.Cu', 
	nu = setStoic(x, c('P0702.Cu', 'P0702', 'Cu'), c(-1, +1, +1)))

rxn[['P0702.Cu degradation']] = list(
	a = '0.00015*P0702.Cu', 
	nu = setStoic(x, c('P0702.Cu', 'Cu'), c(-1, +1)))

rxn[['P0702 degradation']] = list(
	a = 'mu*P0702', 
	nu = setStoic(x, 'P0702', -1))

rxn[['P0702 translation']] = list(
	a = '0.0289855*M0702', 
	nu = setStoic(x, 'P0702', +1))

rxn[['P1179.Cu by chap +1']] = list(
	a = '0.01*P1179*P0702.Cu', 
	nu = setStoic(x, c('P1179', 'P0702.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)))

rxn[['P1179.Cu by chap -1']] = list(
	a = '0.1*P0702.Cu.P1179', 
	nu = setStoic(x, c('P1179.Cu.P0702', 'P0702.Cu', 'P1179'), c(-1, +1, +1)))

rxn[['P1179.Cu by chap 2']] = list(
	a = '0.1*P0702.Cu.P1179', 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P1179.Cu', 'P0702'), c(-1, +1, +1)))

rxn[['P1179.Cu dissociation']] = list(
	a = '0.0001*P1179.Cu', 
	nu = setStoic(x, c('P1179.Cu', 'P1179', 'Cu'), c(-1, +1, +1)))

rxn[['P1179.Cu non-specific']] = list(
	a = '0.001*P1179*Cu', 
	nu = setStoic(x, c('P1179', 'Cu', 'P1179.Cu'), c(-1, -1, +1)))

rxn[['Pgfp degradation']] = list(
	a = 'mu*Pgfp', 
	nu = setStoic(x, 'Pgfp', -1))

rxn[['Pgfp translation']] = list(
	a = '0.008333*Mgfp', 
	nu = setStoic(x, 'Pgfp', +1))

rxn[['Q.Cu by chap +1']] = list(
	a = '0.01*P0702.Cu*Q', 
	nu = setStoic(x, c('P0702.Cu', 'Q', 'P0702.Cu.Q'), c(-1, -1, +1)))

rxn[['Q.Cu by chap -1']] = list(
	a = '0.1*P0702.Cu.Q', 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702.Cu', 'Q'), c(-1, +1, +1)))

rxn[['Q.Cu by chap +2']] = list(
	a = '1.0*P0702.Cu.Q', 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702', 'Q.Cu'), c(-1, +1, +1)))

rxn[['Q.Cu non-specific']] = list(
	a = '0.001*Q*Cu', 
	nu = setStoic(x, c('Q', 'Cu', 'Q.Cu'), c(-1, -1, +1)))

rxn[['Q.Cu turnover']] = list(
	a = 'mu*Q.Cu', 
	nu = setStoic(x, c('Q.Cu', 'Q', 'Cu'), c(-1, +1, +1)))

# State-change matrix
# rows = species, cols = reactions
#nu <- matrix(c(-1, -2, +2,  0,
#                0, +1, -1, -1,
#                0,  0,  0, +1),
#                nrow=3,byrow=TRUE)
nu = NULL
for (r in rxn) nu = cbind(nu, r$nu)
colnames(nu) = names(rxn)

# Propensity vector
#a  <- c("c1*s1", 
#				"c2*s1*s1", 
#				"c3*s2", 
#				"c4*s2")
a = NULL
for (r in rxn) a = c(a, r$a)
names(a) = names(rxn)

# Final time
tf <- 10                                      
simName <- "Cu Efflux Model B"
