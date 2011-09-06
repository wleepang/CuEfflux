## Reaction definitions for CuEfflux model
## used for both SSA and ODE simulation

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
	nu = setStoic(x, c('P0700.Cu', 'P0700'), c(-1, +1)))

rxn[['Cu import']] = list(
	a = paste('k.Cu.Import(Cu.out, ',
						paste(names(x)[(regexpr('Cu', names(x)) > 0) & (names(x) != 'Cu.out')], collapse='+'),
						')*V', sep="", collapse = ""), 
	a.s = FALSE, 
	nu = setStoic(x, c('Cu'), c(+1)))

rxn[['D0700 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D0700/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'D0700', 'P1179.Cu.D0700'), c(-1, -1, +1)))

rxn[['D0700 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700'), c(-1, +1, +1)))

rxn[['D0702 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D0702/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'D0702', 'P1179.Cu.D0702'), c(-1, -1, +1)))

rxn[['D0702 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D0702', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702'), c(-1, +1, +1)))

rxn[['D2581 P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*D2581/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'D2581', 'P1179.Cu.D2581'), c(-1, -1, +1)))

rxn[['D2581 P1179.Cu dissociation']] = list(
	a = '0.1*P1179.Cu.D2581', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D2581', 'P1179.Cu', 'D2581'), c(-1, +1, +1)))

rxn[['Dgfp P1179.Cu binding']] = list(
	a = '0.01*P1179.Cu*Dgfp/V^2', 
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
	a = '0.207918*P1179.Cu.D0700', 
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

rxn[['M2581 degradation']] = list(
	a = '0.003262*M2581', 
	a.s = FALSE, 
	nu = setStoic(x, 'M2581', -1))

rxn[['M2581 transcription']] = list(
	a = '0.424167*P1179.Cu.D2581', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.D2581', 'P1179.Cu', 'D2581', 'M2581'), c(-1, +1, +1, +1)))

rxn[['Mgfp degradation']] = list(
	a = '0.001*Mgfp', 
	a.s = FALSE, 
	nu = setStoic(x, 'Mgfp', -1))

rxn[['Mgfp transcription']] = list(
	a = '0.138889*P1179.Cu.Dgfp', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp', 'Mgfp'), c(-1, +1, +1, +1)))

rxn[['OE0702 transcription']] = list(
	a = '0.646571*OE0702/100', # in order to get OE to be like expt. divide by factor of 100
	a.s = FALSE, 
	nu = setStoic(x, 'M0702', +1))

rxn[['OE2581 transcription']] = list(
	a = '0.424167*OE2581/100', # in order to get OE to be like expt. divide by factor of 100
	a.s = FALSE, 
	nu = setStoic(x, 'M2581', +1))

rxn[['P0700.Cu by 0702 +1']] = list(
	a = '0.01*P0702.Cu*P0700/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu', 'P0700', 'P0702.Cu.P0700'), c(-1, -1, +1)))

rxn[['P0700.Cu by 0702 -1']] = list(
	a = '0.1*P0702.Cu.P0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702.Cu', 'P0700'), c(-1, +1, +1)))

rxn[['P0700.Cu by 0702 +2']] = list(
	a = '1.0*P0702.Cu.P0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702', 'P0700.Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu by 2581 +1']] = list(
	a = '0.01*P2581.Cu*P0700/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu', 'P0700', 'P2581.Cu.P0700'), c(-1, -1, +1)))

rxn[['P0700.Cu by 2581 -1']] = list(
	a = '0.1*P2581.Cu.P0700', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu.P0700', 'P2581.Cu', 'P0700'), c(-1, +1, +1)))

rxn[['P0700.Cu by 2581 +2']] = list(
	a = '0.1*P2581.Cu.P0700', # 10-fold less efficient than 0702
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu.P0700', 'P2581', 'P0700.Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu dissociation']] = list(
	a = '0.0001*P0700.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0700.Cu', 'P0700', 'Cu'), c(-1, +1, +1)))

rxn[['P0700.Cu non-specific']] = list(
	a = '0.001*P0700*Cu/V^2', 
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
	a = '0.01*P0702*Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702', 'Cu', 'P0702.Cu'), c(-1, -1, +1)))

rxn[['P0702.Cu by P1179 +1']] = list(
	a = '0.01*P0702*P1179.Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702', 'P1179.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)))

# rxn[['P0702.Cu by P2581 +1']] = list(
# 	a = '0.01*P0702*P2581.Cu/V^2', 
# 	a.s = FALSE, 
# 	nu = setStoic(x, c('P0702', 'P2581.Cu', 'P0702.Cu.P2581'), c(-1, -1, +1)))

rxn[['P0702.Cu by Q +1']] = list(
	a = '0.01*P0702*Q.Cu/V^2', 
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

rxn[['P2581.Cu bind']] = list(
	a = '0.01*P2581*Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581', 'Cu', 'P2581.Cu'), c(-1, -1, +1)))

rxn[['P2581.Cu by P1179 +1']] = list(
	a = '0.01*P2581*P1179.Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581', 'P1179.Cu', 'P2581.Cu.P1179'), c(-1, -1, +1)))

# rxn[['P2581.Cu by P0702 +1']] = list(
# 	a = '0.01*P2581*P0702.Cu/V^2', 
# 	a.s = FALSE, 
# 	nu = setStoic(x, c('P2581', 'P0702.Cu', 'P0702.Cu.P2581'), c(-1, -1, +1)))

rxn[['P2581.Cu by Q +1']] = list(
	a = '0.01*P2581*Q.Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581', 'Q.Cu', 'P2581.Cu.Q'), c(-1, -1, +1)))

rxn[['P2581.Cu debind']] = list(
	a = '0.001*P2581.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu', 'P2581', 'Cu'), c(-1, +1, +1)))

rxn[['P2581.Cu degradation']] = list(
	a = '0.00015*P2581.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu', 'Cu'), c(-1, +1)))

rxn[['P2581 degradation']] = list(
	a = 'mu*P2581', 
	a.s = FALSE, 
	nu = setStoic(x, 'P2581', -1))

rxn[['P2581 translation']] = list(
	a = '0.030769231*M2581', 
	a.s = FALSE, 
	nu = setStoic(x, 'P2581', +1))

rxn[['P0702.Cu by P2581 -1']] = list(
	a = '0.1*P0702.Cu.P2581',
	a.s = FALSE,
	nu = setStoic(x, c('P0702.Cu.P2581', 'P0702', 'P2581.Cu'), c(-1, +1, +1)))

rxn[['P2581.Cu by P0702 -1']] = list(
	a = '0.1*P0702.Cu.P2581',
	a.s = FALSE,
	nu = setStoic(x, c('P0702.Cu.P2581', 'P0702.Cu', 'P2581'), c(-1, +1, +1)))

rxn[['P1179.Cu by 0702 +1']] = list(
	a = '0.01*P1179*P0702.Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179', 'P0702.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)))

rxn[['P1179.Cu by 0702 -1']] = list(
	a = '0.1*P0702.Cu.P1179', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P0702.Cu', 'P1179'), c(-1, +1, +1)))

rxn[['P1179.Cu by 0702 +2']] = list(
	a = '0.1*P0702.Cu.P1179', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P1179.Cu', 'P0702'), c(-1, +1, +1)))

rxn[['P1179.Cu by 2581 +1']] = list(
	a = '0.01*P1179*P2581.Cu/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179', 'P2581.Cu', 'P2581.Cu.P1179'), c(-1, -1, +1)))

rxn[['P1179.Cu by 2581 -1']] = list(
	a = '0.1*P2581.Cu.P1179', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu.P1179', 'P2581.Cu', 'P1179'), c(-1, +1, +1)))

rxn[['P1179.Cu by 2581 +2']] = list(
	a = '0.01*P2581.Cu.P1179', # 10-fold less effecient than 0702
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu.P1179', 'P1179.Cu', 'P2581'), c(-1, +1, +1)))

rxn[['P1179.Cu dissociation']] = list(
	a = '0.0001*P1179.Cu', 
	a.s = FALSE, 
	nu = setStoic(x, c('P1179.Cu', 'P1179', 'Cu'), c(-1, +1, +1)))

rxn[['P1179.Cu non-specific']] = list(
	a = '0.001*P1179*Cu/V^2', 
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

rxn[['Q.Cu by 0702 +1']] = list(
	a = '0.01*P0702.Cu*Q/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu', 'Q', 'P0702.Cu.Q'), c(-1, -1, +1)))

rxn[['Q.Cu by 0702 -1']] = list(
	a = '0.1*P0702.Cu.Q', 
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702.Cu', 'Q'), c(-1, +1, +1)))

rxn[['Q.Cu by 0702 +2']] = list(
	a = '0.1*P0702.Cu.Q', # 10-fold less efficient than 2581
	a.s = FALSE, 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702', 'Q.Cu'), c(-1, +1, +1)))

rxn[['Q.Cu by 2581 +1']] = list(
	a = '0.01*P2581.Cu*Q/V^2', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu', 'Q', 'P2581.Cu.Q'), c(-1, -1, +1)))

rxn[['Q.Cu by 2581 -1']] = list(
	a = '0.1*P2581.Cu.Q', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu.Q', 'P2581.Cu', 'Q'), c(-1, +1, +1)))

rxn[['Q.Cu by 2581 +2']] = list(
	a = '1.0*P2581.Cu.Q', 
	a.s = FALSE, 
	nu = setStoic(x, c('P2581.Cu.Q', 'P2581', 'Q.Cu'), c(-1, +1, +1)))

rxn[['Q.Cu non-specific']] = list(
	a = '0.001*Q*Cu/V^2', 
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

rxn[['Q growth']] = list(
	a = paste('mu*(', paste(.T('Q'), collapse='+'), ')', sep='', collapse=''),
	a.s = FALSE,
	nu = setStoic(x, 'Q', +1))

rxn[['P1179 growth']] = list(
	a = paste('mu*(', paste(.T('P1179'), collapse='+'), ')', sep='', collapse=''),
	a.s = FALSE,
	nu = setStoic(x, 'P1179', +1))

rxn[['D0700 growth']] = list(
	a = paste('mu*(', paste(.T('D0700'), collapse='+'), ')', sep='', collapse=''),
	a.s = FALSE,
	nu = setStoic(x, 'D0700', +1))

rxn[['D0702 growth']] = list(
	a = paste('mu*(', paste(.T('D0702'), collapse='+'), ')', sep='', collapse=''),
	a.s = FALSE,
	nu = setStoic(x, 'D0702', +1))

rxn[['D2581 growth']] = list(
	a = paste('mu*(', paste(.T('D2581'), collapse='+'), ')', sep='', collapse=''),
	a.s = FALSE,
	nu = setStoic(x, 'D2581', +1))

rxn[['Dgfp growth']] = list(
	a = paste('mu*(', paste(.T('Dgfp'), collapse='+'), ')', sep='', collapse=''),
	a.s = FALSE,
	nu = setStoic(x, 'Dgfp', +1))
