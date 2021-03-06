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
#	nu = setStoic(x, c(), c()))

rxn[['Cu export']] = list(
	a = 'k.Cu.export*P0700.Cu', 
	nu = setStoic(x, c('P0700.Cu', 'P0700'), c(-1, +1)),
  mods= NULL)

rxn[['Cu import']] = list(
# 	a = paste('k.Cu.Import(Cu.out, ',
# 						paste(names(x)[(regexpr('Cu', names(x)) > 0) & (names(x) != 'Cu.out')], collapse='+'),
# 						', Cu.PartCoef, k.Cu.Import.BaseRate)*V', sep="", collapse = ""), 
	a = paste('k.Cu.Import.BaseRate*V*Cu.out', paste('k.Cu.Import.BaseRate*V*Cu.PartCoef*', .T('Cu'), sep='', collapse='-'), sep='-', collapse='-'),
  nu = setStoic(x, c('Cu'), c(+1)),
	mods= c('Cu.out', 'V', .T('Cu')[which(.T('Cu')!='Cu')]))

rxn[['D0700 P1179.Cu binding']] = list(
	a = 'k.D0700.P1179.Cu.binding*P1179.Cu*D0700/V^2', 
	nu = setStoic(x, c('P1179.Cu', 'D0700', 'P1179.Cu.D0700'), c(-1, -1, +1)),
	mods= 'V')

rxn[['D0700 P1179.Cu dissociation']] = list(
	a = 'k.D0700.P1179.Cu.dissociation*P1179.Cu.D0700', 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700'), c(-1, +1, +1)),
	mods= NULL)

rxn[['D0702 P1179.Cu binding']] = list(
	a = 'k.D0702.P1179.Cu.binding*P1179.Cu*D0702/V^2', 
	nu = setStoic(x, c('P1179.Cu', 'D0702', 'P1179.Cu.D0702'), c(-1, -1, +1)),
	mods= 'V')

rxn[['D0702 P1179.Cu dissociation']] = list(
	a = 'k.D0702.P1179.Cu.dissociation*P1179.Cu.D0702', 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702'), c(-1, +1, +1)),
	mods= NULL)

rxn[['D2581 P1179.Cu binding']] = list(
	a = 'k.D2581.P1179.Cu.binding*P1179.Cu*D2581/V^2', 
	nu = setStoic(x, c('P1179.Cu', 'D2581', 'P1179.Cu.D2581'), c(-1, -1, +1)),
	mods= 'V')

rxn[['D2581 P1179.Cu dissociation']] = list(
	a = 'k.D2581.P1179.Cu.dissociation*P1179.Cu.D2581', 
	nu = setStoic(x, c('P1179.Cu.D2581', 'P1179.Cu', 'D2581'), c(-1, +1, +1)),
	mods= NULL)

rxn[['Dgfp P1179.Cu binding']] = list(
	a = 'k.Dgfp.P1179.Cu.binding*P1179.Cu*Dgfp/V^2', 
	nu = setStoic(x, c('P1179.Cu', 'Dgfp', 'P1179.Cu.Dgfp'), c(-1, -1, +1)),
	mods= 'V')

rxn[['Dgfp P1179.Cu dissociation']] = list(
	a = 'k.Dgfp.P1179.Cu.dissociation*P1179.Cu.Dgfp', 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp'), c(-1, +1, +1)),
	mods= NULL)

rxn[['M0700 degradation']] = list(
	a = 'k.M0700.degradation*M0700', 
	nu = setStoic(x, 'M0700', -1),
	mods= NULL)

rxn[['M0700 transcription']] = list(
	a = 'k.M0700.transcription*P1179.Cu.D0700', 
	nu = setStoic(x, c('P1179.Cu.D0700', 'P1179.Cu', 'D0700', 'M0700'), c(-1, +1, +1, +1)),
	mods= NULL)

rxn[['M0702 degradation']] = list(
	a = 'k.M0702.degradation*M0702', 
	nu = setStoic(x, 'M0702', -1),
	mods= NULL)

rxn[['M0702 transcription']] = list(
	a = 'k.M0702.transcription*P1179.Cu.D0702', 
	nu = setStoic(x, c('P1179.Cu.D0702', 'P1179.Cu', 'D0702', 'M0702'), c(-1, +1, +1, +1)),
	mods= NULL)

rxn[['M2581 degradation']] = list(
	a = 'k.M2581.degradation*M2581', 
	nu = setStoic(x, 'M2581', -1),
	mods= NULL)

rxn[['M2581 transcription']] = list(
	a = 'k.M2581.transcription*P1179.Cu.D2581', 
	nu = setStoic(x, c('P1179.Cu.D2581', 'P1179.Cu', 'D2581', 'M2581'), c(-1, +1, +1, +1)),
	mods= NULL)

rxn[['Mgfp degradation']] = list(
	a = 'k.Mgfp.degradation*Mgfp', 
	nu = setStoic(x, 'Mgfp', -1),
	mods= NULL)

rxn[['Mgfp transcription']] = list(
	a = 'k.Mgfp.transcription*P1179.Cu.Dgfp', 
	nu = setStoic(x, c('P1179.Cu.Dgfp', 'P1179.Cu', 'Dgfp', 'Mgfp'), c(-1, +1, +1, +1)),
	mods= NULL)

rxn[['OE0702 transcription']] = list(
	a = 'k.OE0702.transcription*OE0702', 
	nu = setStoic(x, 'M0702', +1),
	mods= 'OE0702')

rxn[['OE2581 transcription']] = list(
	a = 'k.OE2581.transcription*OE2581', 
	nu = setStoic(x, 'M2581', +1),
	mods= 'OE2581')

rxn[['P0700.Cu by P0702 F1']] = list(
	a = 'k.P0700.Cu.by.P0702.F1*P0702.Cu*P0700/V^2', 
	nu = setStoic(x, c('P0702.Cu', 'P0700', 'P0702.Cu.P0700'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P0700.Cu by P0702 R1']] = list(
	a = 'k.P0700.Cu.by.P0702.R1*P0702.Cu.P0700', 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702.Cu', 'P0700'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0700.Cu by P0702 F2']] = list(
	a = 'k.P0700.Cu.by.P0702.F2*P0702.Cu.P0700', 
	nu = setStoic(x, c('P0702.Cu.P0700', 'P0702', 'P0700.Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0700.Cu by P2581 F1']] = list(
	a = 'k.P0700.Cu.by.P2581.F1*P2581.Cu*P0700/V^2', 
	nu = setStoic(x, c('P2581.Cu', 'P0700', 'P2581.Cu.P0700'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P0700.Cu by P2581 R1']] = list(
	a = 'k.P0700.Cu.by.P2581.R1*P2581.Cu.P0700', 
	nu = setStoic(x, c('P2581.Cu.P0700', 'P2581.Cu', 'P0700'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0700.Cu by P2581 F2']] = list(
	a = 'k.P0700.Cu.by.P2581.F2*P2581.Cu.P0700', 
	nu = setStoic(x, c('P2581.Cu.P0700', 'P2581', 'P0700.Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0700.Cu dissociation']] = list(
	a = 'k.P0700.Cu.dissociation*P0700.Cu', 
	nu = setStoic(x, c('P0700.Cu', 'P0700', 'Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0700.Cu non-specific']] = list(
	a = 'k.P0700.Cu.non.specific*P0700*Cu/V^2', 
	nu = setStoic(x, c('P0700', 'Cu', 'P0700.Cu'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P0700 degradation']] = list(
	a = 'k.P0700.degradation*P0700', 
	nu = setStoic(x, 'P0700', -1),
	mods= NULL)

rxn[['P0700 translation']] = list(
	a = 'k.P0700.translation*M0700', 
	nu = setStoic(x, 'P0700', +1),
	mods= 'M0700')

rxn[['P0702.Cu bind']] = list(
	a = 'k.P0702.Cu.bind*P0702*Cu/V^2', 
	nu = setStoic(x, c('P0702', 'Cu', 'P0702.Cu'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P0702.Cu by P1179 F1']] = list(
	a = 'k.P0702.Cu.by.P1179.F1*P0702*P1179.Cu/V^2', 
	nu = setStoic(x, c('P0702', 'P1179.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P0702.Cu by P2581 F1']] = list(
  a = 'k.P0702.Cu.by.P2581.F1*P0702*P2581.Cu/V^2', 
	nu = setStoic(x, c('P0702', 'P2581.Cu', 'P0702.Cu.P2581'), c(-1, -1, +1)),
  mods= 'V')

rxn[['P0702.Cu by P2581 R1']] = list(
	a = 'k.P0702.Cu.by.P2581.R1*P0702.Cu.P2581',
	nu = setStoic(x, c('P0702.Cu.P2581', 'P0702', 'P2581.Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0702.Cu by Q F1']] = list(
	a = 'k.P0702.Cu.by.Q.F1*P0702*Q.Cu/V^2', 
	nu = setStoic(x, c('P0702', 'Q.Cu', 'P0702.Cu.Q'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P0702.Cu debind']] = list(
	a = 'k.P0702.Cu.debind*P0702.Cu', 
	nu = setStoic(x, c('P0702.Cu', 'P0702', 'Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P0702.Cu degradation']] = list(
	a = 'k.P0702.Cu.degradation*P0702.Cu', 
	nu = setStoic(x, c('P0702.Cu', 'Cu'), c(-1, +1)),
	mods= NULL)

rxn[['P0702 degradation']] = list(
	a = 'k.P0702.degradation*P0702', 
	nu = setStoic(x, 'P0702', -1),
	mods= NULL)

rxn[['P0702 translation']] = list(
	a = 'k.P0702.translation*M0702', 
	nu = setStoic(x, 'P0702', +1),
	mods= 'M0702')

rxn[['P2581.Cu bind']] = list(
	a = 'k.P2581.Cu.bind*P2581*Cu/V^2', 
	nu = setStoic(x, c('P2581', 'Cu', 'P2581.Cu'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P2581.Cu by P1179 F1']] = list(
	a = 'k.P2581.Cu.by.P1179.F1*P2581*P1179.Cu/V^2', 
	nu = setStoic(x, c('P2581', 'P1179.Cu', 'P2581.Cu.P1179'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P2581.Cu by P0702 F1']] = list(
	a = 'k.P2581.Cu.by.P0702.F1*P2581*P0702.Cu/V^2', 
	nu = setStoic(x, c('P2581', 'P0702.Cu', 'P0702.Cu.P2581'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P2581.Cu by P0702 R1']] = list(
	a = 'k.P2581.Cu.by.P0702.R1*P0702.Cu.P2581',
	nu = setStoic(x, c('P0702.Cu.P2581', 'P0702.Cu', 'P2581'), c(-1, +1, +1)),
	mods= NULL)
	
rxn[['P2581.Cu by Q F1']] = list(
	a = 'k.P2581.Cu.by.Q.F1*P2581*Q.Cu/V^2', 
	nu = setStoic(x, c('P2581', 'Q.Cu', 'P2581.Cu.Q'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P2581.Cu debind']] = list(
	a = 'k.P2581.Cu.debind*P2581.Cu', 
	nu = setStoic(x, c('P2581.Cu', 'P2581', 'Cu'), c(-1, +1, +1)),
	mods= NULL)
  
rxn[['P2581.Cu degradation']] = list(
	a = 'k.P2581.Cu.degradation*P2581.Cu', 
	nu = setStoic(x, c('P2581.Cu', 'Cu'), c(-1, +1)),
	mods= NULL)

rxn[['P2581 degradation']] = list(
	a = 'k.P2581.degradation*P2581', 
	nu = setStoic(x, 'P2581', -1),
	mods= NULL)
  
rxn[['P2581 translation']] = list(
	a = 'k.P2581.translation*M2581', 
	nu = setStoic(x, 'P2581', +1),
	mods= 'M2581')

rxn[['P1179.Cu by P0702 F1']] = list(
	a = 'k.P1179.Cu.by.P0702.F1*P1179*P0702.Cu/V^2', 
	nu = setStoic(x, c('P1179', 'P0702.Cu', 'P0702.Cu.P1179'), c(-1, -1, +1)),
	mods= 'V')

rxn[['P1179.Cu by P0702 R1']] = list(
	a = 'k.P1179.Cu.by.P0702.R1*P0702.Cu.P1179', 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P0702.Cu', 'P1179'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P1179.Cu by P0702 F2']] = list(
	a = 'k.P1179.Cu.by.P0702.F2*P0702.Cu.P1179', 
	nu = setStoic(x, c('P0702.Cu.P1179', 'P1179.Cu', 'P0702'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P1179.Cu by P2581 F1']] = list(
  a = 'k.P1179.Cu.by.P2581.F1*P1179*P2581.Cu/V^2', 
	nu = setStoic(x, c('P1179', 'P2581.Cu', 'P2581.Cu.P1179'), c(-1, -1, +1)),
  mods= 'V')

rxn[['P1179.Cu by P2581 R1']] = list(
	a = 'k.P1179.Cu.by.P2581.R1*P2581.Cu.P1179', 
	nu = setStoic(x, c('P2581.Cu.P1179', 'P2581.Cu', 'P1179'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P1179.Cu by P2581 F2']] = list(
  a = 'k.P1179.Cu.by.P2581.F2*P2581.Cu.P1179', 
	nu = setStoic(x, c('P2581.Cu.P1179', 'P1179.Cu', 'P2581'), c(-1, +1, +1)),
  mods= NULL)
  
rxn[['P1179.Cu dissociation']] = list(
	a = 'k.P1179.Cu.dissociation*P1179.Cu', 
	nu = setStoic(x, c('P1179.Cu', 'P1179', 'Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['P1179.Cu non-specific']] = list(
	a = 'k.P1179.Cu.non.specific*P1179*Cu/V^2', 
	nu = setStoic(x, c('P1179', 'Cu', 'P1179.Cu'), c(-1, -1, +1)),
	mods= 'V')

rxn[['Pgfp degradation']] = list(
	a = 'k.Pgfp.degradation*Pgfp', 
	nu = setStoic(x, 'Pgfp', -1),
	mods= NULL)

rxn[['Pgfp translation']] = list(
	a = 'k.Pgfp.translation*Mgfp', 
	nu = setStoic(x, 'Pgfp', +1),
	mods= 'Mgfp')

rxn[['Q.Cu by P0702 F1']] = list(
  a = 'k.Q.Cu.by.P0702.F1*P0702.Cu*Q/V^2', 
	nu = setStoic(x, c('P0702.Cu', 'Q', 'P0702.Cu.Q'), c(-1, -1, +1)),
  mods= 'V')

rxn[['Q.Cu by P0702 R1']] = list(
	a = 'k.Q.Cu.by.P0702.R1*P0702.Cu.Q', 
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702.Cu', 'Q'), c(-1, +1, +1)),
	mods= NULL)

rxn[['Q.Cu by P0702 F2']] = list(
  a = 'k.Q.Cu.by.P0702.F2*P0702.Cu.Q',
	nu = setStoic(x, c('P0702.Cu.Q', 'P0702', 'Q.Cu'), c(-1, +1, +1)),
  mods= NULL)

rxn[['Q.Cu by P2581 F1']] = list(
	a = 'k.Q.Cu.by.P2581.F1*P2581.Cu*Q/V^2', 
	nu = setStoic(x, c('P2581.Cu', 'Q', 'P2581.Cu.Q'), c(-1, -1, +1)),
	mods= 'V')

rxn[['Q.Cu by P2581 R1']] = list(
	a = 'k.Q.Cu.by.P2581.R1*P2581.Cu.Q', 
	nu = setStoic(x, c('P2581.Cu.Q', 'P2581.Cu', 'Q'), c(-1, +1, +1)),
	mods= NULL)

rxn[['Q.Cu by P2581 F2']] = list(
	a = 'k.Q.Cu.by.P2581.F2*P2581.Cu.Q', 
	nu = setStoic(x, c('P2581.Cu.Q', 'P2581', 'Q.Cu'), c(-1, +1, +1)),
	mods= NULL)

rxn[['Q.Cu non-specific']] = list(
	a = 'k.Q.Cu.non.specific*Q*Cu/V^2', 
	nu = setStoic(x, c('Q', 'Cu', 'Q.Cu'), c(-1, -1, +1)),
	mods= 'V')
	
rxn[['V growth']] = list(
	a = 'mu*V',
	nu = setStoic(x, 'V', +1),
	mods= NULL)

rxn[['Q growth']] = list(
#	a = paste('mu * ( ', paste(.T('Q'), collapse=' + '), ' )', sep='', collapse=''),
	a = paste('mu*', .T('Q'), sep='', collapse='+'),
	nu = setStoic(x, 'Q', +1),
	mods= .T('Q')[which(.T('Q')!='Q')])

rxn[['P1179 growth']] = list(
#  a = paste('mu*(', paste(.T('P1179'), collapse='+'), ')', sep='', collapse=''),
  a = paste('mu*', .T('P1179'), sep='', collapse='+'),
	nu = setStoic(x, 'P1179', +1),
	mods= .T('P1179')[which(.T('P1179')!='P1179')])

rxn[['D0700 growth']] = list(
#  a = paste('mu*(', paste(.T('D0700'), collapse='+'), ')', sep='', collapse=''),
  a = paste('mu*', .T('D0700'), sep='', collapse='+'),
	nu = setStoic(x, 'D0700', +1),
	mods= .T('D0700')[which(.T('D0700')!='D0700')])

rxn[['D0702 growth']] = list(
#  a = paste('mu*(', paste(.T('D0702'), collapse='+'), ')', sep='', collapse=''),
  a = paste('mu*', .T('D0702'), sep='', collapse='+'),
	nu = setStoic(x, 'D0702', +1),
	mods= .T('D0702')[which(.T('D0702')!='D0702')])

rxn[['D2581 growth']] = list(
#  a = paste('mu*(', paste(.T('D2581'), collapse='+'), ')', sep='', collapse=''),
  a = paste('mu*', .T('D2581'), sep='', collapse='+'),
	nu = setStoic(x, 'D2581', +1),
	mods= .T('D2581')[which(.T('D2581')!='D2581')])

rxn[['Dgfp growth']] = list(
#  a = paste('mu*(', paste(.T('Dgfp'), collapse='+'), ')', sep='', collapse=''),
  a = paste('mu*', .T('Dgfp'), sep='', collapse='+'),
	nu = setStoic(x, 'Dgfp', +1),
	mods= .T('Dgfp')[which(.T('Dgfp')!='Dgfp')])

rxn[['OE0702 growth']] = list(
  a = 'mu*OE0702',
	nu = setStoic(x, 'OE0702', +1),
  mods= NULL)

rxn[['OE2581 growth']] = list(
	a = 'mu*OE2581',
	nu = setStoic(x, 'OE2581', +1),
	mods= NULL)
