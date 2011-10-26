## global parameters
modelName = 'full' # options: full, onechap, twochaps.nb, twochaps.lp
kGenomeCopy = 25
isKO0702 = FALSE
isKO2581 = FALSE
isOE0702 = FALSE
isOE2581 = FALSE
mu = log(2)/(7*60*60) # cell growth rate 1/s

## reaction parameters
parms = c(
    k.Cu.export = 0.1,
    Cu.PartCoef = 1000,
    k.Cu.Import.BaseRate = 1e-5,
    k.D0700.P1179.Cu.binding = 0.01,
    k.D0700.P1179.Cu.dissociation = 0.1,
    k.D0702.P1179.Cu.binding = 0.01,
    k.D0702.P1179.Cu.dissociation = 0.1,
    k.D2581.P1179.Cu.binding = 0.01,      # twochap.xx
    k.D2581.P1179.Cu.dissociation = 0.1,  # twochap.xx
    k.Dgfp.P1179.Cu.binding = 0.01,
    k.Dgfp.P1179.Cu.dissociation = 0.1,
    k.M0700.degradation = 0.00180197,
    k.M0700.transcription = 0.207918,
    k.M0702.degradation = 0.00136725,
    k.M0702.transcription = 0.646571,
    k.M2581.degradation = 0.003262,     # twochap.xx
    k.M2581.transcription = 0.424167,   # twochap.xx
    k.Mgfp.degradation = 0.001,
    k.Mgfp.transcription = 0.138889,
    k.OE0702.transcription = 0.646571/100, # in order to get OE to be like expt. divide by factor of 100
    k.OE2581.transcription = 0.424167/100, # twochap.xx, in order to get OE to be like expt. divide by factor of 100
    k.P0700.Cu.by.P0702.F1 = 0.01,
    k.P0700.Cu.by.P0702.R1 = 0.1,
    k.P0700.Cu.by.P0702.F2 = 1.0,
    k.P0700.Cu.by.P2581.F1 = 0.01,       # twochap.nb
    k.P0700.Cu.by.P2581.R1 = 0.1,        # twochap.nb
    k.P0700.Cu.by.P2581.F2 = 1.0,        # nb = 0.1, twochap.nb, 10-fold less efficient than 0702
    k.P0700.Cu.dissociation = 0.0001,
    k.P0700.Cu.non.specific = 0.001,
    k.P0700.degradation = 0,
    k.P0700.translation = 0.002331,
    k.P0702.Cu.bind = 0.01,             # onechap, twochap.nb
    k.P0702.Cu.by.P1179.F1 = 0.01,      # onechap, twochap.nb
    k.P0702.Cu.by.Q.F1 = 0.01,          # onechap, twochap.nb
    k.P0702.Cu.by.P2581.F1 = 0.01,      # twochap.lp
    k.P0702.Cu.debind = 0.001,
    k.P0702.Cu.degradation = 0,
    k.P0702.degradation = 0,
    k.P0702.translation = 0.0289855,
    k.P2581.Cu.bind = 0.01,             # twochap.xx
    k.P2581.Cu.by.P1179.F1 = 0.01,      # twochap.xx
    k.P2581.Cu.by.P0702.F1 = 0.01,      # twochap.xx
    k.P2581.Cu.by.Q.F1 = 0.01,          # twochap.xx
    k.P2581.Cu.debind = 0.001,          # twochap.xx
    k.P2581.Cu.degradation = 0,
    k.P2581.degradation = 0, 
 		k.P2581.translation = 0.030769231,  # twochap.xx
    k.P0702.Cu.by.P2581.R1 = 0.1,       # twochap.xx
    k.P2581.Cu.by.P0702.R1 = 0.1,       # twochap.xx
    k.P1179.Cu.by.P0702.F1 = 0.01,       # onechap, twochap.nb
    k.P1179.Cu.by.P0702.R1 = 0.1,        # onechap, twochap.nb
    k.P1179.Cu.by.P0702.F2 = 0.1,        # onechap, twochap.nb
    k.P1179.Cu.by.P2581.F1 = 0.01,       # twochap.xx
    k.P1179.Cu.by.P2581.R1 = 0.1,        # twochap.xx
    k.P1179.Cu.by.P2581.F2 = 0.1,    # twochap: nb = 0.01, lp = 0.1
    k.P1179.Cu.dissociation = 0.0001,
    k.P1179.Cu.non.specific = 0.001,
    k.Pgfp.degradation = 0,
    k.Pgfp.translation = 0.008333,
    k.Q.Cu.by.P0702.F1 = 0.01,           # onechap, twochap.nb
    k.Q.Cu.by.P0702.R1 = 0.1,            # onechap, twochap.nb
    k.Q.Cu.by.P0702.F2 = 1.0,            # nb = 0.1, onechap, twochap.nb
    k.Q.Cu.by.P2581.F1 = 0.01,           # twochap.xx
    k.Q.Cu.by.P2581.R1 = 0.1,            # twochap.xx
    k.Q.Cu.by.P2581.F2 = 1.0,            # twochap.xx
    k.Q.Cu.non.specific = 0.001
  )
  
if (modelName %in% c('onechap')) {
	# disable all reactions involving 2581
	parms[regexpr('2581', names(parms)) > 0] = 0
}
 
if (modelName %in% c('twochaps.nb')) {
	# weak interaction between P0702 and Q
	parms['k.Q.Cu.by.P0702.F2'] = 0.1
	
	# weak interaction between P2581 and P1179
	parms['k.P1179.Cu.by.P2581.F2'] = 0.01
	
	# weak interaction between P2581 and P0700
	parms['k.P0700.Cu.by.P2581.F2'] = 0.1
	
	# disable interaction between P2581 and P0702
	parms[regexpr('(.*P0702.*P2581)|(.*P2581.*P0702)', names(parms)) > 0] = 0
}
 
if (modelName %in% c('twochaps.lp')) {
	# disable direct Cu binding by P0702
	parms['k.P0702.Cu.bind'] = 0
	
	# disable interactions between P0702 and P1179, Q
	parms[regexpr('(.*P0702.*P1179)|(.*P1179.*P0702)|(.*P0702.*Q)|(.*Q.*P0702)', names(parms)) > 0] = 0
	
	# disable interactions between P2581 and P0700
	parms[regexpr('(.*P0700.*P2581)|(.*P2581.*P0700)', names(parms)) > 0] = 0
}
 