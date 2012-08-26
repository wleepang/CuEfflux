## Initial state vector
x0 <- c(Cu							= 0, 
				Cu.out					= 500000,
				D0700           = kGenomeCopy,
				Dgfp            = kGenomeCopy,
				M0700           = 0,
				M0702           = 0,
				M2581           = 0,
				Mgfp            = 0,
				P0700           = 0,
				P0700.Cu        = 0,
				P0702.Cu        = 0,
				P0702.Cu.P0700  = 0,
				P0702.Cu.P1179  = 0,
				P0702.Cu.Q      = 0,
				P2581.Cu        = 0,
				P2581.Cu.P0700  = 0,
				P2581.Cu.P1179  = 0,
				P2581.Cu.Q      = 0,
				P0702.P0702     = 0,
				P0702.P2581			= 0,
				P2581.P2581     = 0,
				P0702.Cu.P0702  = 0,
				P0702.Cu.P2581  = 0,
				P2581.Cu.P2581  = 0,
				P1179           = kGenomeCopy*6, # at least two regulators per controlled gene (excluding gfp)
				P1179.Cu        = 0,
				P1179.Cu.D0700  = 0,
				P1179.Cu.D0702  = 0,
				P1179.Cu.D2581  = 0,
				P1179.Cu.Dgfp   = 0,
				Pgfp            = 0,
				Q               = 1000,
				Q.Cu            = 0, 
				X               = 0,
				X.Cu            = 0,
				P2581.Cu.X      = 0,
				P0700.Cu.X      = 0,
				V               = 1)

if (modelVersion == 0) {
  x0 = c(x0,
  			D0702           = kGenomeCopy*2*(!isKO0702),
				D2581           = 0,
  			OE0702          = kGenomeCopy*isOE0702,
				OE2581          = 0,
  			P0702           = 0*(!isKO0702),
  			P2581           = 0)
} else if (modelVersion %in% c(1.1, 1.2, 1.3)) {
  x0 = c(x0,
    		D0702           = kGenomeCopy*(!isKO0702),
				D2581           = kGenomeCopy*(!isKO2581),
  			OE0702          = kGenomeCopy*isOE0702,
				OE2581          = kGenomeCopy*isOE2581,
  			P0702           = 350*(!isKO0702),
  			P2581           = 350*(!isKO2581))
}

# Initialize parameter set based on model name
if (modelVersion %in% c(0)) {
	# disable all reactions involving 2581
	parms[regexpr('2581', names(parms)) > 0] = 0
}
 
if (modelVersion %in% c(1.1)) {
	# weak interaction between P0702 and Q
	parms['k.Q.Cu.by.P0702.F2'] = 0.1
	
	# weak interaction between P2581 and P1179
	parms['k.P1179.Cu.by.P2581.F2'] = 0.01
	
	# weak interaction between P2581 and P0700
	parms['k.P0700.Cu.by.P2581.F2'] = 0.1
	
	# disable interaction between P2581 and P0702
	parms[regexpr('(.*P0702.*P2581)|(.*P2581.*P0702)', names(parms)) > 0] = 0
	
}
 
if (modelVersion %in% c(1.2)) {
	# disable direct Cu binding by P0702
	parms['k.P0702.Cu.bind'] = 0
	
	# disable interactions between P0702 and P1179, Q
	parms[regexpr('(.*P0702.*P1179)|(.*P1179.*P0702)', names(parms)) > 0] = 0
  parms[regexpr('(.*P0702.*Q)|(.*Q.*P0702)', names(parms)) > 0] = 0
	
	# disable interactions between P2581 and P0700
	parms[regexpr('(.*P0700.*P2581)|(.*P2581.*P0700)', names(parms)) > 0] = 0
	
}

if (modelVersion %in% c(1.3)) {
	# basal binding of Cu to P0700
	parms['k.P0702.Cu.bind'] = 0.0002
	
	# basal P2581 transfer of Cu to P0700
	parms['k.P0700.Cu.by.P2581.F1'] = 0.001
	
  # disable interactions between P0702 and P1179, Q
	parms[regexpr('(.*P0702.*P1179)|(.*P1179.*P0702)', names(parms)) > 0] = 0
  parms[regexpr('(.*P0702.*Q)|(.*Q.*P0702)', names(parms)) > 0] = 0
		
	# P2581.Cu.P0702 equillibrium
	parms['k.P2581.Cu.by.P0702.R1'] = 1
	parms['k.P0702.Cu.by.P2581.R1'] = 1
	
	# increase P1179.Cu by P2581 activation rate
	parms['k.P1179.Cu.by.P2581.R1'] = 0.01
	parms['k.P1179.Cu.by.P2581.F2'] = 0.05
	
	# make P2581 expression the same as P0702
	parms['k.M2581.transcription'] = parms['k.M0702.transcription']
	parms['k.M2581.degradation'] = parms['k.M0702.degradation']
	parms['k.OE2581.transcription'] = parms['k.OE0702.transcription']
	parms['k.P2581.translation'] = parms['k.P0702.translation']
	
	# chaperone degradation, needs to be symmetric, only apo chaperones
	parms['k.P0702.degradation'] = 0.0002
	parms['k.P2581.degradation'] = 0.0005
}