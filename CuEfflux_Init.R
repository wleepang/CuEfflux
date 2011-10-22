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
				P0702.Cu.P2581  = 0,
				P1179           = kGenomeCopy*6, # at least two regulators per controlled gene (excluding gfp)
				P1179.Cu        = 0,
				P1179.Cu.D0700  = 0,
				P1179.Cu.D0702  = 0,
				P1179.Cu.D2581  = 0,
				P1179.Cu.Dgfp   = 0,
				Pgfp            = 0,
				Q               = 1000,
				Q.Cu            = 0, 
				V               = 1)

if (modelName == 'onechap') {
  x0 = c(x0,
  			D0702           = kGenomeCopy*2*(!isKO0702),
				D2581           = 0,
  			OE0702          = kGenomeCopy*isOE0702,
				OE2581          = 0,
  			P0702           = 0*(!isKO0702),
  			P2581           = 0)
} else if (modelName %in% c('twochaps.nb', 'twochaps.lp')) {
  x0 = c(x0,
    		D0702           = kGenomeCopy*(!isKO0702),
				D2581           = kGenomeCopy*(!isKO2581),
  			OE0702          = kGenomeCopy*isOE0702,
				OE2581          = kGenomeCopy*isOE2581,
  			P0702           = 0*(!isKO0702),
  			P2581           = 0*(!isKO2581))
}