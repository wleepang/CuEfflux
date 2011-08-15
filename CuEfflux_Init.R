## initialize Cu Efflux model

# Define parameters
#parms <- c(c1=1.0, c2=0.002, c3=0.5, c4=0.04)

kGenomeCopy = 15
isKO = FALSE
isOE = FALSE
mu = 1/(7*60*60) # cell growth rate 1/s

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
				Q.Cu            = 0, 
				V               = 1000)
