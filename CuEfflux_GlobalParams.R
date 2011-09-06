## global parameters
modelName = 'onechap' # options: onechap, twochaps.nb, twochaps.lp
kGenomeCopy = 25
isKO0702 = FALSE
isKO2581 = FALSE
isOE0702 = FALSE
isOE2581 = FALSE
isKO = isKO0702 || isKO2581
isOE = isOE0702 || isOE2581
mu = log(2)/(7*60*60) # cell growth rate 1/s