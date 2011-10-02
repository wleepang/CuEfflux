## Utility functions for CuEfflux Model

# ODE response function
# uses state change matrix and propensity vector to calculate a vector of 
# differential rates vector for ODE simulation
dxdt = function(time, x, parms, ...) {
  
  dxdt = rep(0, length(x))
  names(dxdt) = names(x)
 
  with(as.list(c(x, parms)), {
    eval_a = rep(0, length(a))
    names(eval_a) = names(a)
    parse_a = parse(text = a)
    for (.i in 1:length(parse_a)) {
      eval_a[.i] = eval(parse_a[.i])
    }
#     if (any(eval_a < 0)) {
#       warning("negative propensities coerced to zero")
#       eval_a[any(eval_a)] = 0
#     }
    
    dxdt = drop(sign(nu) %*% eval_a)
    return(list(dxdt))
  })
}


# returns names of all variant of a species for statistical operations
.T = function(species, data=NULL, fun=sum) {
	
	if (species == 'Cu') {
		col.ix = names(x)[(regexpr(species, names(x)) > 0) & (names(x) != 'Cu.out')]
	} else {
		col.ix = names(x)[(regexpr(species, names(x)) > 0)]
	}
	
	if (!is.null(data) && !is.null(fun)) {
		return(apply(data[, col.ix, drop=F], 1, fun))
	} else if (!is.null(data) && is.null(fun)) {
		return(data[, col.ix])
	} else {
		return(col.ix)
	}
}


###
# Defines the Cu import rate via the gradient between extra- and intra- cellular
# Cu abundance
k.Cu.Import = function(Cu.out, Cu.T, Cu.PartCoef = 1000, k.Cu.Import.BaseRate = 1e-5) {
	grad = Cu.out - Cu.T*Cu.PartCoef
	rate = 0
	if (grad >= 0) {
		rate = grad*k.Cu.Import.BaseRate
	}
	
	return(rate)
}

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

# creates the state-change matrix
set.nu = function(rxn) {
	nu = NULL
	for (r in rxn) nu = cbind(nu, r$nu)
	colnames(nu) = names(rxn)
	return(nu)
}

# creates the reaction propensity vector
set.a = function(rxn) {
	a = NULL
	for (r in rxn) a = c(a, r$a)
	names(a) = names(rxn)
	return(a)
}

# simpler console printing
printf = function(fmt, ...) {
	cat(sprintf(fmt, ...), sep='')
}

## performs scans of model
model.scan = function(ScanVar, ScanLevels, file.prefix='scan', t.final=36000, t.step=100) {
	printf('%s\n', file.prefix)
	printf('Scanning %d levels in %s:\n', length(ScanLevels), ScanVar)
	printf('%5s%15s%15s\n', '#', 'Value', 'Sim.Time (s)')
	
	iter = 1
	scan.result = NULL
	for (Level in ScanLevels) {
		x0[[ScanVar]] = Level
		
		printf('%5d%15e', iter, Level)
		
		tic()
		out = ode(x0, seq(0,t.final,by=t.step), dxdt, c(list(nu=nu, a=a), parms), method='daspk')
		t.done = toc(F)
		
		printf('%15.2f\n', t.done)
		
		# save this sim run for later processing if necessary
		save(out, file=sprintf('./RData/%s_iter-%05d.RData', file.prefix, iter))
		
		# collect endpoint values
		t.end = dim(out)[1]
		scan.result = rbind(scan.result, 
												c(id=iter, Level=Level, out[t.end, -1]))
		
		iter = iter + 1
	}
	printf('DONE\n')
	
	save(scan.result, nu, a, x0, file=sprintf('./RData/%s_result.RData', file.prefix))
	
	return(scan.result)

}

## perform 2d scan of model
model.scan.2d = function(ScanVars, ScanLevels, file.prefix='scan', t.final=36000, t.step=100) {
	if (length(ScanVars) < 2) stop('ScanVars must be of length 2')
	
	printf('%s\n', file.prefix)
	printf('Scanning %d levels in %s:\n', length(ScanLevels[[ScanVars[1]]]), ScanVars[1])
	printf('Scanning %d levels in %s:\n', length(ScanLevels[[ScanVars[2]]]), ScanVars[2])
	printf('%5s%15s%15s%15s\n', '#', 'Value-1', 'Value-2', 'Sim.Time (s)')
	
	iter = 1
	scan.result = NULL
	for (iLevel in ScanLevels[[ScanVars[1]]]) {
		x0[[ScanVars[1]]] = iLevel
		for (jLevel in ScanLevels[[ScanVars[2]]]) {
			x0[[ScanVars[2]]] = jLevel
			
			printf('%5d%15e%15e', iter, iLevel, jLevel)
			
			tic()
			out = ode(x0, seq(0,t.final,by=t.step), dxdt, list(nu=nu, a=a), method='daspk')
			t.done = toc(F)
			
			printf('%15.2f\n', t.done)
			
			# save this sim run for later processing if necessary
			save(out, file=sprintf('./RData/%s_iter-%05d.RData', file.prefix, iter))
			
			# collect endpoint values
			t.end = dim(out)[1]
			scan.result = rbind(scan.result, 
													c(id=iter, iLevel=iLevel, jLevel=jLevel, out[t.end, -1]))
			
			iter = iter + 1
		}
	}
	printf('DONE\n')
	
	save(scan.result, nu, a, x0, file=sprintf('./RData/%s_result.RData', file.prefix))
	
	return(scan.result)

}