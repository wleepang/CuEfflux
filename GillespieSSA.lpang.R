ssa <- function (x0 = stop("undefined 'x0'"), a = stop("undefined 'a'"), 
    nu = stop("undefined 'nu'"), parms = NULL, ti = 0, tf = stop("undefined 'tf'"), 
    method = "D", simName = "", tau = 0.3, f = 10, epsilon = 0.03, 
    nc = 10, hor = NaN, dtf = 10, nd = 100, ignoreNegativeState = TRUE, 
    consoleInterval = 0, censusInterval = 0, verbose = FALSE, 
    maxWallTime = Inf) 
{
    ssa.check.args(x0, a, nu, tf, method, tau, f, epsilon, nc, 
        hor, dtf, nd, ignoreNegativeState, consoleInterval, censusInterval, 
        verbose)
    
    method = toupper(method)
    
    ssa.check.method(x0, a, nu, method, tau, f)
    if ((length(a)/dim(nu)[2] > 1) && (length(x0)/dim(nu)[1]) > 
        1) {
        if (method == "D") 
            method <- "D.diag"
        if (method == "ETL") 
            method <- "ETL.diag"
        if (method == "BTL") 
            method <- "BTL.diag"
        if (method == "OTL") 
            method <- "OTL.diag"
    }
    args <- list(x0 = x0, a = a, nu = nu, parms = parms, ti = ti, tf = tf, 
        method = method, tau = tau, f = f, epsilon = epsilon, 
        nc = nc, hor = hor, dtf = dtf, nd = nd, ignoreNegativeState = ignoreNegativeState, 
        consoleInterval = consoleInterval, censusInterval = censusInterval, 
        verbose = verbose, simName = simName)
    out.rxn <- ssa.run(x0, a, nu, parms, ti, tf, method, tau, f, 
        epsilon, nc, hor, dtf, nd, ignoreNegativeState, consoleInterval, 
        censusInterval, verbose, maxWallTime)
    out.summary <- ssa.terminate(args, out.rxn, tf, method, maxWallTime, 
        verbose)
    return(out.summary)
}

ssa.run <- function (x0, a, nu, parms, ti, tf, method, tau, f, epsilon, nc, 
    hor, dtf, nd, ignoreNegativeState, consoleInterval, censusInterval, 
    verbose, maxWallTime) 
{
    .x0 <- x0
    .a <- a
    .nu <- nu
    .tau <- tau
    .f <- f
    .epsilon <- epsilon
    .nc <- nc
    .x <- .x0
    varNames <- names(.x)
    for (.i in seq(length(varNames))) assign(varNames[.i], .x[[.i]])
    if (!is.null(parms)) {
        parmsNames <- names(parms)
        for (.i in seq(length(parmsNames))) assign(parmsNames[.i], 
            parms[[.i]])
    }
    
    # define .t as simulation time
    # define .T as model time
    .t <- 0
    .T <- ti 			# model initial time
    TF <- tf      # model final time
    tf <- tf - ti
    
    # difference between model time and simulation time
    .dTt <- .T - .t
    
    timeOfNextCensus <- .t + censusInterval
    timeForConsole <- .t + consoleInterval
    timeToTerminate <- FALSE
    timeSeries <- c(.T, .x)
    numCols <- length(timeSeries)
    timeSeries <- rbind(timeSeries, matrix(nrow = 1000, ncol = (numCols)))
    .M <- length(.a)
    eval_a <- rep(0, .M)
    parse_a <- parse(text = .a)
    for (.i in seq(length(parse_a))) eval_a[.i] <- eval(parse_a[.i])
    if (any(eval_a < 0)) 
        stop("negative propensity function")
    if (method == "OTL") {
        if (any(is.na(hor))) 
            hor <- rep(2, length(.x0))
        else if (length(hor) != length(.x0)) 
            stop("length of hor vector is different from length of 'x0'")
        else if (any(hor != 1 & hor != 2 & hor != 22)) 
            stop("wrong value(s) in hor vector (can only be 1, 2, or 22)")
    }
    procTimeStart <- proc.time()
    elapsedWallTime <- 0
    startWallTime <- format(Sys.time())
    if (verbose) {
        cat("Running ", method, " method with console output every ", 
            consoleInterval, " time step\n", sep = "")
        cat("Start wall time: ", startWallTime, "...\n", sep = "")
        flush.console()
    }
    if ((verbose) & (consoleInterval > 0)) {
        cat("t (sim | model) = ", .t, " | ", .T, " : ", sep = "")
        cat(.x, sep = ",")
        cat("\n")
        flush.console()
    }
    stepSize <- NULL
    currentRow <- 2
    suspendedTauLeapMethod <- FALSE
    nSuspendedTauLeaps <- 0
    while ((.t < tf) & (any(.x > 0)) & (all(.x >= 0)) & (any(eval_a > 0)) & (elapsedWallTime <= maxWallTime)) {
        doCalc <- TRUE
        if ((verbose) & (timeForConsole <= .t)) {
            cat("(", elapsedWallTime, "s) t (sim | model) = ", .t, " | ", .T," : ", sep = "")
            cat(.x, sep = ",")
            cat("\n")
            flush.console()
            timeForConsole <- timeForConsole + consoleInterval
        }
        switch(method, D = {
            out <- ssa.d(eval_a, .nu)
            if (suspendedTauLeapMethod) {
                suspendedTauLeapMethod <- suspendedTauLeapMethod - 1
                nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                if (!suspendedTauLeapMethod) method <- "OTL"
            }
        }, D.ETG = {
        		out <- ssa.d.etg(eval_a, a.s, .nu, .T, .x, partitionedSpecies)
        }, ETL = {
            out <- ssa.etl(eval_a, .nu, .tau)
        }, BTL = {
            out <- ssa.btl(.x, eval_a, .nu, .f)
        }, OTL = {
            out <- ssa.otl(.x, eval_a, .nu, hor, .nc, .epsilon, dtf, nd)
            suspendedTauLeapMethod <- out$suspendedTauLeapMethod
            if (suspendedTauLeapMethod) {
                method <- "D"
                doCalc <- FALSE
            }
        }, D.diag = {
            out <- ssa.d.diag(eval_a, .nu)
            if (suspendedTauLeapMethod) {
                suspendedTauLeapMethod <- suspendedTauLeapMethod - 1
                nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                if (!suspendedTauLeapMethod) method <- "OTL.diag"
            }
        }, ETL.diag = {
            out <- ssa.etl.diag(eval_a, .nu, .tau)
        }, BTL.diag = {
            out <- ssa.btl.diag(.x, eval_a, .nu, .f)
        }, OTL.diag = {
            out <- ssa.otl.diag(.x, eval_a, .nu, hor, .nc, .epsilon, dtf, nd)
            suspendedTauLeapMethod <- out$suspendedTauLeapMethod
            if (suspendedTauLeapMethod) {
                method <- "D.diag"
                doCalc <- FALSE
            }
        }, stop("unknown SSA method"))
        if (doCalc) {
            .t <- .t + out$tau
            .T <- .T + out$tau
            .x <- .x + out$nu_j
            if ((any(.x < 0)) & (!ignoreNegativeState)) {
                cat("at least one population in 'x' is negative. Bailing to browser...\n")
                browser()
            }
            stepSize <- c(stepSize, out$tau)
            if (timeOfNextCensus <= .t) {
                timeSeries[currentRow, ] <- c(.T, .x)
                currentRow <- currentRow + 1
                timeOfNextCensus <- .t + censusInterval
                if (currentRow > dim(timeSeries)[1]) 
                  timeSeries <- rbind(timeSeries, matrix(nrow = 1000, 
                    ncol = (numCols)))
            }
            
            for (.i in seq(length(varNames))) assign(varNames[.i], .x[[.i]])
            for (.i in seq(length(parse_a))) eval_a[.i] <- eval(parse_a[.i])
            
            eval_a[is.na(eval_a)] <- 0
            
            if (any(eval_a < 0)) 
                warning("negative propensity function - coersing to zero\n")
            
            eval_a[eval_a < 0] <- 0
        }
        procTimeEnd <- proc.time()
        elapsedWallTime <- procTimeEnd[3] - procTimeStart[3]
    }
    if (verbose) {
        cat("t (sim | model) = ", .t, " | ", .T," : ", sep = "")
        cat(.x, sep = ",")
        cat("\n")
        flush.console()
    }
    timeSeries <- timeSeries[!is.na(timeSeries[, 1]), ]
    timeSeries <- rbind(timeSeries, c(.T, .x))
    endWallTime <- format(Sys.time())
    return(list(timeSeries = timeSeries, eval_a = eval_a, elapsedWallTime = elapsedWallTime, 
        startWallTime = startWallTime, endWallTime = endWallTime, 
        stepSize = stepSize, nSuspendedTauLeaps = nSuspendedTauLeaps))
}

# new Exact Time-dependent Gillespie algorithm as described by
# Lu et.al Syst. Biol. (stevenage) Vol 1, No 1, 2004
ssa.d.etg <- function (a   = stop("missing propensity vector (a)"), 
																			a.s = stop("missing time dependent propensity logical vector (a.s)"), 
																			nu  = stop("missing state-change matrix (nu)"), 
																			.T  = stop("missing model time value"), 
																			x   = stop("missing populations vector"),
																			x.p = stop("missing partitioned species vector")) 
{
  A.s = sum(a[a.s])
  A.q = sum(a[!a.s])
  c = log(2)
  
  u1 = runif(1)
  
  if (A.q != 0 && A.s != 0) {
  	# time dependent and time independent reactions
  	alpha = A.s/A.q
  	beta = c*log(u1)/A.q
  	tau = (1/c)*lambertW(alpha*exp(alpha + beta))-alpha-beta
  	
  } else if (A.q == 0 && A.s != 0) {
  	# only time dependent reactions
    if (u1 > 1-exp(-A.s/c)) {
    	# advance to next division
    	tau = 1- (.T %% 1)
    } else {
    	tau = (1/c)*log(A.s / (A.s + c*log(1 - u1)))
    }
  	
  } else if (A.q != 0 && A.s == 0) {
  	# only time independent reactions (standard Gillespie)
  	tau = -(1/A.q)*log(u1)
  }
  
  # cell division is triggered by model time not simulation time
  # thus, to randomize starting volume ti (initial time) must be reflective of
  # initial starting volume according to V(T) = V0*exp(c*T)
  # (note: model time T should be normalized by the cell division time s.t.:
  #   t.d = log(2)/mu
  #   T = .t / t.d
  # )
  cellDivision = FALSE
  if (((.T %% 1) + tau) < 1) {
  	# select reaction channel
		j <- sample(seq(length(a)), size = 1, prob = a)
	  nu_j <- nu[, j]
  	# update time exactly
  } else {
    # update time to next cell division time
  	tau = 1 - (.T %% 1)
  	
  	# non-reaction based particle change
  	# nu_j is used to divide cell populations in half
  	# note: V can be a partitionable species
  	nu_j = rep(0, dim(nu)[1])
  	nu_j[x.p] = -floor(x[x.p] / 2)
  	
  	# reset volume
  	# partition "free" species
  	cellDivision = TRUE
  	
  }
  
  print(c(.T=.T, A.s=A.s, A.q=A.q, u1=u1, tau=tau, cellDivision=cellDivision))
  
  
  return(list(tau = tau, nu_j = nu_j, cellDivision = cellDivision))
}

# Lambert W function from:
# https://stat.ethz.ch/pipermail/r-help/2003-November/042793.html
lambertW = function(z, b=0, maxiter=10, eps=.Machine$double.eps, min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
        ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}