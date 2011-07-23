GillespieSSA:::ssa <- function (x0 = stop("undefined 'x0'"), a = stop("undefined 'a'"), 
    nu = stop("undefined 'nu'"), parms = NULL, tf = stop("undefined 'tf'"), 
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
    args <- list(x0 = x0, a = a, nu = nu, parms = parms, tf = tf, 
        method = method, tau = tau, f = f, epsilon = epsilon, 
        nc = nc, hor = hor, dtf = dtf, nd = nd, ignoreNegativeState = ignoreNegativeState, 
        consoleInterval = consoleInterval, censusInterval = censusInterval, 
        verbose = verbose, simName = simName)
    out.rxn <- ssa.run(x0, a, nu, parms, tf, method, tau, f, 
        epsilon, nc, hor, dtf, nd, ignoreNegativeState, consoleInterval, 
        censusInterval, verbose, maxWallTime)
    out.summary <- ssa.terminate(args, out.rxn, tf, method, maxWallTime, 
        verbose)
    return(out.summary)
}

GillespieSSA:::ssa.run <- function (x0, a, nu, parms, tf, method, tau, f, epsilon, nc, 
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
    .t <- 0
    timeOfNextCensus <- .t + censusInterval
    timeForConsole <- .t + consoleInterval
    timeToTerminate <- FALSE
    timeSeries <- c(.t, .x)
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
        cat("t=", .t, " : ", sep = "")
        cat(.x, sep = ",")
        cat("\n")
        flush.console()
    }
    stepSize <- NULL
    currentRow <- 2
    suspendedTauLeapMethod <- FALSE
    nSuspendedTauLeaps <- 0
    while ((.t < tf) & (any(.x > 0)) & (all(.x >= 0)) & (any(eval_a > 
        0)) & (elapsedWallTime <= maxWallTime)) {
        doCalc <- TRUE
        if ((verbose) & (timeForConsole <= .t)) {
            cat("(", elapsedWallTime, "s) t=", .t, " : ", sep = "")
            cat(.x, sep = ",")
            cat("\n")
            flush.console()
            timeForConsole <- timeForConsole + consoleInterval
        }
        switch(method, D = {
            out <- ssa.d(eval_a, .nu)
            if (suspendedTauLeapMethod) {
                suspendedTauLeapMethod <- suspendedTauLeapMethod - 
                  1
                nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                if (!suspendedTauLeapMethod) method <- "OTL"
            }
        }, D.ETG = {
        		out <- ssa.d.etg(eval_a, a.s, .nu, .t, .x, partitionedSpecies)
        		
        }, ETL = {
            out <- ssa.etl(eval_a, .nu, .tau)
        }, BTL = {
            out <- ssa.btl(.x, eval_a, .nu, .f)
        }, OTL = {
            out <- ssa.otl(.x, eval_a, .nu, hor, .nc, .epsilon, 
                dtf, nd)
            suspendedTauLeapMethod <- out$suspendedTauLeapMethod
            if (suspendedTauLeapMethod) {
                method <- "D"
                doCalc <- FALSE
            }
        }, D.diag = {
            out <- ssa.d.diag(eval_a, .nu)
            if (suspendedTauLeapMethod) {
                suspendedTauLeapMethod <- suspendedTauLeapMethod - 
                  1
                nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                if (!suspendedTauLeapMethod) method <- "OTL.diag"
            }
        }, ETL.diag = {
            out <- ssa.etl.diag(eval_a, .nu, .tau)
        }, BTL.diag = {
            out <- ssa.btl.diag(.x, eval_a, .nu, .f)
        }, OTL.diag = {
            out <- ssa.otl.diag(.x, eval_a, .nu, hor, .nc, .epsilon, 
                dtf, nd)
            suspendedTauLeapMethod <- out$suspendedTauLeapMethod
            if (suspendedTauLeapMethod) {
                method <- "D.diag"
                doCalc <- FALSE
            }
        }, stop("unknown SSA method"))
        if (doCalc) {
            .t <- .t + out$tau
            .x <- .x + out$nu_j
            if ((any(.x < 0)) & (!ignoreNegativeState)) {
                cat("at least one population in 'x' is negative. Bailing to browser...\n")
                browser()
            }
            stepSize <- c(stepSize, out$tau)
            if (timeOfNextCensus <= .t) {
                timeSeries[currentRow, ] <- c(.t, .x)
                currentRow <- currentRow + 1
                timeOfNextCensus <- .t + censusInterval
                if (currentRow > dim(timeSeries)[1]) 
                  timeSeries <- rbind(timeSeries, matrix(nrow = 1000, 
                    ncol = (numCols)))
            }
            for (.i in seq(length(varNames))) assign(varNames[.i], 
                .x[[.i]])
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
        cat("t=", .t, " : ", sep = "")
        cat(.x, sep = ",")
        cat("\n")
        flush.console()
    }
    timeSeries <- timeSeries[!is.na(timeSeries[, 1]), ]
    timeSeries <- rbind(timeSeries, c(.t, .x))
    endWallTime <- format(Sys.time())
    return(list(timeSeries = timeSeries, eval_a = eval_a, elapsedWallTime = elapsedWallTime, 
        startWallTime = startWallTime, endWallTime = endWallTime, 
        stepSize = stepSize, nSuspendedTauLeaps = nSuspendedTauLeaps))
}

# new Exact Time-dependent Gillespie algorithm as described by
# Lu et.al Syst. Biol. (stevenage) Vol 1, No 1, 2004
GillespieSSA:::ssa.d.etg <- function (a   = stop("missing propensity vector (a)"), 
																			a.s = stop("missing time dependent propensity logical vector (a.s)"), 
																			nu  = stop("missing state-change matrix (nu)"), 
																	    t   = stop("missing time value"), 
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
  	tau = (1/c)*log(A.s / (A.s + c*log(1 - u1)))
  	
  } else if (A.q != 0 && A.s == 0) {
  	# only time independent reactions (standard Gillespie)
  	tau = -(1/A.q)*log(u1)
  }
  
  cellDivision = FALSE
  if ((t %% 1) + tau < 1) {
  	# select reaction channel
		j <- sample(seq(length(a)), size = 1, prob = a)
	  nu_j <- nu[, j]
  	# update time exactly
  } else {
    # update time to next cell division time
  	tau = 1 - (t %% 1)
  	
  	# non-reaction based particle change
  	# nu_j is used to divide cell populations in half
  	nu_j = rep(0, dim(nu)[1])
  	nu_j[x.p] = -(x[x.p] / 2)
  	
  	# reset volume
  	# partition "free" species
  	cellDivision = TRUE
  	
  }
  
  
  return(list(tau = tau, nu_j = nu_j, cellDivision = cellDivision))
}