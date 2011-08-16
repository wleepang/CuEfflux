# convert state change matrix and propensity vector to a vector of eval()-able
# differential vector for lsoda/ODE simulation
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
    if (any(eval_a < 0)) {
      warning("negative propensities coerced to zero")
      eval_a[any(eval_a)] = 0
    }
    
    # cell division is a special case that allows a negative propensity:
    # if (V >= V.max) {
    #  eval_a['V growth'] = -(V-V.min)
    # }
    
    dxdt = drop(sign(nu) %*% eval_a)
    return(list(dxdt))
  })
}

