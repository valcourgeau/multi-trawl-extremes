
diff.val <- function(value.p, value.m, epsilon){
  # centered differentiation
  
  return((value.p-value.m)/(2*epsilon))
}

sec.diff.val <- function(value.p, value.c, value.m, epsilon){
  temp <- value.p + value.m
  temp <- temp - 2.0*value.c
  temp <- temp / epsilon
  return(temp / epsilon)
}

sec.diff.mix.val <- function(value.c,
                             value.x.p, value.x.m,
                             value.y.p, value.y.m,
                             value.xy.p, value.xy.m,
                             epsilon){
  temp <- + value.xy.p + value.xy.m  + 2.0*value.c
  temp <- temp - value.x.m - value.y.m - value.x.p - value.y.p
  temp <- temp / epsilon
  return(temp/(2.0*epsilon))
}

grad.f <- function(f,
                   params,
                   epsilon=1e-6){
  params.fixed <- params
  d <- length(params)
  answer <- rep(0, d)
  
  for(index_par in 1:d){
    params.fixed[index_par] <- params[index_par] + epsilon
    value.p <- f(params.fixed)
    
    params.fixed[index_par] <- params[index_par] - epsilon
    value.m <- f(params.fixed)
    
    answer[index_par] <- diff.val(value.p = value.p,
                                  value.m = value.m,
                                  epsilon = epsilon)
    params.fixed <- params
  }
  
  return(answer)
}

evaluate.f <- function(f,
                       params,
                       epsilon){
  d <- length(params)
  eval.f <- matrix(0, d, d)
  params.fixed <- params
  
  for(main in 1:d){
    params.fixed[main] <- params[main] + epsilon
    for(second in main:d){
      if(second > main){
        params.fixed[second] <- params[second] + epsilon
        eval.f[main, second] <- f(params.fixed)
        params.fixed[second] <- params[second]
      }else{
        eval.f[main, second] <- f(params.fixed)
      }
      
    }
    
    params.fixed <- params
  }
  
  return(eval.f)
}

hessian.f <- function(f,
                      params,
                      epsilon=1e-6){
  if(epsilon < 0.0) stop("Epsilon must be positive.")
  
  d <- length(params)
  eval.f.p <- evaluate.f(f, params, epsilon = epsilon)
  eval.f.m <- evaluate.f(f, params, epsilon = -epsilon)
  hess.f <- matrix(0.0, d, d)
  f.value <- f(params)
  
  for(main in 1:d){
    for(second in main:d){
      if(main == second){
        hess.f[main, second] <- sec.diff.val(value.p = eval.f.p[main, second],
                                             value.m = eval.f.m[main, main],
                                             value.c = f.value,
                                             epsilon = epsilon)
      }else{
        hess.f[main, second] <- sec.diff.mix.val(value.x.p = eval.f.p[main, main],
                                                 value.x.m = eval.f.m[main, main],
                                                 value.y.p = eval.f.p[second, second],
                                                 value.y.m = eval.f.m[second, second],
                                                 value.xy.p = eval.f.p[main, second],
                                                 value.xy.m = eval.f.m[main, second],
                                                 value.c = f.value,
                                                 epsilon = epsilon)
      }
    }
  }
  
  hess.f <- hess.f + t(hess.f)
  diag(hess.f) <- 0.5*diag(hess.f)
  return(hess.f)
}

# test
## product of squares
f.t <- function(params){
  return(prod(params^2))
}

grad.f(f.t, params = c(0.5, 7), epsilon = 1e-6)
hessian.f(f.t, rep(1,2), epsilon = 1e-6)

## sum of squares
f.t <- function(params){
  return(sum(params^2))
}

grad.f(f.t, params = c(0.5, 7), epsilon = 1e-6)
hessian.f(f.t, rep(1,2), epsilon = 1e-6)

