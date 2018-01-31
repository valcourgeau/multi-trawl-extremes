setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('gaver_lewis_simulation.R')

pw_lp1 <- function(x, alpha, beta){
  return((beta / (beta + x))^alpha)
}

pw_lp2_gl <- function(x1, x2, j, alpha, beta, rho){
  temp <- (beta + rho^j*x2) * beta
  temp <- (temp / (temp + x1))^alpha * pw_lp1(x1, alpha, beta)
  
  return(temp)
}

pw_lp2_warren <- function(x1, x2, j, alpha, beta, rho){
 temp <- 1 + (x1+x2)/beta + (1-rho^{j})*x1*x2/beta^2
 return(temp^{-alpha})
}

pw_lp1_dx <- function(x, alpha, beta){
  return(-alpha*(beta+x)*pw_lp1(x, alpha, beta))
}

pw_lp2_gl_dx2 <- function(x1, x2, j, alpha, beta, rho){
  temp <- alpha * x1 * (beta + x1) * rho^{j}
  temp <- temp / (beta * (beta + x2 * rho^{j})^2)
  temp <- temp * pw_lp2_gl(x1, x2, j, alpha, beta, rho)
  temp <- temp * pw_lp2_gl(x1, x2, j, 1.0, beta, rho)
  
  return(temp)
}

pw_lp2_gl_dx1 <- function(x1, x2, j, alpha, beta, rho){
  temp <- -alpha * x1 * (2*beta + x2*rho^{j} + 2*x1)
  temp <- temp / (beta * (beta + x2 * rho^{j}))
  temp <- temp * pw_lp2_gl(x1, x2, j, alpha, beta, rho)
  temp <- temp * pw_lp2_gl(x1, x2, j, 1.0, beta, rho)
  
  return(temp)
}

pw_lp2_gl_dx1dx2 <- function(x1, x2, j, alpha, beta, rho){
  temp <- (beta+x1)*(beta+x2*rho^{j})*(beta+x2*rho^{j}+x1)^2
  temp <- alpha * rho^{j} / temp
  temp <- temp * (beta^2+beta*x2*rho^{j}-2*alpha*beta*x1+beta*x1+(1-alpha)*x1*x2-2*alpha*x1^2)
  temp <- temp * pw_lp2_gl(x1, x2, j, alpha, beta, rho)
  
  return(temp)
}

pairwise_gl <- function(x1, x2, j, alpha, beta, rho, kappa, tol=1e-16){
  if(x1 < tol & x2 < tol){
    
    temp <- pw_lp2_gl(x1 = kappa,
                      x2 = kappa,
                      j = j-1,
                      alpha = alpha,
                      beta = beta,
                      rho = rho)
    return(1.0-2*pw_lp1(kappa, alpha = alpha, beta = beta)+temp)
  }else{
    if(x2 < tol){
      
      return(-pw_lp1_dx(x1+kappa, alpha, beta)+pw_lp2_gl_dx1(x1 = x1+kappa,
                                                       x2 = kappa,
                                                       j = j-1,
                                                       alpha = alpha,
                                                       beta = beta,
                                                       rho = rho))
    }else{
      if(x1 < tol){
       
        return(-pw_lp1_dx(x2+kappa, alpha, beta)+pw_lp2_gl_dx2(x1 = kappa,
                                                         x2 = x2+kappa,
                                                         j = j-1,
                                                         alpha = alpha,
                                                         beta = beta,
                                                         rho = rho))
      }else{
        # both not zero
        return(pw_lp2_gl_dx1dx2(x1 = x1+kappa,
                             x2 = x2+kappa,
                             j = j-1,
                             alpha = alpha,
                             beta = beta,
                             rho = rho))
      }
    }
  }
}

pl_gl_full <- function(data, params, delta){
  pl <- 0.0
  n <- length(data)
  
  # Decode params
  alpha <- params[1]
  beta <- params[2]
  rho <- params[3]
  kappa <- params[4]
  i <- 0
  for(first_index in 1:(n-1)){
    for(second_index in (first_index+1):min(n, first_index+delta)){
     #cat(data[first_index], data[second_index], "\n")
      temp <- pairwise_gl(x1 = data[first_index],
                          x2 = data[second_index],
                          j = second_index,
                          alpha = alpha,
                          beta = beta,
                          rho = rho,
                          kappa = kappa)
      #print(temp)
      if(is.na(temp) | is.nan(temp)){
        temp <- 0
      }
      if(temp > 1e-16){
        pl <- pl + log(temp)
      }else{
        i <- i + 1
        #cat(data[first_index], data[second_index], temp, "\n")
        pl <- pl - 1000
      }
    }
  }
  cat("Not accepted:", i/length(data), "\n")
  return(pl)
}

transform_to_alpha_beta <- function(shape, scale=NA){
  if(is.na(scale)){
    scale <- shape[2]
    shape <- shape[1]
  }
  alpha <- 1.0/shape
  beta <- scale/shape
  return(c("alpha"=alpha, "beta"=scale))
}

library(evir)
lgp_sim <- rlgprocess(alpha = 3.00,
                      beta = 1.0,
                      kappa = 1.7,
                      rho = 0.4,
                      timesteps = 10000,
                      n=1)

# 
# plot(density(rgpd(n = 10000, xi = 1/169, beta = (2168+27)/169)))
# lines(density(lgp_sim[lgp_sim > 0.0]))
params_init <- c(3.00, 1.0, 0.4, 1.7) * 0.8
k <- pl_gl_full(lgp_sim, params_init, 4)
-pl_gl_full(lgp_sim, c(3.00, 1.0, 0.4, 1.7), 4)

fn_to_optim <- function(params){
  print(params)
  -pl_gl_full(lgp_sim, params, 4) + 300*sum((params)^2)+150*sum((1/params)^2)
}
fn_to_optim(c(3.00, 1.0, 0.4, 1.7))
fn_to_optim(c(1.65, 0.96, 1.00, 1.51))

optim(fn = fn_to_optim, par=params_init, method="L-BFGS-B", control=list(trace=2, 
                                                                     parscale=c(0.5,0.5,0.5,0.5),
                                                                     pgtol = 1e-6),
      lower = rep(0.01, 4), upper = c(5, 5, 0.99, 5))

# 2.3751860 1.1335884 0.8419516 2.6835953


### Cross-validation

perform_cv <- function(data, params, kfold){
  n <- length(data)
  seq_stop <- seq(1, n+1, kfold)
  lambda_array <- c(10,100,400,500,600,1000,2000)
  for(index_without  in 1:(kfold-1)){
    fn_to_optim <- function(params){
      -pl_gl_full(data[-(seq_stop[index_without]):(seq_stop[index_without+1])], params, 4) + 300*sum((params)^2)+150*sum((1/params)^2)
    }
  }
  
}

