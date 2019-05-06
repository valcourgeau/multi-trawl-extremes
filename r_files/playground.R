alpha <- 9
beta <- 4
rho <- 0.2
time <- 3
kappa <- 0.3

(1+kappa/beta)^{-alpha}

k <- 5
s <- 2

moment_gamma <- function(k, alpha, beta){
  prod(1:k + alpha-1)/beta^k
}

link_matrix <- function(k, s, alpha, beta){
  link_ks <- moment_gamma(k = k+s, alpha = alpha, beta = beta)
  link_k <- moment_gamma(k = k, alpha = alpha, beta = beta)
  link_s <- moment_gamma(k = s, alpha = alpha, beta = beta)
  coeff <- factorial(k) * factorial(s)
  
  return((link_ks - link_k * link_s)/coeff)
}

link_matrix(k, s, alpha, beta)
link_matrix(s, k, alpha, beta)

cross_integrated <- function(k, alpha, beta, kappa){
  if(alpha-1 <= k | k < 0) stop('Wrong: need alpha-1 > k >= 0')
  result <- beta/(alpha-k-1)*(1+kappa/beta)^{1-(alpha-k)}
  last <- result
  if(k >= 1){
    for(index in k:1){
      last <- beta/(alpha-index)*(kappa^{index}*(1+kappa/beta)^{-(alpha-index)} + index * last)
      #result <- cbind(result, last)
    }
  }
  return(last)
}

cross_integrated(k, alpha, beta, kappa)

cross_integrated(0, alpha, beta, kappa)
beta/(alpha-1)*(1+kappa/beta)^{-alpha + 1}

cross_integrated(1, alpha, beta, kappa)
beta/(alpha-1)*(kappa*(1+kappa/beta)^{-alpha + 1} + 1 * cross_integrated(0, alpha-1, beta, kappa))

covariance_exceed <- function(time, alpha, beta, kappa){
  alpha_exc <- alpha * (1 - exp(-rho*time))
  alpha_inter <- alpha * exp(-rho*time)
  k_max <- 3
  
  A_values <- vapply(1:k_max, function(k){return((-1)^(k)*cross_integrated(k, alpha_exc, beta, kappa))}, 1.0) #integrals u^k
  link_values <- vapply(1:k_max, function(k){return(vapply(1:k_max, function(s){link_matrix(k, s, alpha_inter, beta)}, 1.0))}, rep(0, k_max))
  #return(link_values)
  return(A_values)
  return((A_values) %*% link_values %*% A_values)
} 

covariance_exceed(time, alpha, beta, kappa)/((beta/alpha)^2/((1-1/alpha)^2*(1-2/alpha)))

alpha/beta^2*exp(-rho*time)                      

