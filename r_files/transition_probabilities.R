### 
# Implement transition probabilities for univariate trawl-latent model
###

transi.p.geq.given.geq <- function(v, u, alpha, beta, kappa, B1, B_inter, A){
  # compute P(X_{t+h} > v | X_{t} > u]
  if(abs(A-B1-B_inter)>1e-3)stop("Wrong values for A, B1, B_inter.")

  b1 <- -alpha*B1/A
  b_inter <- -alpha*B_inter/A
  b3 <- -alpha*B3/A
  res <- (1+(2*kappa+u+v)/beta)^b_inter
  res <- res*(1+(kappa+u))^{b3+alpha}
  res <- res*(1+(kappa+v))^{b1}
  
  return(res)
}

transi.p.zero.given.geq <- function(u, alpha, beta, kappa, B1, B_inter, A){
  # compute P(X_{t+h} = 0 | X_{t} > u]
  return(1-transi.p.geq.given.geq(v = 0, u, alpha,beta,rho,kappa,B1,B_inter,B3,A))
}

# transi.p.geq.given.geq(v = 0, u = 1, 
#                        alpha = 8, beta = 4.37,
#                        kappa = 1.96,
#                        B1 = 1-exp(-0.26*1), B_inter = exp(-0.26*1), A = 1)

transi.density.geq.given.geq <- function(v, u, alpha, beta, kappa, B1, B_inter, A){
  # compute f(X_{t+h} > v | X_{t} > u]
  if(abs(A-B1-B_inter)>1e-3)stop("Wrong values for A, B1, B_inter.")
  
  b1 <- -alpha*B1/A
  b_inter <- -alpha*B_inter/A
  b3 <- -alpha*B3/A
  res <- (1+(2*kappa+u+v)/beta)^{b_inter-1}
  res <- res*(1+(kappa+u))^{b3+alpha}
  res <- res*(1+(kappa+v))^{b1-1}
  res <- res/beta^2
  res <- res * (b1*(1+(kappa+v)/beta) + b_inter*(1+(kappa+u+v)/beta))
  
  return(-res)
}