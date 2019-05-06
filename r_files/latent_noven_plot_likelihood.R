setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")
source('pairwise_latent_trawl.R')


library('lattice')

t1 <- 1
t2 <- 2
rho <- 0.3
kappa <- 12
x1 <- 0.1
x2 <- 0.3

n <- 20
res_00 <- matrix(0, n, n)
res_10 <- matrix(0, n, n)
res_11 <- matrix(0, n, n)

i <- 1
j <- 1

alpha_seq <- seq(0.01, 25.00, length.out = n)
beta_seq <- seq(0.01, 70, length.out = n)

# alpha_seq <- seq(0.01, 8.00, length.out = n)
# beta_seq <- seq(10, 25, length.out = n)

for(alpha in alpha_seq){
  for(beta in beta_seq){
    res_00[i, j] <- pairwise_00_exp(t1=t1,t2=t2,alpha=alpha, beta = beta, rho = rho, kappa = kappa)
    res_10[i, j] <- pairwise_10_exp(t1=t1,x1=x1,t2=t2,alpha=alpha, beta = beta, rho = rho, kappa = kappa, transformation = F)
    res_11[i, j] <- pairwise_11_exp(t1=t1,x1=x1,t2=t2,x2=x2,alpha=alpha, beta = beta, rho = rho, kappa = kappa, transformation = F)
    
    j <- j + 1
    print(j)
  }
  i <- i + 1
  j <- 1
}

persp(alpha_seq, beta_seq, 
      res_00, phi = 10, theta = 210,
      xlab = "alpha", ylab = "beta",
      zlab = "likelihood",
      main = "log-likelihood 00",
      shade=1.5, col="red", border="white"
)


persp(alpha_seq, beta_seq, 
      res_10, phi = 30, theta = -30,
      xlab = "alpha", ylab = "beta",
      zlab = "likelihood",
      main = "log-likelihood 10",
      shade=1.5, col="red", border="white"
)


persp(alpha_seq, beta_seq, 
      res_11, phi = 40, theta = 190,
      xlab = "alpha", ylab = "beta",
      zlab = "likelihood",
      main = "log-likelihood 11",
      shade=1.5, col="red", border="white"
)
