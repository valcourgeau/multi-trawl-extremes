setwd("~/GitHub/multi-trawl-extremes/r_files/")
source('latent_novel_simulation.R')
source('pairwise_latent_trawl.R')
source("imp_latent_noven_final.R")
source("latent_novel_simulation.R")
source("imp_warren.R")

n_sims <- 10000
n_timesteps <- 30

times <- 1:(n_timesteps+2*30)
alpha <- 4
beta <- 1
rho <- 0.05
kappa <- 0.75

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
exceeds <- matrix(0, nrow = length(times)-2*30, ncol = n_sims)
for(i in 1:n_sims){
  exceeds[,i] <- rltrawl(alpha = alpha,
                          beta = beta,
                          kappa = 0.0,
                          trawl_fs = trawl_1,
                          trawl_fs_prim = trawl_1_prim,
                          times = times,
                          n = 1,
                          transformation = F)[1:n_timesteps]
  print(i)
}

write.table(x = t(exceeds), file = "acf_estimation_sim_10000_small.csv", row.names = F, col.names = F)



### LARGE RHO 0.2

n_sims <- 10000
n_timesteps <- 30

times <- 1:(n_timesteps+2*30)
alpha <- 4
beta <- 4
rho <- 0.2
kappa <- 3.1131

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
exceeds <- matrix(0, nrow = length(times)-2*30, ncol = n_sims)
for(i in 1:n_sims){
  exceeds[,i] <- rltrawl(alpha = alpha,
                         beta = beta,
                         kappa = 0.0,
                         trawl_fs = trawl_1,
                         trawl_fs_prim = trawl_1_prim,
                         times = times,
                         n = 1,
                         transformation = F)
  print(i)
}

write.table(x = t(exceeds), file = "acf_estimation_sim_10000.csv", row.names = F, col.names = F)




compute_acf <- function(data_mat, alpha, beta, kappa, p = 0.1){
  d <- dim(data_mat)[1]
  return(vapply(1:d, function(i){
                      lagg <- exp(-kappa*(data_mat[1,]+data_mat[i,]))/(data_mat[1,]*data_mat[i,])
                      
                      return(c(mean(lagg)-mean(exp(-kappa*data_mat[1,]))^2, sd(lagg)))
    
                      },c(1.0,1.0)))
}

compute_acf_trawl <- function(data_mat, alpha, beta, kappa){
  d <- dim(data_mat)[1]
  return(vapply(1:d, function(i){
    #lagg <- data_mat[1,]*data_mat[i,]
    return(c(cov(data_mat[1,], data_mat[i,]), sd(data_mat[1,]*data_mat[i,])))
    #return(c(mean(lagg) - mean(data_mat[1,])^2, sd(lagg)))
    
  },c(1.0,1.0)))
}

exceeds <- t(read.table(file = "acf_estimation_sim_10000.csv", header = F))
exceeds_small <- t(read.table(file = "acf_estimation_sim_10000_small.csv", header = F))

par(mfrow = (c(1,2)))

alpha <- 4
beta <- 4
rho <- 0.2
kappa <- 3.1131
plotted <- compute_acf(exceeds, alpha, beta, kappa)
plot(0:29, plotted[1,]/plotted[1,1], type = 'l', xlab = "Lags", ylab = "ACF", main = "Exceedance ACF", col = "black", xlim = c(0,29), ylim=c(0,1))
lines(0:29, plotted[1,]/plotted[1,1] + 1.96 * plotted[2,]/sqrt(n_sims)/plotted[1,1], lty = 2, col = "blue")
lines(0:29, plotted[1,]/plotted[1,1] - 1.96 * plotted[2,]/sqrt(n_sims)/plotted[1,1], lty = 2, col = "blue")

alpha <- 4
beta <- 1
rho <- 0.05
kappa <- 0.75
plotted_small <- compute_acf(exceeds_small, alpha, beta, kappa, p = 0.05)
lines(0:29, plotted_small[1,]/plotted_small[1,1], type = 'l', xlab = "Lags", ylab = "ACF", main = "Exceedance ACF", col = "black", xlim = c(0,29))
lines(0:29, plotted_small[1,]/plotted_small[1,1] + 1.96 * plotted_small[2,]/sqrt(n_sims)/plotted_small[1,1], lty = 2, col = "blue")
lines(0:29, plotted_small[1,]/plotted_small[1,1] - 1.96 * plotted_small[2,]/sqrt(n_sims)/plotted_small[1,1], lty = 2, col = "blue")



## second plot
alpha <- 4
beta <- 4
rho <- 0.2
kappa <- 3.1131
plotted <- compute_acf_trawl(exceeds, alpha, beta, kappa)
plot(0:29, exp(-0.2*0:29), col = "red", xlab = "Lags", ylab = "ACF", main = paste("Trawl ACF"), type = "l", lwd = 2, ylim = c(-0.05,1))
lines(0:29, plotted[1,]/plotted[1,1], type = 'l', lwd = 1)
lines(0:29, plotted[1,]/plotted[1,1] + 1.96 * plotted[2,]/sqrt(n_sims)/plotted[1,1], lty = 2, col = "blue")
lines(0:29, plotted[1,]/plotted[1,1] - 1.96 * plotted[2,]/sqrt(n_sims)/plotted[1,1], lty = 2, col = "blue")

alpha <- 4
beta <- 1
rho <- 0.05
kappa <- 0.75
plotted_small <- compute_acf_trawl(exceeds_small, alpha, beta, kappa)
lines(0:29, exp(-0.1*0:29), col = "red", xlab = "Lags", ylab = "ACF", main = paste("Trawl ACF"), type = "l", lwd = 2, ylim = c(-0.05,1))
lines(0:29, plotted_small[1,]/plotted_small[1,1], type = 'l', lwd = 1)
lines(0:29, plotted_small[1,]/plotted_small[1,1] + 1.96 * plotted_small[2,]/sqrt(n_sims)/plotted_small[1,1], lty = 2, col = "blue")
lines(0:29, plotted_small[1,]/plotted_small[1,1] - 1.96 * plotted_small[2,]/sqrt(n_sims)/plotted_small[1,1], lty = 2, col = "blue")

par(mfrow=c(1,1))

plotted
