# This files contains code to generate the data for simulation study
setwd("~/GitHub/multi-trawl-extremes/r_files/")
source("latent_novel_simulation.R")

# Example
set.seed(42)
n_sims <- 50
times <- 1:800
kappa <- 0.32485
alpha <- 4
beta <- 2
rho <- 0.3
n_moments <- 4

(1+kappa/beta)^{-alpha}

## Find offset scale
offset_shape <- n_moments + 1
kappa / ((1+kappa/beta)^{alpha/offset_shape} - 1)
trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
offset_scale  <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)

cat("Prob non zero for non-trf",(1+kappa/beta)^{-alpha}, "\n")
cat("Prob non zero for trf",(1+kappa/offset_scale)^(-offset_shape), "\n")

## Trawl process simulation
library(gPdtest)

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

system.time(exceed_p <- rlexceed(alpha = alpha,
                              beta = beta,
                              kappa = kappa,
                              times = times,
                              trawl_fs = trawl_1,
                              trawl_fs_prim = trawl_1_prim,
                              n = 10,
                              transformation = F))
exceed_p
compute_mean_std <- function(values_array){
  # values_array contains the time series with first axis as time and second as # of time series
  n_dim <- length(values_array[1,])
  n_values <- length(values_array[,1])
  results <- rep(0, n_dim)
  for(index in 1:n_dim){
    results[index] <- length(which(values_array[,index] > 0))/ n_values
  }
  
  return(list(mean=mean(results), sd=sd(results), prob=results))
}

mom_gpd <- function(values_array){
  # workds under the assumption that alpha > 2
  # values_array contains the time series with first axis as time and second as # of time series
  n_dim <- length(values_array[1,])
  n_values <- length(values_array[,1])
  alphas_mom <- rep(0, n_dim)
  betas_mom <- rep(0, n_dim)
  
  for(index in 1:n_dim){
    var_mom <- var(values_array[,index][values_array[,index]>0])
    mean_mom <- mean(values_array[,index][values_array[,index]>0])
    alphas_mom[index] <- 2*var_mom/(var_mom-mean_mom^2)
    betas_mom[index] <- mean_mom*(alpha_mom-1)
  }
  
  return(list(alpha=alphas_mom, beta=betas_mom, mean_alpha=mean(alphas_mom), mean_beta=mean(betas_mom), sd_alpha=sd(alphas_mom), sd_beta=sd(betas_mom)))
}

compute_mean_std(values_array = exceed_p)
alpha_values <- seq(alpha/5, alpha*3, length.out = 100)
paras_mom <- mom_gpd(exceed_p)

fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c(),
                                                           fixed_params = c(),
                                                           params = c(exp(x), beta+kappa),
                                                           model_vars_names = c("alpha", "beta"),
                                                           logscale = T,
                                                           transformation = F,
                                                           n_moments = 4))}
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c(),
                                                           fixed_params = c(),
                                                           params = c(exp(x[1]), exp(x[2])),
                                                           model_vars_names = c("alpha", "beta"),
                                                           logscale = T,
                                                           transformation = F,
                                                           n_moments = 4))}
res <- optim(par = c(log(paras_mom$alpha[10]), log(paras_mom$beta[10])), fn = fn_to_optim, method = "BFGS")
exp(res$par)
alpha
beta

plot(alpha_values, vapply(alpha_values,
                          function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                                      fixed_names = c("beta"),
                                                                      fixed_params = c(beta),
                                                                      params = c(x),
                                                                      model_vars_names = c("alpha", "beta"),
                                                                      logscale = T,
                                                                      transformation = F,
                                                                      n_moments = 4))},
                          1.0), type = "l", ylab="log-likelihood")
abline(v=alpha, col="red")
abline(v=exp(res$par[1]), col = "blue")


