setwd("~/GitHub/multi-trawl-extremes/r_files/")
source('latent_novel_simulation.R')
source('pairwise_latent_trawl.R')
source("imp_latent_noven_final.R")
source("latent_novel_simulation.R")
source("imp_warren.R")

warren_1 <- read.table("latent_4_4_3.113118_0.2_1000_v3.csv", header = F)
acf(warren_1[,1])


mom_l1 <- initial_guess_trawl(warren_1)
index<-2
n_timestamps <- 4400  
n_sims <- 500
mat_res <- matrix(0, nrow = n_sims, ncol = 4)
dim(mat_res)

for(index in 1:n_sims){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(warren_1[1:n_timestamps,index][warren_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(warren_1[1:n_timestamps,index] > 0))/length(warren_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*params_to_work_with[1]
  params_to_work_with[3] <- if(mom_l1$rho[index] != Inf) mom_l1$rho[index] else mom_l1$mean_rho
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- function(params){return(-pl_univ_warren(times = 1:n_timestamps,
                                values = warren_1[1:n_timestamps, index],
                                delta = 4,
                                params = params,
                                model_vars_names = univ_model_vars_names,
                                fixed_names = c(),
                                fixed_params = c(),
                                logscale = T,
                                transformation = F))}
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(3,0.5,0.01,exp(0.2))
  upper_b <- c(10,10,exp(-0.01),5)
  
  params_to_work_with
  lower_b
  upper_b
  print(index)
  print(params_to_work_with)
  system.time(res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, factr=1e12),
                           method = "L-BFGS-B", lower = lower_b, upper=upper_b))
  
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}
