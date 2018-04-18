setwd("~/GitHub/multi-trawl-extremes/r_files/")
source('latent_novel_simulation.R')
source('pairwise_latent_trawl.R')
source("imp_latent_noven_final.R")
source("latent_novel_simulation.R")
source("imp_warren.R")

# NO TRF
## CASE 1
# General parameters
alpha_true <- 4
beta_true <- 4
rho_true <- 0.2
kappa_true <- 3.113118

latent_1 <- read.csv("latent_4_4_3.113118_0.2_30.csv", sep = " ")
# latent_1 <- latent_1[,1:2]
mom_l1 <- initial_guess_trawl(latent_1)
index<-1
n_timestamps <- 4400 - 1 
n_sims <- 30
mat_res <- matrix(0, nrow = n_sims, ncol = 4)

for(index in 1:n_sims){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(latent_1[1:n_timestamps,index][latent_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(latent_1[1:n_timestamps,index] > 0))/length(latent_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*params_to_work_with[1]
  params_to_work_with[3] <- if(mom_l1$rho[index] != Inf) mom_l1$rho[index] else mom_l1$mean_rho
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = latent_1[1:n_timestamps,index],
                                          delta = 4,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = F)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(3,0.5,exp(-2.5),exp(0.2))
  upper_b <- c(10,10,exp(-0.01),exp(3))
  
  params_to_work_with
  lower_b
  upper_b
  
  res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, pgtol=1e5),#, factr=1e10),
               method = "L-BFGS-B", lower = lower_b, upper=upper_b)
  print(index)
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}



## CASE 2
# General parameters
alpha_true <- 9
beta_true <- 1
rho_true <- 0.05
kappa_true <- 0.75

latent_1 <- read.csv("latent_9_1_0.75_0.05_30.csv", sep = " ")
# latent_1 <- latent_1[,1:2]
mom_l1 <- initial_guess_trawl(latent_1)
index<-1
n_timestamps <- 4400 - 1 
n_sims <- 30
mat_res <- matrix(0, nrow = n_sims, ncol = 4)

for(index in 1:n_sims){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(latent_1[1:n_timestamps,index][latent_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(latent_1[1:n_timestamps,index] > 0))/length(latent_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*params_to_work_with[1]
  params_to_work_with[3] <- if(mom_l1$rho[index] != Inf) mom_l1$rho[index] else mom_l1$mean_rho
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  #print(params_to_work_with)
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = latent_1[1:n_timestamps,index],
                                          delta = 4,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = F)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(4,
               0.1,
               0.01,
               0.2
               )
  upper_b <- c(12,
               15,
               0.8,
               1.0
               )
  
  params_to_work_with
  lower_b
  upper_b
  
  res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, pgtol=1),#, factr=1e10),
               method = "L-BFGS-B", lower = lower_b, upper=upper_b)
  print(index)
  print(params_to_work_with)
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}






#library(DEoptim)
#DEoptim::DEoptim(fn_to_optim, lower = lower_b, upper_b, control = DEoptim.control(trace=T, itermax=200))

fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:500,
                                        values = latent_1[1:500,index],
                                        delta = 4,
                                        lambda = 0.0,
                                        model_vars_names = univ_model_vars_names,
                                        fixed_names = c("beta", "rho"),
                                        fixed_params = c(4, -1.6),#c(params_to_work_with[c(2,3,4)]),
                                        logscale = T,
                                        transformation = F)
plot(seq(0.1,2,length.out = 20), vapply(seq(0.1,2,length.out = 20), fn_to_optim, 1))
n_points_grid <- 20
alpha_s <- seq(2.5,10,length.out = n_points_grid)
beta_s <- seq(0.1,2,length.out = n_points_grid)
vis_mat <- matrix(0, ncol=n_points_grid, nrow=n_points_grid)

for(index_a in 1:n_points_grid){
  for(index_b in 1:n_points_grid){
    vis_mat[index_a,index_b] <- fn_to_optim(c(alpha_s[index_a], beta_s[index_b]))
  }
  print(index_a)
}
library(plot3D)
persp3D(alpha_s, beta_s, vis_mat, theta = 10, phi = 40)

# estimating the ACF
trawl_list_grp <- collection_trawl(times = 1:50, params = list(rho=rho), type = "exp", prim = F)
trawl_list_prim_grp <- collection_trawl(times = 1:50, params = list(rho=rho), type = "exp", prim = T)
tp1 <- 0.0
tp2 <- 0.0
tp3 <- 0.0
for(i in 1:200){
  gen_exc_grp <- rltrawl(alpha = alpha,
                         beta = beta,
                         kappa = kappa,
                         times = 1:50,
                         trawl_fs = trawl_list_grp,
                         trawl_fs_prim = trawl_list_prim_grp,
                         n = 1,
                         transformation = F)
  tp1 <- tp1 + exp(-kappa*(gen_exc_grp[10]+gen_exc_grp[11]))/(gen_exc_grp[10]*gen_exc_grp[11])
  tp2 <- tp2 + exp(-kappa*(gen_exc_grp[10]+gen_exc_grp[12]))/(gen_exc_grp[10]*gen_exc_grp[12])
  tp3 <- tp3 + exp(-kappa*(gen_exc_grp[10]+gen_exc_grp[13]))/(gen_exc_grp[10]*gen_exc_grp[13])
}
tp1/200
tp2/200
tp3/200

rho_univ <- optim(par=params_to_work_with[4],
                 fn = fn_to_optim,
                 control = list(trace=4, maxit=10, pgtol=1e-3, parscale=rep(0.5, 1)),
                 method = "L-BFGS-B",
                 lower=-3,
                 upper=1)
rho_univ$par


gen_exc <- rlexceed(alpha = alpha,
                    beta = beta,
                    kappa = kappa,
                    times = times,
                    trawl_fs = trawl_list,
                    trawl_fs_prim = trawl_list_prim,
                    n = 1,
                    transformation = F)

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

compute_mean_std(values_array = exceed_p)
alpha_values <- seq(alpha/5, alpha*3, length.out = 100)
beta_values <- seq(beta/5, beta*3, length.out = 100)
kappa_values <- seq(kappa/5, kappa*3, length.out = 100)
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
                                                           fixed_names = c("alpha", "beta"),
                                                           fixed_params = c(alpha, beta),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha","beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res <- optim(par = c(0.3*0.8), fn = fn_to_optim, method = "BFGS")
exp(res$par)
alpha
beta

# TRF

## alpha
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("beta", "kappa"),
                                                           fixed_params = c(beta, kappa),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha", "beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res_alpha <- optim(par = log(paras_mom$alpha[10]), fn = fn_to_optim, method = "BFGS")
exp(res_alpha$par)
alpha

plot(alpha_values, vapply(alpha_values,
                          function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                                      fixed_names = c("beta", "kappa"),
                                                                      fixed_params = c(beta, kappa),
                                                                      params = c(x),
                                                                      model_vars_names = c("alpha", "beta", "kappa"),
                                                                      logscale = T,
                                                                      transformation = T,
                                                                      n_moments = 4))},
                          1.0), type = "l", ylab="log-likelihood")
abline(v=alpha, col="red")
abline(v=exp(res_alpha$par[1]), col = "blue")

## beta
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("alpha", "kappa"),
                                                           fixed_params = c(alpha, kappa),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha", "beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res_beta <- optim(par = log(paras_mom$beta[10]), fn = fn_to_optim, method = "BFGS")
exp(res_beta$par)
beta
plot(beta_values, vapply(beta_values,
                          function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                                      fixed_names = c("alpha", "kappa"),
                                                                      fixed_params = c(alpha, kappa),
                                                                      params = c(x),
                                                                      model_vars_names = c("alpha", "beta", "kappa"),
                                                                      logscale = T,
                                                                      transformation = T,
                                                                      n_moments = 4))},
                          1.0), type = "l", ylab="log-likelihood")
abline(v=beta, col="red")
abline(v=exp(res_beta$par[1]), col = "blue")

## kappa
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("alpha", "beta"),
                                                           fixed_params = c(alpha, beta),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha", "beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res_kappa <- optim(par = log(0.8*0.3), fn = fn_to_optim, method = "BFGS")
exp(res_kappa$par)
kappa
plot(kappa_values, vapply(kappa_values,
                         function(x){return(-marginal_gpd_likelihood(values = trf_inv_g(z=trawl_p[,10][trawl_p[,10] > 0.0], alpha = alpha,
                                                                                        beta = beta, kappa = 0.25/0.75, 
                                                                                        offset_scale = trf_find_offset_scale(alpha = alpha, beta = beta, 
                                                                                              kappa = 0.25/0.75, offset_shape = 4)+x,
                                                                                        offset_shape = 4),
                                                                     fixed_names = c("alpha"),
                                                                     fixed_params = c(4),
                                                                     params = c(trf_find_offset_scale(alpha = alpha, beta = beta, 
                                                                                                      kappa = 0.25/0.75, offset_shape = 4)+x),
                                                                     model_vars_names = c("alpha", "beta"),
                                                                     logscale = T,
                                                                     transformation = F,
                                                                     n_moments = 4))},
                         1.0), type = "l", ylab="log-likelihood")
abline(v=kappa, col="red")
abline(v=exp(res_kappa$par[1]), col = "blue")

# Initial guess
## MoM on GPD(alpha, beta+kappa)
## Use this and proba > 0 to get kappa: kappa ~ (beta+kappa)*(1-p^{1/alpha})
## Now we have alpha, beta and kappa individually


setwd("~/GitHub/multi-trawl-extremes/data/simulations/")
e_notrf_4000 <- read.csv(file = "exceed_4000a4b2k0.32485.csv")
dim(e_notrf_4000)
delta <- 4

ig_t_notrf_4000 <- initial_guess_trawl(e_notrf_4000)
params_to_work_with <- c((ig_t_notrf_4000$alpha[1]),
                         (ig_t_notrf_4000$beta[1]),
                         log(ig_t_notrf_4000$rho[1]),
                         log(ig_t_notrf_4000$kappa[1]))

fn_to_optim <- loglikelihood_pl_univ_ic(times = (1:4000),
                                        values = e_notrf_4000[1:4000,1],
                                        delta = delta,
                                        lambda = 0.0,
                                        model_vars_names = univ_model_vars_names,
                                        fixed_names = c("alpha"),
                                        fixed_params = c(4.09),
                                        logscale = T,
                                        transformation = F)

system.time(fn_to_optim(params_to_work_with))

optim(params_to_work_with[c(1)], fn_to_optim, method = "BFGS", 
      hessian = F, control = list(trace=5))$par

optim(params_to_work_with[c(2,3,4)], fn_to_optim, method = "L-BFGS-B", 
      hessian = F, control = list(trace=5), lower = c(0.5*2.5,-3,-3), upper = c(1.5*2.5,2,2))

# diagnostic plots
trial_alpha <- seq(0.01, 1, length.out = 15)
system.time(plot(trial_alpha, vapply(trial_alpha, fn_to_optim, 1)))

xs <- seq(0.1, 2, length.out = 50)
plot(xs, evir::dgpd(xs, xi = 1/5, beta = (0.1)/5))
lines(xs, evir::dgpd(xs, xi = 1/alpha, beta = (beta+kappa)/alpha))


# Warren process
lw_sim <- matrix(0, nrow = length(times), ncol = 10)
for(index in 1:10){
 lw_sim[,index] <- rwprocess(alpha = alpha,
                             beta = beta,
                             kappa = kappa,
                             rho = rho,
                             timesteps = length(times),
                             n=1)
}

lw_mom <- mom_gpd(lw_sim)
acf(lw_sim[,5])
