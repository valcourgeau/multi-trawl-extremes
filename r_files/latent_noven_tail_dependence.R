setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('pairwise_latent_trawl.R')
source('latent_novel_simulation.R')

library('evd')

compute_B3_exp_time <- function(rho){
  return(function(t1, t2){
    return(compute_B3_exp(rho = rho, t1 = t1, t2 = t2))
    })
}

compute_B1_exp_time <- function(rho){
  return(function(t1, t2){
    return(compute_B1_exp(rho = rho, t1 = t1, t2 = t2))
  })
}

compute_B_inter_exp_time <- function(rho){
  return(function(t1, t2){
    return(compute_B_inter_exp(rho = rho, t1 = t1, t2 = t2))
  })
}

inv_F_2e <- function(thres_data, u){
  data_to_work_on <- thres_data[thres_data > 0.0]
  return(as.numeric(quantile(x = data_to_work_on, u)))
}

tail_dep_function <- function(thres_data, h, u1, u2, alpha, beta, kappa, compute_B3, compute_B_inter){
  result <- 1+inv_F_2e(thres_data = thres_data, u2) / (beta+2*kappa+inv_F_2e(thres_data = thres_data, u1))
  result <- result^{compute_B_inter(t1 = 0, t2 = h)}
  result <- result * (1+inv_F_2e(thres_data = thres_data, u2) / (beta+kappa))^{compute_B3(t1 = 0, t2 = h)}
  return(result)
}

tail_dep_function_exp <- function(thres_data, u, alpha, beta, kappa, rho){
  return(tail_dep_function(thres_data = thres_data,
                        h = 1,
                        u1 = u,
                        u2 = u,
                        alpha = alpha,
                        beta = beta,
                        kappa = kappa,
                        compute_B3 = compute_B3_exp_time(rho),
                        compute_B_inter = compute_B_inter_exp_time(rho)))
}

# Example

bl_thres[,1]
u <- 0:100/100
alpha <- 6.33
beta <- 20.12
rho <- 0.27
kappa <- 12.2

slsl <- trawl_slice_sets(alpha = alpha,
                         beta = beta,
                         times = 1:1000,
                         trawl_1,
                         trawl_1_prim,
                         1)

plot(inv_F_2e(thres_data = bl_thres[,1], u), inv_F_2e(thres_data = bl_thres[,1], tail_dep_function_exp(thres_data = bl_thres[,1],
                   u = u,
                   alpha = alpha,
                   beta = beta,
                   kappa = kappa,
                   rho = rho)))


evd::clusters(slsl)
