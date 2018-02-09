setwd("~/GitHub/multi-trawl-extremes/r_files/")
source("imp_latent_noven_final.R")
source("latent_novel_simulation.R")

alpha <- 6
beta <- 20
rho <- 0.4
kappa <- 20
times_val <- 1:1000
trawl_val <- trawl_exp(t, rho)
trawl_prim <- trawl_exp_primitive(t, rho)

rltrawl(alpha = alpha, 
         beta = beta, 
         kappa = kappa,
         times = times_val,
         n = 1,
         trawl_f = trawl_val,
         trawl_f_prim = trawl_prim,
        transformation = T)
