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

system.time(exceed_p_trf <- rlexceed(alpha = alpha,
                                 beta = beta,
                                 kappa = kappa,
                                 times = times,
                                 trawl_fs = trawl_1,
                                 trawl_fs_prim = trawl_1_prim,
                                 n = 10,
                                 transformation = T))



