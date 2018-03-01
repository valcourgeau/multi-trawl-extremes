setwd("~/GitHub/multi-trawl-extremes/r_files/")
source("imp_latent_noven_final.R")
source("latent_novel_simulation.R")
source("imp_warren.R")

# General parameters
set.seed(42)
kappa <- 0.9
alpha <- 5
beta <- 5
rho <- 0.4
n_moments <- 4
times <- 1:1000

(1+kappa/beta)^{-alpha}

## Trawl process simulation
library(gPdtest)

### Generating the functions
trawl_list <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_list_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

# no trf
gen_trawl <- rltrawl(alpha = alpha,
                     beta = beta,
                     times = times,
                     n = 1,
                     trawl_fs = trawl_list,
                     trawl_fs_prim = trawl_list_prim,
                     kappa = kappa,
                     transformation = F)
gen_exc <- rlexceed(alpha = alpha,
                    beta = beta,
                    kappa = kappa,
                    times = times,
                    trawl_fs = trawl_list,
                    trawl_fs_prim = trawl_list_prim,
                    n = 1,
                    transformation = F)

params_to_work_with <- c(alpha, beta, log(rho), log(kappa))
fn_to_optim <- loglikelihood_pl_univ_ic(times = times,
                                        values = gen_exc,
                                        delta = 4,
                                        lambda = 1.0,
                                        model_vars_names = univ_model_vars_names,
                                        fixed_names = c("alpha", "rho", "kappa"),
                                        fixed_params = params_to_work_with[c(1,3,4)],
                                        logscale = T,
                                        transformation = F)
#sample_d <- vapply(X = seq(0.5, 8, length.out = 30), FUN = function(x){fn_to_optim((x))}, 1.0)
# plot(seq(0.5, 8, length.out = 30), sample_d, type ="l",
#      main="PL", xlab="beta", ylab="PL")
# abline(v = beta, col = "red")
plot(seq(0.5, 8, length.out = 30), vapply(X = seq(0.5, 8, length.out = 30), FUN = function(x){fn_to_optim((x))}, 1.0), type ="l",
    main="PL", xlab="rho", ylab="PL")
plot(seq(0.3, 0.5, length.out = 30), vapply(X = seq(0.3, 0.5, length.out = 30), FUN = function(x){fn_to_optim(log(x))}, 1.0) )
plot(density(gen_exc))
hist(gen_exc[gen_exc > 0.1], breaks=20, probability = T)
points(dlgpd(c(0.1,0.8), alpha = 8, beta = 6+kappa))
plot(density(gen_exc[gen_exc > 0]))
lines(seq(0.001,10,length.out = 100), dlgpd(x = seq(0.001,10,length.out = 100), alpha = alpha, beta = beta+kappa))
gPdtest::gpd.fit(gen_exc[gen_exc > 0], method = "amle")

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

# Warren process
lw_sim <- rwprocess(alpha = alpha,
                    beta = alpha,
                    kappa = alpha,
                    rho = rho,
                    timesteps = length(times),
                    n=1)
