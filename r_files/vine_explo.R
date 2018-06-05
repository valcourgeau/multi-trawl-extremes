setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")

pd <- read.csv("bloomsbury_1994_1998.csv")
pd <- pd[, c(1,3,5,7,13,15,17)] # extract only numerical values
pd <- cbind(1:length(pd[,1]), pd) # adding timestamps
pd <- pd[rowSums(is.na(pd))==0,] # keeps only the full rows without NAs
colnames(pd) <- c("time", "date", "O3", "NO", "NO2", "SO2", "CO", "PM10")
dim(pd)

# % of NAs
1 - 39009 / 43825 # 10.9%

# With STL
stlpd <- pd
freq_pd <- c(7500, rep(7500, 5))
for(i_agent in 3:8){
  print(i_agent)
  stl_d <- stl(ts(pd[,i_agent], frequency = freq_pd[i_agent-2]), s.window = "periodic", robust = T)
  stlpd[,i_agent] <- stlpd[,i_agent] - stl_d$time.series[,1] - stl_d$time.series[,2]
  fit <- lm(pd[,i_agent] ~ cos(2*pi*1:length(pd[,i_agent])/(340*24))
     + sin(2*pi*1:length(pd[,i_agent])/(340*24)))
}



par(mfrow=c(6,1), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 3:8){
  plot(as.numeric(stlpd[,i_agent]), type = "l")
}
par(mfrow=c(1,1))

# Selection of threshold
library(evir)

pd <- stlpd

par(mfrow=c(1,6))
for(i_agent in 3:8){
  meplot(pd[,i_agent])
  print(quantile(pd[,i_agent], probs = c(0.7, 0.8, 0.9, 0.95, .97)))
}
par(mfrow=c(1,1))

# 95% 56
# 90% 111
# 90% 96
# 95% 69
# 95% 1.7
# 95% 74

thr_pd <- c(56, 111, 96, 69, 1.7, 74)
thr_pd_stl <- c(38, 66, 26, 54, 0.95, 43)

epd <- apply(as.matrix(pd[,-c(1,2)]), 
             FUN = function(x){
               (x - thr_pd_stl) * (x > thr_pd_stl)
             }, MARGIN = 1)
epd <- t(epd)

par(mfrow=c(3,2))
for(i_agent in 3:8){
  qplot(pd[,i_agent][pd[,i_agent] > thr_pd_stl[i_agent-2]])
}
par(mfrow=c(1,1))


par(mfrow=c(6,2), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 3:8){
  plot(pd[,i_agent], type = 'l', ylab = colnames(pd)[i_agent], xlab = "Timestamp")
  plot(epd[,i_agent-2], type = 'l', ylab = colnames(pd)[i_agent], xlab = "Timestamp")
}
par(mfrow=c(1,1))


require(hypergeo)
zeta <- function(alpha, beta, kappa){
  res.zeta <- (alpha-1) * beta^alpha / ((beta+kappa)^{alpha-1})
  res.zeta <- res.zeta * hypergeo(A = alpha-1, B = alpha-1, C = alpha, z = - beta/(beta+kappa))
  return(Re(res.zeta))
}

get_t1 <- function(alpha, beta, kappa, rho){
  # TODO adapt to other trawl functions
  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta)^{1-alpha} + (2*alpha-3)/((alpha-2)*(alpha-1)))
  d_times <- beta * log(1+kappa/beta)/(alpha*(alpha-1)-1)*(alpha*beta/(alpha-2)*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + alpha/(alpha-2)*zeta(alpha-1,beta,kappa) + zeta(alpha+1,beta,kappa))
  return(-rho * (d_plus - d_times))
}

get_estimate_rho <- function(alpha, beta, kappa, index, data){
  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+2*kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta)^{1-alpha} + (2*alpha-3)/((alpha-2)*(alpha-1)))
  #d_times <- beta * log(1+kappa/beta)/(alpha*(alpha-1)-1)*(alpha*beta/(alpha-2)*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + alpha/(alpha-2)*zeta(alpha-1,beta,kappa) + zeta(alpha+1,beta,kappa))
  
  d_times <- 2*beta^2/((alpha-2)*(alpha-1))*((1+2*kappa/beta)^{2-alpha}*log(1+kappa/beta)+zeta(alpha = alpha-1, beta = beta, kappa = kappa))
  return( - var(data) / (index * alpha * (d_plus-d_times)))
}

generate_parameters <- function(data, cluster.size){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(data[data > 0], u= 0)@fit$fit$par
  p_nz <- length(which(data > 0))/length(data)
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*params_to_work_with[1]
  
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  
  # if(params_to_work_with[1] > 4.5){
  #   al <- params_to_work_with[1]
  #   ratio <- (params_to_work_with[2]+params_to_work_with[4])/params_to_work_with[1]
  #   params_to_work_with[1] <- 4 
  #   params_to_work_with[2] <- ratio * 4 - params_to_work_with[4]
  #   params_to_work_with[4] <- (ratio - params_to_work_with[2]/al)/4
  # }
  # beta_trf <- params_to_work_with[4]/(p_nz^{-1/4}-1)
  # params_to_work_with[3] <- get_estimate_rho(alpha = 4, beta = beta_trf, kappa = params_to_work_with[4], cluster.size, 
  #                                            trf_inv_g(z = latent_1[,i], alpha = params_to_work_with[1], beta = params_to_work_with[2],
  #                                                      kappa = params_to_work_with[4], offset_scale = beta_trf+params_to_work_with[4], offset_shape = 4))
  # params_to_work_with[3] <- min(params_to_work_with[3], 1.0)
  
  params_to_work_with[3] <- get_estimate_rho(alpha = params_to_work_with[1],
                                             beta =  params_to_work_with[2],
                                             kappa = params_to_work_with[4],
                                             cluster.size, 
                                             data)
  #params_to_work_with[3] <- (min(params_to_work_with[3], 1.0))
  
  if(is.na(params_to_work_with[3])){
    params_to_work_with[3] <- runif(min = 0.1,0.5, n = 1)
  }
  return(params_to_work_with)
}

library(forecast)
par(mfrow=c(3,2))
for(i_agent in 1:6){
  Acf(epd[,i_agent])
  #(Pacf(epd[,i_agent]))
}
par(mfrow=c(1,1))

clusters_pd <- c(5, 6, 5, 4, 3, 11)

val_params <- matrix(0, nrow = length(epd[1,]), ncol = 4)
par(mfrow=c(3,2))
for(i_agent in 1:6){
  val_params[i_agent,] <- generate_parameters(epd[,i_agent], cluster.size = clusters_pd[i_agent])
  qplot(epd[,i_agent][epd[,i_agent] > 0], xi = 1/val_params[i_agent,1], labels = T, main=colnames(epd)[i_agent])
  print((1+val_params[i_agent,4]/val_params[i_agent,2])^{-val_params[i_agent,1]})
}
par(mfrow=c(1,1))
val_params

### FITTING MODEL
mat_res <- val_params
delta_params <- rep(4, 7, 4, 4, 4, 6)

for(i_agent in 1:6){
  params_to_work_with <- val_params[i_agent,]
  fn_to_optim <- loglikelihood_pl_univ_ic(times = pd[,1],
                                          values = as.numeric(epd[,i_agent]),
                                          delta = delta_params[i_agent],
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          # fixed_names = c("alpha"),
                                          # fixed_params = c(val_params[i_agent,1]),
                                          logscale = T,
                                          transformation = F)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(2,0.2,exp(-4),0.2)
  upper_b <- c(10,150,exp(-0.01),100)
  
  system.time(res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, factr=1e12),
                           method = "L-BFGS-B", lower = lower_b, upper=upper_b))
  print(i_agent)
  print(res$par)
  mat_res[i_agent,] <- c(val_params[i_agent,1], res$par)
}

# [1] 0.05244943
# [1] 0.1250993
# [1] 0.1038735
# [1] 0.1004896
# [1] 0.0504499
# [1] 0.1063088

p.zeroes <- 1-colSums(epd > 0)/colSums(epd >= 0 | epd < 0)

set.seed(42)

plgpd <- function(x, p.zero, alpha, beta, kappa){
  if(p.zero < 0 | p.zero > 1) stop("p.zero should be between 0 and 1.")
  if(x == 0)
    return(runif(n = 1, min = 0, max = p.zero))
  else{
    return(p.zero + (1-p.zero)*(1-max(0, (1+sign(alpha)*x/(beta+kappa))^{-alpha})))
  }
}

plgpd.row <- function(xs, p.zeroes, params.mat){
  res <- rep(0, length(xs))
  for(i in 1:length(xs)){
    res[i] <- plgpd(x = xs[i],
                    p.zero = p.zeroes[i],
                    alpha = params.mat[i,1],
                    beta = params.mat[i,2],
                    kappa = params.mat[i,4])
  }
  
  return(res)
}

epd_cdf <- apply(X = epd, MARGIN = 1, FUN = function(x){return(plgpd.row(xs = x, p.zeroes = p.zeroes, params.mat = val_params))})
epd_cdf <- t(epd_cdf)

library("corrplot")
corrplot(cor(epd), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

### VINE COPULA
library(VineCopula)
library(CDVine)
epd_vc <- as.copuladata(epd_cdf)
colnames(epd_vc) <- colnames(epd)


fit0 <- RVineStructureSelect(epd_vc, familyset = 3, type = 0,
                             selectioncrit = "logLik", indeptest = TRUE, level = 0.05,
                             trunclevel = NA, progress = TRUE, weights = NA, treecrit = "tau",
                             se = FALSE, rotations = TRUE, method = "mle", cores = 7)
RVineBIC(RVM = fit0, data = epd_vc)
fit1 <- RVineStructureSelect(epd_vc, familyset = 3, type = 1,
                            selectioncrit = "logLik", indeptest = TRUE, level = 0.05,
                            trunclevel = NA, progress = TRUE, weights = NA, treecrit = "tau",
                            se = FALSE, rotations = TRUE, method = "mle", cores = 7)
RVineBIC(RVM = fit1, data = epd_vc)
fit_34_rot <- RVineStructureSelect(epd_vc, familyset = c(3, 4), type = 0,
                            selectioncrit = "logLik", indeptest = TRUE, level = 0.05,
                            trunclevel = NA, progress = TRUE, weights = NA, treecrit = "tau",
                            se = FALSE, rotations = TRUE, method = "mle", cores = 7)
RVineBIC(RVM = fit_34_rot, data = epd_vc)
summary(fit_34_rot)
plot(fit_34_rot)

#### FITTING THE VINES
c_indices_tot <- rowSums(epd > 0) > 0
ci_indices_tot <- which(epd_cdf == apply(epd_cdf, 1, max), arr.ind = TRUE)
ci_indices_tot <- ci_indices_tot[order(ci_indices_tot[,1]),]

# TODO FILTER USING V_X^i
get_ci_indices <- function(i, c_indices, ci_indices, q){
  temp <- c_indices[-c(1:q)][ci_indices[-c(1:q),2] == i & 
                                       c_indices[1:(length(c_indices)-q)]]
  return(temp)
}

# for i = 1
par(mfrow=c(3,2), mar=c(5.1,4.1,4.1,2.1))
for(index in 1:6){
  plot(vapply(1:15, 
              FUN = function(q){return(sum(1*get_ci_indices(i = index,
                                                      c_indices = c_indices_tot,
                                                      ci_indices = ci_indices_tot,
                                                      q = q))
                                       )
              }, FUN.VALUE = 1),
       ylab = "Numbers of extremes",
       xlab = "Lag q", main=paste("Exceedances at t+q given an", colnames(epd)[index], "-dominated extreme at t."))
}
par(mfrow=c(1,1))

q_horizon <- 1
vc_fit <- list()
for(index in 1:1){
  ci_ind <- get_ci_indices(i = index, c_indices = c_indices_tot, ci_indices = ci_indices_tot, q = q_horizon)
  ci_data <- epd_cdf[ci_ind,]
  colnames(ci_data) <- colnames(epd)
  vc_data <- as.copuladata(ci_data)
  fit_temp <- RVineStructureSelect(vc_data, familyset = c(3, 4), type = 0,
                                     selectioncrit = "BIC", indeptest = TRUE, level = 0.05,
                                     trunclevel = NA, progress = TRUE, weights = NA, treecrit = "tau",
                                     se = FALSE, rotations = TRUE, method = "mle", cores = 7)
  
  vc_fit[colnames(epd)[index]] <- list(fit_temp)
}
par(mfrow=c(3,2), mar=c(5.1,4.1,4.1,2.1))
plot(fit_temp)
summary(fit_temp)
RVineBIC(RVM = fit_temp, data = vc_data)
