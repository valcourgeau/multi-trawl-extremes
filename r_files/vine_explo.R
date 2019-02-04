setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("pairwise_latent_trawl.R")
source("imp_latent_noven_final.R")

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")
pd <- read.csv("bloomsbury_1994_1998.csv")
pd <- pd[, c(1,3,5,7,13,15,17)] # extract only numerical values
pd <- cbind(1:length(pd[,1]), pd) # adding timestamps
pd <- pd[which(rowSums(is.na(pd))==0),] # keeps only the full rows without NAs

colnames(pd) <- c("time", "date", "O3", "NO", "NO2", "SO2", "CO", "PM10")
dim(pd)

# % of NAs
1 - 39009 / 43825 # 10.9%

library(lubridate)

# With STL
stlpd <- pd
freq_pd <- c(8700, rep(8700, 5))
for(i_agent in 3:8){
  if(i_agent != 1){
    fitting_matrix <- cbind(cos(2*pi*pd$time/(8700)),
                            sin(2*pi*pd$time/(8700)),
                            cos(2*pi*pd$time/(24)),
                            sin(2*pi*pd$time/(24)),
                            #as.numeric(isWeekend(pd$date)==T),
                            vapply(1:12, FUN = function(i){month(as.Date(pd$date)) == i}, FUN.VALUE = wday(as.Date(pd$date))),
                            vapply(1:7, FUN = function(i){wday(as.Date(pd$date)) == i}, FUN.VALUE = wday(as.Date(pd$date))),
                            vapply(0:23, FUN = function(i){pd$time %% 24 == i}, FUN.VALUE = pd$time))
  }else{
    fitting_matrix <- cbind(cos(2*pi*1:length(pd[,i_agent])/(7500)),
                            sin(2*pi*1:length(pd[,i_agent])/(7500)),
                            cos(2*pi*1:length(pd[,i_agent])/(24)),
                            sin(2*pi*1:length(pd[,i_agent])/(24)))
    
  }
  fit <- lm(pd[,i_agent] ~ fitting_matrix)
  fitting_indices <- which(summary(fit)$coefficients[,4] < 0.05)
  if(1 %in% fitting_indices){
    fitting_indices <- fitting_indices[-1]
  }
  fitting_matrix <- fitting_matrix[,fitting_indices-1]
  print(fitting_indices)
  stlpd[,i_agent] <- lm(pd[,i_agent] ~ fitting_matrix)$residuals
}


par(mfrow=c(6,1), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 3:8){
  plot(as.numeric(stlpd[,i_agent]), type = "l")
}
par(mfrow=c(1,1))

# Selection of threshold
library(evir)

pd <- stlpd

par(mfrow=c(3,2))
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
thr_pd_fit <- c(33, 54, 28, 46.5, 0.89, 39) # 95 90 90 95 95 95

epd <- apply(as.matrix(pd[,-c(1,2)]), 
             FUN = function(x){
               (x - thr_pd_fit) * (x > thr_pd_fit)
             }, MARGIN = 1)
epd <- t(epd)
epd <- apply(X = epd, MARGIN = 2, FUN = function(x){return(x/sd(x))})

par(mfrow=c(3,2))
for(i_agent in 3:8){
  qplot(pd[,i_agent][pd[,i_agent] > thr_pd_fit[i_agent-2]])
}
par(mfrow=c(1,1))

par(mfrow=c(6,2), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 3:8){
  plot(pd[,i_agent], type = 'l', ylab = colnames(pd)[i_agent], xlab = "Timestamp", ylim=c(0,max(pd[,i_agent]*1.4)))
  plot(epd[,i_agent-2], type = 'l', ylab = colnames(pd)[i_agent], xlab = "Timestamp")
}
par(mfrow=c(1,1))



library(forecast)
par(mfrow=c(3,2))
for(i_agent in 1:6){
  Acf(epd[,i_agent], lag.max = 70)
  #(Pacf(epd[,i_agent]))
}
par(mfrow=c(1,1))

clusters_pd <- c(4, 6, 7, 6, 5, 7)

val_params <- matrix(0, nrow = length(epd[1,]), ncol = 4)
par(mfrow=c(3,2), mar=c(5.1,4.1,2.1,2.1))
for(i_agent in 1:6){
  val_params[i_agent,] <- generate_parameters(epd[,i_agent], cluster.size = clusters_pd[i_agent])
  qplot(epd[,i_agent][epd[,i_agent] > 0], xi = round(1/val_params[i_agent,1],3), labels = T, main=colnames(epd)[i_agent])
  print((1+val_params[i_agent,4]/val_params[i_agent,2])^{-val_params[i_agent,1]})
}
par(mfrow=c(1,1))
val_params

### FITTING MODEL
mat_res <- val_params
delta_params <- c(5, 6, 4, 4, 5, 6)

for(i_agent in 2:2){
  params_to_work_with <- val_params[i_agent,]
  fn_to_optim <- loglikelihood_pl_univ_ic(times = pd[,1],
                                          values = as.numeric(epd[,i_agent]),
                                          delta = delta_params[i_agent],
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          # fixed_names = c("beta"),
                                          # fixed_params = c(params_to_work_with[2]),
                                          # fixed_names = c("alpha"),
                                          # fixed_params = c(val_params[i_agent,1]),
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = F)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(0.1,0.2,1e-3,0.2)
  upper_b <- c(10,10,0.2,10)
  
  system.time(res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, ndeps=rep(1e-6,4),
                                                                                                   parscale=c(1,1,100,1), 
                                                                                                   factr=1e14, REPORT=1),
                           method = "L-BFGS-B", lower = lower_b, upper=upper_b))
  print(i_agent)
  print(res$par)
  #mat_res[i_agent,] <- c(val_params[i_agent,1], res$par)
  mat_res[i_agent,] <- res$par
}


zo <- vapply(1:10/10000, function(x){fn_to_optim(c(1.5, 6.45, x, 6.46))}, 1.0)

# [1] 0.1250993
# [1] 0.1038735
# [1] 0.1004896
# [1] 0.0504499
# [1] 0.1063088

p.zeroes <- 1-colSums(epd > 0)/colSums(epd >= 0 | epd < 0)

set.seed(42)

plgpd <- function(x, p.zero, alpha, beta, kappa){
  print(p.zero)
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

