library(VineCopula)
library(CDVine)
setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("pairwise_latent_trawl.R")

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")
pd <- read.csv("bloomsbury_1994_1998.csv")
pd <- pd[, c(1,3,5,7,13,15,17)] # extract only numerical values
pd <- cbind(1:length(pd[,1]), pd) # adding timestamps
pd <- pd[rowSums(is.na(pd))==0,] # keeps only the full rows without NAs
colnames(pd) <- c("time", "date", "O3", "NO", "NO2", "SO2", "CO", "PM10")

#pd_marg <- read.csv("bloombury_marg_params.csv")


stlpd <- pd
freq_pd <- c(7500, rep(7500, 5))
for(i_agent in 3:8){
  fitting_matrix <- cbind(cos(2*pi*1:length(pd[,i_agent])/(7500)),
                          sin(2*pi*1:length(pd[,i_agent])/(7500)),
                          cos(2*pi*1:length(pd[,i_agent])/(24)),
                          sin(2*pi*1:length(pd[,i_agent])/(24)))
  fit <- lm(pd[,i_agent] ~ fitting_matrix)
  fitting_indices <- which(summary(fit)$coefficients[,4] > 0.05)
  print(fitting_indices)
  if(length(fitting_indices) > 0){
    fitting_matrix <- fitting_matrix[,-fitting_indices]
  }
  
  stlpd[,i_agent] <- lm(pd[,i_agent] ~ fitting_matrix)$residuals
}

thr_pd_fit <- c(33, 54, 28, 46.5, 0.89, 39) # 95 90 90 95 95 95

pd <- stlpd
epd <- apply(as.matrix(pd[,-c(1,2)]), 
             FUN = function(x){
               (x - thr_pd_fit) * (x > thr_pd_fit)
             }, MARGIN = 1)
epd <- t(epd)
epd <- apply(X = epd, MARGIN = 2, FUN = function(x){return(x/sd(x[x>0]))})

p.zeroes <- 1-colSums(epd > 0)/colSums(epd >= 0 | epd < 0)
p.zeroes

clusters_pd <- c(6, 11, 12, 10, 7, 13)
delta_params <- c(5, 7, 5, 4, 5, 6)

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
RVineSim(N = 10, RVM = fit0)

c_indices_tot <- rowSums(epd > 0) > 0
ci_indices_tot <- which(epd_cdf == apply(epd_cdf, 1, max), arr.ind = TRUE)
ci_indices_tot <- ci_indices_tot[order(ci_indices_tot[,1]),]

# TODO FILTER USING V_X^i
get_ci_indices <- function(i, c_indices, ci_indices, q){
  temp <- c_indices[-c(1:q)][ci_indices[-c(1:q),2] == i & 
                               c_indices[1:(length(c_indices)-q)]]
  return(temp)
}

q_horizon <- 1
vc_fit <- list()
for(index in 3:3){
  ci_ind <- get_ci_indices(i = index, c_indices = c_indices_tot, ci_indices = ci_indices_tot, q = q_horizon)
  ci_data <- epd_cdf[ci_ind,]
  colnames(ci_data) <- colnames(epd)
  vc_data <- as.copuladata(ci_data)
  fit_temp <- RVineStructureSelect(vc_data, familyset = c(3, 4), type = 0,
                                   selectioncrit = "BIC", indeptest = TRUE, level = 0.05,
                                   trunclevel = NA, progress = TRUE, weights = NA, treecrit = "tau",
                                   se = FALSE, rotations = TRUE, method = "mle", cores = 7)
  sim_temp <- RVineSim(N = 10^6, RVM = fit_temp)
  print(length(which(rowSums(sim_temp > p.zeroes)>0))/10^6)
  vc_fit[colnames(epd)[index]] <- list(fit_temp)
}



transition.lgpd <- function(b,a,alpha,beta,kappa,b1,b_inter,b3){
  # b1 b2 b3 should be chosen with correct h = |t2-t1|
  res <- (1+(2*kappa+a+b)/beta)^{b2}
  res <- res * (1+(kappa+a)/beta)^{alpha+b3}
  res <- res * (1+(kappa+b)/beta)^{b1}
}

simulataneous.lgpd <-function(b,a,alpha,beta,kappa,b1,b_inter,b3){
  # b1 b2 b3 should be chosen with correct h = |t2-t1|
  res <- (1+(2*kappa+a+b)/beta)^{b2}
  res <- res * (1+(kappa+a)/beta)^{b3}
  res <- res * (1+(kappa+b)/beta)^{b1}
}


