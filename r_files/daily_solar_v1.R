setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/")
solar <- read.csv("daily_solar_ams.csv")
station_info <- read.csv("station_info.csv")
solar$Date <- format(as.Date(as.character(solar$Date), "%Y%m%d"), format = "%Y-%m-%d")
wh_stations <- station_info$stid[which(abs(station_info$nlat-35.46)<0.7 
                                       & abs(station_info$elon+97.51)<0.7)]
wh_stations <- c(as.character(wh_stations), "SPEN")
wh_stations
# solar <- solar[,colnames(solar) %in% c("Date",
#                               as.character(wh_stations))]

solar <- solar[,c(1,6:10)]
head(solar)

plot(solar$ACME, type="l")
plot(solar$MARE, type="l")

library(corrplot)
corrplot(cor(solar[,-c(1,2)]), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
dim(solar)

names_vars <- colnames(solar)
solar <- cbind(1:length(solar[,1]), solar)
colnames(solar) <- c("index", names_vars)
  
library(evir)
library(forecast)
library(lubridate)
solcl <- solar
solcl[,-c(1,2)] <- apply(X = solcl[,-c(1,2)], MARGIN = 2, FUN = function(x){return(x/sd(x))})
n_vars_sol <- length(solcl[1,])-2
for(i_agent in 3:(n_vars_sol+1)){
  fitting_matrix <- cbind(cos(2*pi*1:length(solcl$index)/150),
                          sin(2*pi*1:length(solcl$index)/150),
                          cos(2*pi*1:length(solcl$index)/365),
                          sin(2*pi*1:length(solcl$index)/365),
                          cos(2*pi*1:length(solcl$index)/720),
                          sin(2*pi*1:length(solcl$index)/720),
                          #as.numeric(isWeekend(stlpd$date)==T))
                          vapply(1:3, FUN = function(i){quarter(solcl$Date) == i}, FUN.VALUE = quarter(as.Date(solcl$Date))),
                          vapply(1:11, FUN = function(i){month(solcl$Date) == i}, FUN.VALUE = month(as.Date(solcl$Date))))
                          #vapply(1:6, FUN = function(i){wday(as.Date(solcl$Date)) == i}, FUN.VALUE = wday(as.Date(solcl$Date))))
  
  fit <- lm(solcl[,i_agent] ~ fitting_matrix)
  summary(fit)
  fitting_indices <- which(summary(fit)$coefficients[,4] < 0.05)
  if(1 %in% fitting_indices){
    fitting_indices <- fitting_indices[-1]
  }
  fitting_matrix <- fitting_matrix[,fitting_indices-1]
  print(fitting_indices)
  solcl[,i_agent] <- lm(solcl[,i_agent] ~ fitting_matrix)$residuals
}
#solcl[,-c(1,2)] <- apply(X = solcl[,-c(1,2)], MARGIN = 2, FUN = function(x){return(x/sd(x))})
solcl[,-c(1,2)] <- apply(X = solcl[,-c(1,2)], MARGIN = 2, FUN = function(x){return(x/sd(x))})

par(mfrow=c(1,1))
plot(solcl[1:1000,3], type="l")
plot(-solcl[,5], type="l")

library(corrplot)
corrplot(cor(solcl[,-c(1,2)]), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

par(mfrow=c(5,3))
for(i_agent in 3:(5+2)){
  Acf(solcl[,i_agent], lag.max = 20)
  plot(solcl[,i_agent], type = "l")
  (Pacf(solcl[,i_agent]))
}
par(mfrow=c(1,1))

q.s.sol <- rep(0.80, n_vars_sol) #80% everywhere
thr_stl.sol <- rep(0, n_vars_sol)

par(mfrow=c(3,2))
for(i_agent in 3:(n_vars_sol+2)){
  #meplot(solcl[,i_agent])
  #print(quantile(stlpd[,i_agent], probs = c(0.7, 0.8, 0.9, 0.95, .97)))
  thr_stl.sol[i_agent-2] <- quantile(solcl[,i_agent], probs = q.s.sol[i_agent-2])[[1]]
}
par(mfrow=c(1,1))

esol <- apply(as.matrix(solcl[,-c(1:2)]), 
             FUN = function(x){
               (x - thr_stl.sol) * (x > thr_stl.sol)
             }, MARGIN = 1)
esol <- t(esol)
esol <- apply(X = esol, MARGIN = 2, FUN = function(x){return(x/sd(x[x>0]))})

library(forecast)
par(mfrow=c(6,3))
for(i_agent in 3:(6+2)){
  acf(esol[,i_agent-2], lag.max = 40)
  plot(esol[,i_agent-2], type="l")
  (pacf(esol[,i_agent-2]))
}
par(mfrow=c(1,1))

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")
n_vars_sol <- length(colnames(esol))
s.clusters.sol <- rep(5,n_vars_sol)
val_params.sol <- matrix(0, nrow = length(esol[1,]), ncol = 4)
val_params.sol.trf <- val_params.sol

for(i_agent in 1:n_vars_sol){
  val_params.sol[i_agent,] <- generate_parameters(esol[,i_agent], cluster.size = s.clusters.sol[i_agent])
  #print((1+val_params.sol[i_agent,4]/val_params.sol[i_agent,2])^{-val_params.sol[i_agent,1]})
  if(val_params.sol[i_agent,1] < 0){
    val_params.sol.trf[i_agent,] <- val_params.sol[i_agent,]
    val_params.sol.trf[i_agent, 1] <- 4
    val_params.sol.trf[i_agent, 2] <- 1
  }
}

par(mfrow=c(1,2), mar=c(5.1,4.1,2.1,2.1))
for(i_agent in c(1,4)){
  evir::qplot(esol[,i_agent][esol[,i_agent] > 0], xi = round(1/val_params.sol[i_agent,1],3), labels = T, main=colnames(esol)[i_agent])
  #print((1+val_params.sol[i_agent,4]/val_params.sol[i_agent,2])^{-val_params.sol[i_agent,1]})
}
par(mfrow=c(1,1))
val_params.sol
val_params.sol.trf

n_timestamps <- 5113
#n_sims <- 1000
mat_res <- matrix(0, nrow = 5, ncol = 4)

for(index in c(1,3)){
  params_to_work_with <- val_params.sol[index,]

  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = esol[1:n_timestamps,index],
                                          delta = 4,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = T)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(-10,0.1,exp(-2.5),exp(-1.2))
  upper_b <- c(-0.5,10,exp(-0.01),exp(3))
  
  params_to_work_with
  lower_b
  upper_b
  
  system.time(res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)],
                           control = list(trace=3, factr=5e13),
                           method = "L-BFGS-B", lower = lower_b, upper=upper_b))
  print(index)
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}

library(rlist)
list.save(x = list(mat_res), file = "results_univ.RData")

plot(ts(esol[,1], frequency=365.25, start=c(1994,1,1)), type="l", col = alpha("#00ba38", 0.6)
     ,ylab="Daily incoming solar energy exceedances (J.m-2)")
lines(ts(esol[,3], frequency=365.25, start=c(1994,1,1)), type="l", pch = 8, lty=2, col = alpha("#f8766d", 0.6))

par(mfrow=c(1,2))
evir::qplot(esol[,1][esol[,1]>0], xi = round(1/mat_res[1,1],3))
evir::qplot(esol[,3][esol[,3]>0], xi = round(1/mat_res[3,1],3))
par(mfrow=c(1,1))


acf(esol[,1][esol[,1]>0], lag.max = 23)
acf(esol[,3][esol[,3]>0], lag.max = 23)



source("infer_latent_value.R")
esol.latent <- get.latent.values.mat(esol, val_params = val_params.sol.trf[,-3], randomise=F)






### Fitting Vine Copulas

# latent
library(viridis)
library(ggplot2)
library(ggalt)


s.sample <- 5000
vars_names_sol <- colnames(esol)
n_vars_sol <- 6
par(mfrow=c(n_vars_sol,n_vars_sol), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars_sol){
  for(j in 1:n_vars_sol){
    # plot(pgamma(esol.latent[1:s.sample,i], shape = val_params[i,1], rate = val_params[i,2]),
    # pgamma(esol.latent[1:s.sample,j], shape = val_params[j,1], rate = val_params[j,2]), pch=20)
    smoothScatter(pgamma(esol.latent[1:s.sample,i], shape = val_params.sol.trf[i,1], rate = val_params.sol.trf[i,2]),
                  pgamma(esol.latent[1:s.sample,j], shape = val_params.sol.trf[j,1], rate = val_params.sol.trf[j,2]), 
                  colramp=viridis, xlab = (vars_names_sol[i]), ylab=(vars_names_sol[j]))
  }
}

s.sample <- 5001
vars_names <- colnames(esol)
par(mfrow=c(n_vars_sol,n_vars_sol))
for(i in 1:6){
  for(j in 1:6){
    # plot(pgamma(esol.latent[1:s.sample,i], shape = val_params[i,1], rate = val_params[i,2]),
    # pgamma(esol.latent[1:s.sample,j], shape = val_params[j,1], rate = val_params[j,2]), pch=20)
    smoothScatter(pgamma(esol.latent[solcl[2:s.sample-1,3+1] > thr_stl.sol[1], i], shape = val_params.sol.trf[i,1], rate = val_params.sol.trf[i,2]),
                  pgamma(esol.latent[solcl[2:s.sample-1,3+1] > thr_stl.sol[1], j], shape = val_params.sol.trf[j,1], rate = val_params.sol.trf[j,2]), 
                  colramp=viridis, xlab = (vars_names_sol[i]), ylab=(vars_names_sol))
  }
}
par(mfrow=c(1,1))


# Observed
n_vars_sol <- length(colnames(solar))-2
p.zeroes.sol <- 1-(1+val_params.sol[,4]/val_params.sol[,2])^{-val_params.sol[,1]}
p.zeroes.sol <- rep(0.80, n_vars_sol)
esol_cdf <- apply(X = esol, MARGIN = 1, FUN = function(x){return(plgpd.row(xs = x, p.zeroes = p.zeroes.sol, params.mat = val_params.sol))})
esol_cdf <- t(esol_cdf)

esol_cdf_ecdf <- esol_cdf
for(i in 1:length(esol[1,])){
  esol_cdf_ecdf[which(esol[,i]==0), i] <- ecdf(solar[which(esol[,i]==0), 2+i])(solar[which(esol[,i]==0), 2+i])*p.zeroes.sol[i]
  print(summary(ecdf(solar[which(esol[,i]==0), 2+i])(solar[which(esol[,i]==0), 2+i])))
}

library("corrplot")
par(mfrow=c(1,1))
corrplot(cor(esol_cdf), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(cor(esol_cdf_ecdf), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(cor(esol), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

s.sample <- 5000
par(mfrow=c(6,6), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:6){
  for(j in 1:6){
    smoothScatter(esol_cdf[1:s.sample, j],
                  esol_cdf[1:s.sample, i],
                  colramp=viridis, xlab = (colnames(esol)[j]), ylab=(colnames(esol)[i]))
  }
}
par(mfrow=c(1,1))

s.sample <- 5000
par(mfrow=c(6,6), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:6){
  for(j in 1:6){
    smoothScatter(esol_cdf_ecdf[1:s.sample, j],
                  esol_cdf_ecdf[1:s.sample, i],
                  colramp=viridis, xlab = (colnames(esol)[j]), ylab=(colnames(esol)[i]))
  }
}
par(mfrow=c(1,1))

# conditional
{
  

s.sample <- 5001
par(mfrow=c(6,6), mar=c(4.1,4.1,0.5,0.5))
cditonal_on <- 1
for(i in 1:6){
  for(j in 1:6){
    data_j <- esol_cdf[which(esol[1:s.sample, j] > 0), j]
    data_i <- esol_cdf[which(esol[1:s.sample, j] > 0), i]
    min_j <- min(data_j)
    max_j <- max(data_j)
    min_i <- min(data_i)
    max_i <- max(data_i)
    smoothScatter((data_j-min_j)/(max_j-min_j),
                  (data_i-min_i)/(max_i-min_i),
                  colramp=viridis, xlab = (colnames(esol)[j]), ylab=(colnames(esol)[i]))
  }
}
par(mfrow=c(1,1))

s.sample <- 5001
par(mfrow=c(6,6), mar=c(4.1,4.1,0.5,0.5))
cditonal_on <- 1
for(i in 1:6){
  for(j in 1:6){
    data_j <- esol_cdf_ecdf[which(esol[1:s.sample, j] > 0), j]
    data_i <- esol_cdf_ecdf[which(esol[1:s.sample, j] > 0), i]
    min_j <- min(data_j)
    max_j <- max(data_j)
    min_i <- min(data_i)
    max_i <- max(data_i)
    smoothScatter((data_j-min_j)/(max_j-min_j),
                  (data_i-min_i)/(max_i-min_i),
                  colramp=viridis, xlab = (colnames(esol)[j]), ylab=(colnames(esol)[i]))
  }
}
par(mfrow=c(1,1))

s.sample <- 5001
par(mfrow=c(6,6), mar=c(4.1,4.1,0.5,0.5))
cditonal_on <- 1
for(i in 1:6){
  for(j in 1:6){
    data_j <- esol_cdf_ecdf[which(esol[1:s.sample, j] > 0), j]
    data_i <- esol_cdf_ecdf[which(esol[1:s.sample, j] > 0), i]
    min_j <- min(data_j)
    max_j <- max(data_j)
    min_i <- min(data_i)
    max_i <- max(data_i)
    smoothScatter(ecdf((data_j-min_j)/(max_j-min_j))((data_j-min_j)/(max_j-min_j)),
                  ecdf((data_i-min_i)/(max_i-min_i))((data_i-min_i)/(max_i-min_i)),
                  colramp=viridis, xlab = (colnames(esol)[j]), ylab=(colnames(esol)[i]))
  }
}
par(mfrow=c(1,1))
}

par(mfrow=c(9,9), mar=c(4.4,4.1,0.5,0.5))

horizon <- 1
for(i in 1:9){
  for(j in 1:9){
    #if(j < i){
    #plot.new()
    #}else{
    {
      data_j <- esol_cdf_ecdf[which(esol[1:(s.sample-horizon), j] > 0)+horizon, j]
      data_i <- esol_cdf_ecdf[which(esol[1:(s.sample-horizon), j] > 0)+horizon, i]
      
      data_j <- ecdf(data_j)(data_j)
      data_i <- ecdf(data_i)(data_i)
      
      
      smoothScatter(data_j, data_i,
                    colramp=inferno, xlab = (colnames(esol)[j]), ylab=(colnames(esol)[i]))  
    }
  }
}
par(mfrow=c(1,1))

### CREATION OF MATRICES
list_of_list_horizons_sol <- list()
horizon_sol <- c(1,2,3,4,7)
s.sample <- 5000
#s.sample <- 10000

for(h in horizon_sol){
  list_of_matrices_conditional <- list()
  quantile.update.values <- matrix(0, nrow = length(esol[1,]), ncol = length(esol[1,]))
  for(i in 1:n_vars_sol){
    mat_temp <- matrix(0,
                       nrow = length(which(esol[1:(s.sample-h), i] > 0)),
                       ncol = n_vars_sol+1)
    temp <- esol_cdf_ecdf[which(esol[1:(s.sample-h), i] > 0), i]
    mat_temp[,n_vars_sol+1] <- ecdf(temp)(temp)
    
    for(j in 1:n_vars_sol){
      data_j <- esol_cdf_ecdf[which(esol[1:(s.sample-h), i] > 0)+h, j]
      quantile.update.values[i, j] <- mean(data_j <= q.s.sol[j])
      data_j <- ecdf(data_j)(data_j)
      mat_temp[,j] <- data_j 
    }
    
    colnames(mat_temp) <- c(colnames(esol), colnames(esol)[i])
    list_of_matrices_conditional[[i]] <- mat_temp
  }
  colnames(quantile.update.values) <- colnames(esol)
  rownames(quantile.update.values) <- colnames(esol)
  
  list_of_list_horizons_sol[[h]] <- list(unif.values=list_of_matrices_conditional,
                                     quantiles.values=quantile.update.values)
}

require(rlist)
list.save(list_of_list_horizons_sol, file="daily-solar-12347.RData")

library(VineCopula)
list_of_list_horizons_sol <- list.load(file = "daily-solar-12347.RData")
list_of_list_horizons_vines_sol <- list()

for(h in horizon_sol){
  list_of_vines_mat <- list()
  cat("Horizon: ", h, "\n")
  for(i in 1:n_vars_sol){
    list_of_vines_mat[[i]] <- RVineStructureSelect(
      data = list_of_list_horizons_sol[[h]]$unif.values[[i]], familyset = c(3,4), type = 0,
      selectioncrit = "AIC", indeptest = TRUE, level = 0.05,
      trunclevel = NA, progress = FALSE, weights = NA, treecrit = "tau",
      se = FALSE, rotations = TRUE, method = "mle", cores = 7)
    cat("--->", colnames(esol)[i], "DONE\n")
  }
  list_of_list_horizons_vines_sol[[h]] <- list_of_vines_mat
}
list.save(list_of_list_horizons_vines_sol, file = "daily-solar-vines-12347.RData")

list_of_list_horizons_vines_loaded_sol <- list.load("daily-solar-vines-12347.RData")
list_of_list_horizons_sol <- list.load(file = "daily-solar-12347.RData")

tron_probabilities_sol <- list()
set.seed(42)
for(h in horizon_sol){
  tron_proba_matrix <- matrix(0, nrow = length(esol[1,]), ncol = length(esol[1,]))
  colnames(tron_proba_matrix) <- colnames(esol)
  rownames(tron_proba_matrix) <- colnames(esol)
  tron_proba_matrix_sd <- tron_proba_matrix
  
  for(i in 1:n_vars_sol){
    te.st <- RVineSim(RVM = list_of_list_horizons_vines_loaded_sol[[h]][[i]], N = 1000)
    qq.values <- list_of_list_horizons_sol[[h]]$quantiles.values[i,]
    qq.values <- c(qq.values, NA)
    #print(length(qq.values))
    te.st <- t(apply(te.st, 1, function(x){x>qq.values}))
    #te.st <- te.st[,1:13]
    #print(te.st)
    #print(apply(te.st, 2, mean))
    tron_proba_matrix[i,] <- apply(te.st, 2, mean)[1:13]
    tron_proba_matrix_sd[i,] <- (apply(te.st, 2, sd)/sqrt(length(te.st[,1])))[1:13]
  }
  tron_probabilities_sol[[h]] <- list(mean=tron_proba_matrix, sd=tron_proba_matrix_sd)
}
tron_probabilities_sol[[1]]$mean
round(tron_probabilities_sol[[1]]$sd, 3)

list.save(tron_probabilities_sol, file = "daily-solar-tron-12347.RData")