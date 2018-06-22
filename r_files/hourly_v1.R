setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/")
pdbl <- read.csv("hourly_bloomsbury_2000_2017.csv")

# Reducing sample size
pdbl <- pdbl

library(evir)
library(forecast)
library(lubridate)
stlpd <- pdbl
n_vars <- length(pdbl[1,]) - 3
for(i_agent in 4:(n_vars+3)){
  fitting_matrix <- cbind(cos(2*pi*1:length(stlpd$index)/200),
                          sin(2*pi*1:length(stlpd$index)/200),
                          cos(2*pi*1:length(stlpd$index)/14),
                          sin(2*pi*1:length(stlpd$index)/14),
                          #as.numeric(isWeekend(stlpd$date)==T))
                          vapply(1:3, FUN = function(i){quarter(stlpd$date) == i}, FUN.VALUE = quarter(as.Date(stlpd$date))),
  #vapply(1:12, FUN = function(i){month(as.Date(stlpd$date)) == i}, FUN.VALUE = month(as.Date(stlpd$date))),
  vapply(1:6, FUN = function(i){wday(as.Date(stlpd$date)) == i}, FUN.VALUE = wday(as.Date(stlpd$date))),
  vapply(1:23, FUN = function(i){hour(as.Date(lubridate::hm(stlpd$time))) == i}, FUN.VALUE = wday(as.Date(stlpd$date))))
  
  fit <- lm(stlpd[,i_agent] ~ fitting_matrix)
  summary(fit)
  fitting_indices <- which(summary(fit)$coefficients[,4] < 0.05)
  if(1 %in% fitting_indices){
    fitting_indices <- fitting_indices[-1]
  }
  fitting_matrix <- fitting_matrix[,fitting_indices-1]
  print(fitting_indices)
  stlpd[,i_agent] <- lm(stlpd[,i_agent] ~ fitting_matrix)$residuals
}

par(mfrow=c(n_vars,3))
for(i_agent in 4:(n_vars+3)){
  Acf(stlpd[,i_agent], lag.max = 20)
  plot(stlpd[,i_agent], type = "l")
  (Pacf(stlpd[,i_agent]))
}
par(mfrow=c(1,1))

q.s <- rep(0.95, n_vars) #95% everywhere
thr_stl <- rep(0, n_vars)

par(mfrow=c(3,2))
for(i_agent in 4:(n_vars+3)){
  #meplot(stlpd[,i_agent])
  #print(quantile(stlpd[,i_agent], probs = c(0.7, 0.8, 0.9, 0.95, .97)))
  thr_stl[i_agent-3] <- quantile(stlpd[,i_agent], probs = q.s[i_agent-3])[[1]]
}
par(mfrow=c(1,1))

epd <- apply(as.matrix(stlpd[,-c(1:3)]), 
             FUN = function(x){
               (x - thr_stl) * (x > thr_stl)
             }, MARGIN = 1)
epd <- t(epd)
epd <- apply(X = epd, MARGIN = 2, FUN = function(x){return(x/sd(x))})

library(forecast)
par(mfrow=c(6,3))
for(i_agent in 3:(n_vars+2)){
  Acf(epd[,i_agent-2], lag.max = 40)
  plot(epd[,i_agent-2], type="l")
  (Pacf(epd[,i_agent-2]))
}
par(mfrow=c(1,1))

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")
s.clusters <- c(10, 15, 11, 11, 12, 13)
val_params <- matrix(0, nrow = length(epd[1,]), ncol = 4)
par(mfrow=c(3,2), mar=c(5.1,4.1,2.1,2.1))
for(i_agent in 1:n_vars){
  val_params[i_agent,] <- generate_parameters(epd[,i_agent], cluster.size = s.clusters[i_agent])
  qplot(epd[,i_agent][epd[,i_agent] > 0], xi = round(1/val_params[i_agent,1],3), labels = T, main=colnames(epd)[i_agent])
  print((1+val_params[i_agent,4]/val_params[i_agent,2])^{-val_params[i_agent,1]})
}
par(mfrow=c(1,1))
val_params

source("infer_latent_value.R")
epd.latent <- get.latent.values.mat(epd, val_params = val_params[,-3], randomise=F)
plot(epd.latent[,3])

### Fitting Vine Copulas

# latent
library(viridis)
library(ggplot2)
library(ggalt)


s.sample <- 5000
vars_names <- colnames(epd)
par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars){
  for(j in 1:n_vars){
    # plot(pgamma(epd.latent[1:s.sample,i], shape = val_params[i,1], rate = val_params[i,2]),
         # pgamma(epd.latent[1:s.sample,j], shape = val_params[j,1], rate = val_params[j,2]), pch=20)
    smoothScatter(pgamma(epd.latent[1:s.sample,i], shape = val_params[i,1], rate = val_params[i,2]),
                  pgamma(epd.latent[1:s.sample,j], shape = val_params[j,1], rate = val_params[j,2]), 
                  colramp=viridis, xlab = (vars_names[i]), ylab=(vars_names[j]))
  }
}
par(mfrow=c(1,1))

corrplot(cor(epd.latent), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

# 

p.zeroes <- 1-(1+val_params[,4]/val_params[,2])^{-val_params[,1]}
epd_cdf <- apply(X = epd, MARGIN = 1, FUN = function(x){return(plgpd.row(xs = x, p.zeroes = p.zeroes, params.mat = val_params))})
epd_cdf <- t(epd_cdf)

library("corrplot")
corrplot(cor(epd_cdf), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(cor(epd), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)


s.sample <- 5000
par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars){
  for(j in 1:n_vars){
    smoothScatter(epd_cdf[1:s.sample,j],
                  epd_cdf[1:s.sample,i],
                  colramp=viridis, xlab = (vars_names[j]), ylab=(vars_names[i]))
  }
}
par(mfrow=c(1,1))
