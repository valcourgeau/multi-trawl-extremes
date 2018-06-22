setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/")

pdbl <- read.csv("daily_bloomsbury_2000_2017.csv")

library(lubridate)
stlpd <- pdbl
n_vars <- length(pdbl[1,]) - 2
for(i_agent in 3:(n_vars+2)){
    fitting_matrix <- cbind(I(cos(2*pi*1:length(stlpd$index)/200)),
                            I(sin(2*pi*1:length(stlpd$index)/200)),
                            I(cos(2*pi*1:length(stlpd$index)/14)),
                            I(sin(2*pi*1:length(stlpd$index)/14)),
                            #as.numeric(isWeekend(stlpd$date)==T))
                            vapply(1:4, FUN = function(i){quarter(stlpd$date) == i}, FUN.VALUE = quarter(as.Date(stlpd$date))))
                            #vapply(1:12, FUN = function(i){month(as.Date(stlpd$date)) == i}, FUN.VALUE = month(as.Date(stlpd$date))),
                            #vapply(1:7, FUN = function(i){wday(as.Date(stlpd$date)) == i}, FUN.VALUE = wday(as.Date(stlpd$date))))
  
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

par(mfrow=c(3,2))
for(i_agent in 3:(n_vars+2)){
  Acf(pdbl[,i_agent], lag.max = 20)
  #plot(stlpd[,i_agent], type = "l")
  #(Pacf(stlpd[,i_agent]))
}
par(mfrow=c(1,1))

q.s <- rep(0.95, 6) #95% everywhere
thr_stl <- rep(0, 6)

par(mfrow=c(3,2))
for(i_agent in 3:(n_vars+2)){
  meplot(stlpd[,i_agent])
  print(quantile(stlpd[,i_agent], probs = c(0.7, 0.8, 0.9, 0.95, .97)))
  thr_stl[i_agent-2] <- quantile(stlpd[,i_agent], probs = q.s[i_agent-2])[[1]]
}
par(mfrow=c(1,1))

epd <- apply(as.matrix(stlpd[,-c(1,2)]), 
             FUN = function(x){
               (x - thr_stl) * (x > thr_stl)
             }, MARGIN = 1)
epd <- t(epd)
epd <- apply(X = epd, MARGIN = 2, FUN = function(x){return(x/sd(x[x>0]))})

library(forecast)
par(mfrow=c(6,3))
for(i_agent in 3:(n_vars+2)){
  acf(epd[,i_agent-2], lag.max = 10)
  plot(epd[,i_agent-2], type="l")
  (pacf(epd[,i_agent-2]))
}
par(mfrow=c(1,1))

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")
s.clusters <- c(3, 5, 5, 2, 4, 3)
val_params <- matrix(0, nrow = length(epd[1,]), ncol = 4)
par(mfrow=c(3,2), mar=c(5.1,4.1,2.1,2.1))
for(i_agent in 1:n_vars){
  val_params[i_agent,] <- generate_parameters(epd[,i_agent], cluster.size = s.clusters[i_agent])
  qplot(epd[,i_agent][epd[,i_agent] > 0], xi = round(1/val_params[i_agent,1],3), labels = T, main=colnames(epd)[i_agent])
  print((1+val_params[i_agent,4]/val_params[i_agent,2])^{-val_params[i_agent,1]})
}
par(mfrow=c(1,1))
val_params
