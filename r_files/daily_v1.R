setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/")

pdbl <- read.csv("daily_cromwell_2000_2017.csv")

library(lubridate)

stlpd <- pdbl
for(i_agent in 3:7){
    fitting_matrix <- cbind(cos(2*pi*1:length(stlpd$index)/200),
                            sin(2*pi*1:length(stlpd$index)/200),
                            #as.numeric(isWeekend(stlpd$date)==T))
                            vapply(1:4, FUN = function(i){quarter(stlpd$date) == i}, FUN.VALUE = quarter(as.Date(stlpd$date))))
                            #vapply(1:12, FUN = function(i){month(as.Date(stlpd$date)) == i}, FUN.VALUE = month(as.Date(stlpd$date))))
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

for(i_agent in 3:7){
  Acf(stlpd[,i_agent], lag.max = 70)
  plot(stlpd[,i_agent])
  (Pacf(stlpd[,i_agent]))
}

q.s <- rep(0.85, 5) #95% everywhere
thr_stl <- rep(0, 5)

par(mfrow=c(3,2))
for(i_agent in 3:7){
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
for(i_agent in 3:7){
  acf(epd[,i_agent-2], lag.max = 10)
  plot(epd[,i_agent-2], type="l")
  (pacf(epd[,i_agent-2]))
}
par(mfrow=c(1,1))

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")

s.clusters <- c(3, 5, 5, 2, 4)
for(i in 1:5){
  print(generate_parameters(data = epd[,i], cluster.size = s.clusters[i]))
}

