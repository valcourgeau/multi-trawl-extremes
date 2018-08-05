setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")
pd <- read.csv("bloomsbury_1994_1998.csv")
pd <- pd[, c(1,3,5,7,13,15,17)] # extract only numerical values
pd <- cbind(1:length(pd[,1]), pd) # adding timestamps
pd <- pd[which(rowSums(is.na(pd))==0),] # keeps only the full rows without NAs
colnames(pd) <- c("time", "date", "O3", "NO", "NO2", "SO2", "CO", "PM10")
dim(pd)

stlpd <- pd
freq_pd <- c(7500, rep(7500, 5))
for(i_agent in 3:8){
  fitting_matrix <- cbind(cos(2*pi*pd$time/(8700)),
                          sin(2*pi*pd$time/(8700)),
                          cos(2*pi*pd$time/(24)),
                          sin(2*pi*pd$time/(24)),
                          as.numeric(isWeekend(pd$date)==T),
                          vapply(1:12, FUN = function(i){month(as.Date(pd$date)) == i}, FUN.VALUE = wday(as.Date(pd$date))))
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
pd_daily <- stats::aggregate(pd[,-c(1,2)], by=list(pd[,2]), FUN=mean)

library(forecast)
par(mfrow=c(3,2))
for(i_agent in 1:6){
  Acf(pd_daily[,1+i_agent], lag.max = 70)
  #(Pacf(epd[,i_agent]))
}
par(mfrow=c(1,1))