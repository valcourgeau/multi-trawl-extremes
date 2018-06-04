setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")

pd <- read.csv("bloomsbury_1994_1998.csv")
pd <- pd[, c(1,3,5,7,13,15,17)] # extract only numerical values
pd <- cbind(1:length(pd[,1]), pd) # adding timestamps
pd <- pd[rowSums(is.na(pd))==0,] # keeps only the full rows without NAs
colnames(pd) <- c("time", "date", "O3", "NO", "NO2", "SO2", "CO", "PM10")
dim(pd)

# % of NAs
1 - 39009 / 43825 # 10.9%

# Selection of threshold
library(evir)

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

par(mfrow=c(1,6))
for(i_agent in 3:8){
  qplot(pd[,i_agent][pd[,i_agent] > thr_pd[i_agent-2]])
}
par(mfrow=c(1,1))

epd <- apply(as.matrix(pd[,-c(1,2)]), 
             FUN = function(x){
                (x - thr_pd) * (x > thr_pd)
             }, MARGIN = 1)
epd <- t(epd)


par(mfrow=c(6,2), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 3:8){
  plot(pd[,i_agent], type = 'l', ylab = colnames(pd)[i_agent], xlab = "Timestamp")
  plot(epd[,i_agent-2], type = 'l', ylab = colnames(pd)[i_agent], xlab = "Timestamp")
}
par(mfrow=c(1,1))


par(mfrow=c(6,1), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 1:6){
  plot(epd[,i_agent], type = 'l', ylab = colnames(pd)[i_agent+2], xlab = "Timestamp")
}
par(mfrow=c(1,1))


# With STL
stlpd <- pd

for(i_agent in 3:8){
  stl_d <- stl(ts(pd[,i_agent], frequency = round(45000/5.5)), "periodic")
  stlpd[,i_agent] <- stlpd[,i_agent] - stl_d$time.series[,1] - stl_d$time.series[,2]
}

par(mfrow=c(6,1), mar=c(4.1,4.5,.6,1.5))
for(i_agent in 3:8){
  plot(as.numeric(stlpd[,i_agent]), type = "l")
}
par(mfrow=c(1,1))


stlpd <- apply(as.matrix(stlpd[,-c(1,2)]), 
             FUN = function(x){
               (x - thr_pd) * (x > thr_pd)
             }, MARGIN = 1)
epd <- t(epd)
stl_d <- stl(ts(pd[,4], frequency = 45000/6), "periodic")
plot(pd[,3]-stl_d$time.series[,1])
meplot(pd[,3]-stl_d$time.series[,1])
