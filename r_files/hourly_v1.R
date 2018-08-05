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
par(mfrow=c(6,2), mar=c(5.1,4.1,0.5,1.1))
for(i_agent in 3:(n_vars+2)){
  Acf(epd[,i_agent-2], lag.max = 40,
      ylab=paste("ACF", colnames(epd)[i_agent-2]),
      main="")
  plot(epd[,i_agent-2], type="l", ylab=paste("Exceed", colnames(epd)[i_agent-2]),
       xlab="Time")
  #(Pacf(epd[,i_agent-2]))
}
par(mfrow=c(1,1))

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")
s.clusters <- c(10, 15, 11, 11, 13, 13)
val_params <- matrix(0, nrow = length(epd[1,]), ncol = 4)
par(mfrow=c(3,2), mar=c(5.1,4.1,2.1,2.1))
for(i_agent in 1:n_vars){
  #val_params[i_agent,] <- generate_parameters(epd[,i_agent], cluster.size = s.clusters[i_agent])
  evir::qplot(epd[,i_agent][epd[,i_agent] > 0], xi = round(1/val_params[i_agent,1],3), labels = T, main=(colnames(epd)[i_agent]))
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

s.sample <- 2000
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

s.sample <- 2000
vars_names <- colnames(epd)
par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars){
  for(j in 1:n_vars){
    # plot(pgamma(epd.latent[1:s.sample,i], shape = val_params[i,1], rate = val_params[i,2]),
    # pgamma(epd.latent[1:s.sample,j], shape = val_params[j,1], rate = val_params[j,2]), pch=20)
    smoothScatter(epd.latent[1:s.sample,i],
                  epd.latent[1:s.sample,j], 
                  colramp=viridis, xlab = (vars_names[i]), ylab=(vars_names[j]))
  }
}
par(mfrow=c(1,1))


## CONDITIONAL TEST

s.sample <- 5001
vars_names <- colnames(epd)
par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars){
  for(j in 1:n_vars){
    # plot(pgamma(epd.latent[1:s.sample,i], shape = val_params[i,1], rate = val_params[i,2]),
    # pgamma(epd.latent[1:s.sample,j], shape = val_params[j,1], rate = val_params[j,2]), pch=20)
    smoothScatter(pgamma(epd.latent[which(epd_cdf[2:s.sample,i] > q.s[i]), i], shape = val_params[i,1], rate = val_params[i,2]),
                  pgamma(epd.latent[which(epd_cdf[2:s.sample,i] > q.s[i]), j], shape = val_params[j,1], rate = val_params[j,2]), 
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

epd_cdf_ecdf <- epd_cdf
for(i in 1:length(epd[1,])){
  epd_cdf_ecdf[which(epd[,i]==0), i] <- ecdf(pdbl[which(epd[,i]==0), 3+i])(pdbl[which(epd[,i]==0), 3+i]) * p.zeroes[i]
}

colnames(epd_cdf_ecdf) <- colnames(epd)

library("corrplot")
par(mfrow=c(1,3))
corrplot(cor(epd_cdf), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(cor(epd_cdf_ecdf), type = "upper", order = "hclust", method="ellipse",
         tl.col = "black", tl.srt = 45)
corrplot(cor(epd_cdf_ecdf), order = "AOE", method = "color", addCoef.col="white",
         type = "upper", tl.col = "black")
corrplot(cor(epd), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

par(mfrow=c(1,1))
corrplot(cor(epd_cdf_ecdf), order = "AOE", method = "color", addCoef.col="white",
         type = "upper", tl.col = "black")
par(mfrow=c(1,1))


s.sample <- 5000
par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars){
  for(j in 1:n_vars){
    smoothScatter(epd_cdf[1:s.sample,j],
                  epd_cdf[1:s.sample,i],
                  colramp=viridis, xlab = (colnames(epd)[j]), ylab=(colnames(epd)[i]))
  }
}
par(mfrow=c(1,1))

vars_ordered <- c(1,2,6,5,3,4)
s.sample <- 5000
par(mfrow=c(n_vars,n_vars), mar=c(4.4,4.1,0.5,0.5))

horizon <- 12
for(i in 1:n_vars){
  for(j in 1:n_vars){
    #if(j < i){
      #plot.new()
    #}else{
    {
      data_j <- epd_cdf_ecdf[which(epd[1:(s.sample-horizon), vars_ordered[j]] > 0)+horizon, vars_ordered[j]]
      data_i <- epd_cdf_ecdf[which(epd[1:(s.sample-horizon), vars_ordered[j]] > 0)+horizon, vars_ordered[i]]
      
      data_j <- ecdf(data_j)(data_j)
      data_i <- ecdf(data_i)(data_i)
      
      
      smoothScatter(data_j, data_i,
                    colramp=inferno, xlab = (colnames(epd)[vars_ordered[j]]), ylab=(colnames(epd)[vars_ordered[i]]))  
    }
  }
}
par(mfrow=c(1,1))

### CREATION OF MATRICES
list_of_list_horizons <- list()
horizon <- c(1,2,3,6,12,24)
s.sample <- 157710
#s.sample <- 10000

for(h in horizon){
  list_of_matrices_conditional <- list()
  quantile.update.values <- matrix(0, nrow = length(epd[1,]), ncol = length(epd[1,]))
  for(i in 1:n_vars){
    mat_temp <- matrix(0,
                       nrow = length(which(epd[1:(s.sample-h), i] > 0)),
                       ncol = n_vars+1)
    temp <- epd_cdf_ecdf[which(epd[1:(s.sample-h), i] > 0), i]
    mat_temp[,n_vars+1] <- ecdf(temp)(temp)
    
    for(j in 1:n_vars){
        data_j <- epd_cdf_ecdf[which(epd[1:(s.sample-h), i] > 0)+h, j]
        quantile.update.values[i, j] <- mean(data_j <= q.s[j])
        data_j <- ecdf(data_j)(data_j)
        mat_temp[,j] <- data_j 
    }
    
    colnames(mat_temp) <- c(colnames(epd), colnames(epd)[i])
    list_of_matrices_conditional[[i]] <- mat_temp
  }
  colnames(quantile.update.values) <- colnames(epd)
  rownames(quantile.update.values) <- colnames(epd)
  
  list_of_list_horizons[[h]] <- list(unif.values=list_of_matrices_conditional,
                                     quantiles.values=quantile.update.values)
}

require(rlist)
list.save(list_of_list_horizons, file="hourly-bloomsbury-12361224.RData")

horizon <- c(1,2,3,6,12,24)
library(VineCopula)

list_of_list_horizons <- list.load(file = "hourly-bloomsbury-12361224.RData")
list_of_list_horizons_vines <- list()

for(h in horizon){
  list_of_vines_mat <- list()
  cat("Horizon: ", h, "\n")
  for(i in 1:n_vars){
    list_of_vines_mat[[i]] <- RVineStructureSelect(
                                   data = list_of_list_horizons[[h]][[i]], familyset = c(3,4), type = 0,
                                   selectioncrit = "AIC", indeptest = TRUE, level = 0.05,
                                   trunclevel = NA, progress = FALSE, weights = NA, treecrit = "tau",
                                   se = FALSE, rotations = TRUE, method = "mle", cores = 7)
    cat("--->", colnames(epd)[i], "DONE\n")
  }
  list_of_list_horizons_vines[[h]] <- list_of_vines_mat
}
list.save(list_of_list_horizons_vines, file = "hourly-bloomsbury-vines-12361224-v2.RData")

list_of_list_horizons_vines_loaded <- list.load("hourly-bloomsbury-vines-12361224-v2.RData")
list_of_list_horizons <- list.load(file = "hourly-bloomsbury-12361224.RData")

tron_probabilities <- list()
set.seed(42)
for(h in horizon){
  tron_proba_matrix <- matrix(0, nrow = length(epd[1,]), ncol = length(epd[1,]))
  colnames(tron_proba_matrix) <- colnames(epd)
  rownames(tron_proba_matrix) <- colnames(epd)
  tron_proba_matrix_sd <- tron_proba_matrix
  
  for(i in 1:n_vars){
    te.st <- RVineSim(RVM = list_of_list_horizons_vines_loaded[[h]][[i]], N = 100000)
    qq.values <- list_of_list_horizons[[h]]$quantiles.values[i,]
    qq.values <- c(qq.values, NA)
    print(qq.values)
    te.st <- t(apply(te.st, 1, function(x){x>qq.values}))
    te.st <- te.st[,1:6]
    #print(te.st)
    tron_proba_matrix[i,] <- apply(te.st, 2, mean)[1:6]
    tron_proba_matrix_sd[i,] <- (apply(te.st, 2, sd)/sqrt(length(te.st[,1])))[1:6]
  }
  tron_probabilities[[h]] <- list(mean=tron_proba_matrix, sd=tron_proba_matrix_sd)
}
tron_probabilities[[1]]$mean
tron_probabilities[[1]]$sd

list.save(tron_probabilities, file = "hourly-bloomsbury-tron-12361224.RData")

tron_probabilities_N <- list()
N_sims <- 2^(8:17)
for(h in N_sims){
  set.seed(42)
  tron_proba_matrix <- matrix(0, nrow = length(epd[1,]), ncol = length(epd[1,]))
  colnames(tron_proba_matrix) <- colnames(epd)
  rownames(tron_proba_matrix) <- colnames(epd)
  tron_proba_matrix_sd <- tron_proba_matrix
  
  for(i in 1:n_vars){
    te.st <- RVineSim(RVM = list_of_list_horizons_vines_loaded[[6]][[i]], N = h)
    qq.values <- list_of_list_horizons[[6]]$quantiles.values[i,]
    qq.values <- c(qq.values, NA)
    print(qq.values)
    te.st <- t(apply(te.st, 1, function(x){x>qq.values}))
    te.st <- te.st[,1:6]
    #print(te.st)
    tron_proba_matrix[i,] <- apply(te.st, 2, mean)[1:6]
    tron_proba_matrix_sd[i,] <- (apply(te.st, 2, sd)/sqrt(length(te.st[,1])))[1:6]
  }
  tron_probabilities_N[[h]] <- list(mean=tron_proba_matrix, sd=tron_proba_matrix_sd)
}

par(mfrow=c(1,1), mar=c(4.1,4.1,0.5,0.5))
plot(log(N_sims), log(vapply(N_sims, function(i){tron_probabilities_N[[i]]$sd[2,3]}, 1)),
     xlab="Simulations", ylab="Standard deviation")
abline( h = log(seq( -7, -2, 2^{-4})), lty = 3, col = colors()[ 440 ] )
abline( v = seq( 0, 2^18, 2^10), lty = 3, col = colors()[ 440 ] )
line(log(N_sims), log(vapply(N_sims, function(i){tron_probabilities_N[[i]]$sd[2,4]}, 1)))

axis(2, at=x,labels=x, col.axis="red", las=2)



### WORK with TRON probabilities

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/results/MMSEV/")
for(i in 1:n_vars){
  tron_p <- round(t(vapply(1:n_vars, function(h){tron_probabilities[[horizon[h]]]$mean[i,]},
                            rep(0,n_vars))), digits = 3)
  tron_sd <- round(t(vapply(1:n_vars, function(h){tron_probabilities[[horizon[h]]]$sd[i,]},
                            rep(0,n_vars))), digits = 5)
  colnames(tron_p) <- colnames(epd)
  rownames(tron_p) <- horizon
  colnames(tron_sd) <- colnames(epd)
  rownames(tron_sd) <- horizon
  
  write.csv(tron_p, row.names = T, file = paste("tron_",colnames(epd)[i],".csv", sep = ""))
  write.csv(tron_sd, row.names = T, file = paste("tron_",colnames(epd)[i],"_sd.csv", sep = ""))
}

par(mfrow=c(2,3), mar=c(4.1,4.1,0.5,0.5))
plot(list_of_list_horizons_vines_loaded[[1]][[1]], legend.pos="bottomright", type=1, 
     interactive=F, label.bg="white", label.col="black", label.cex = 1.2, edge.lwd=1.15,
     edge.labels=c("family-par"), edge.label.cex=1.5, edge.label.col="blue")


RVineStdError(
 hessian = -RVineHessian(data = RVineSim(N = 10000, RVM = list_of_list_horizons_vines_loaded[[1]][[1]]),
               RVM = list_of_list_horizons_vines_loaded[[1]][[1]])$hessian, 
               RVM=list_of_list_horizons_vines_loaded[[1]][[1]]
  )



par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:n_vars){
  for(j in 1:n_vars){
    data_j <- epd_cdf_ecdf[which(epd[1:(s.sample), vars_ordered[j]] > 0), vars_ordered[j]]
    data_i <- epd_cdf_ecdf[which(epd[1:(s.sample), vars_ordered[j]] > 0), vars_ordered[i]]
    
    data_j <- ecdf(data_j)(data_j)
    data_i <- ecdf(data_i)(data_i)

    if(j < i){
      plot.new()
    }else{
      smoothScatter(data_j,
                    data_i,
                    colramp=inferno, xlab = (colnames(epd)[vars_ordered[j]]), ylab=(colnames(epd)[vars_ordered[i]]))  
    }
  }
}
par(mfrow=c(1,1))

par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
s.sample <- 50000
for(i in 1:n_vars){
  for(j in 1:n_vars){
    data_j <- epd_cdf_ecdf[which(epd[1:s.sample, j] > 0), j]
    data_i <- epd_cdf_ecdf[which(epd[1:s.sample, j] > 0)+2, i]
    print(cor(data_i, data_j, method = "spearman"))
    min_j <- min(data_j)
    max_j <- max(data_j)
    min_i <- min(data_i)
    max_i <- max(data_i)
    smoothScatter((data_j-min_j)/(max_j-min_j),
                  (data_i-min_i)/(max_i-min_i),
                  colramp=viridis, xlab = (colnames(epd)[j]), ylab=(colnames(epd)[i]))
  }
}
par(mfrow=c(1,1))



par(mfrow=c(n_vars,n_vars), mar=c(4.1,4.1,0.5,0.5))
s.sample <- 50000
for(i in 1:n_vars){
  for(j in 1:n_vars){
    data_j <- epd_cdf_ecdf[which(epd[1:s.sample, j] > 0), j]
    data_i <- epd_cdf_ecdf[which(epd[1:s.sample, j] > 0), i]
    #print(cor(data_i, data_j, method = "spearman"))
    min_j <- min(data_j)
    max_j <- max(data_j)
    min_i <- min(data_i)
    max_i <- max(data_i)
    smoothScatter(ecdf((data_j-min_j)/(max_j-min_j))((data_j-min_j)/(max_j-min_j)),
                  ecdf((data_i-min_i)/(max_i-min_i))((data_i-min_i)/(max_i-min_i)),
                  colramp=viridis, xlab = (colnames(epd)[j]), ylab=(colnames(epd)[i]))
  }
}
par(mfrow=c(1,1))