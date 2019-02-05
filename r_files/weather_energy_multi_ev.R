# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)

source("multi_ev.R")

# loading data
merged_dataset_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/merged-datasets/"
setwd(merged_dataset_folder)
load("energy-weather-merge.Rda")

ok <- ignoreColStartingWith(data = energy_weather_merged, tags = c("humidity",
                                                                   "wind_direction"))
dim(ok)
colnames(energy_weather_merged)

cols_to_ignore <- c("datetime", "index.x", "index.y")
types_to_ignore <- c("factor")
tags_to_ignore <- c("humidity",
         "wind_direction")

core_energy_data <- getCoreData(data = energy_weather_merged, 
                    ignore_tags = tags_to_ignore,
                    ignore_cols = cols_to_ignore,
                    ignore_data_type = types_to_ignore)

dim(test)

# dates is a vector of datetime
# data is a vector of data 

deterministicCleaning <- function(dates, data, p.adjust.method="bonferroni"){
  dates <- strptime(energy_weather_merged$datetime, "%Y-%m-%d %H:%M:%S", tz = "GMT")
  fitting_matrix <- cbind(cos(2*pi*1:length(data)/24),
                          sin(2*pi*1:length(data)/24),
                          cos(2*pi*1:length(data)/(24*365)),
                          sin(2*pi*1:length(data)/(24*365)),
                          #as.numeric(isWeekend(dates)==T),
                          vapply(1:3, 
                                 FUN = function(i){quarter(dates) == i}, 
                                 FUN.VALUE = quarter(dates)),
                          vapply(1:12, 
                                 FUN = function(i){month(dates) == i}, 
                                 FUN.VALUE = month(dates)),
                          vapply(1:6, 
                                 FUN = function(i){wday(dates) == i}, 
                                 FUN.VALUE = wday(dates)),
                          vapply(1:23, 
                                 FUN = function(i){hour(dates) == i}, 
                                 FUN.VALUE = hour(dates)))
  fit <- lm(data ~ fitting_matrix)
  
  # adjusting the method to take in account multiple p.values testing
  p.v.adjusted <- summary(fit)$coefficients[,4]
  p.v.adjusted <- p.adjust(p.v.adjusted, method = p.adjust.method)
  # TODO test that!
  # print(p.v.adjusted)
  
  fitting_indices <- which(p.v.adjusted < 0.05)
  print(length(fitting_indices))
  if(1 %in% fitting_indices){
    fitting_indices <- fitting_indices[-1]
  }
  fitting_matrix <- fitting_matrix[,fitting_indices-1]
  return(lm(data ~ fitting_matrix)$residuals)
}

datasetCleaning <- function(data, dates){
  result <- apply(data, MARGIN = 2, 
                  FUN = function(x){
                    deterministicCleaning(data = x, dates = dates)
                  }
  )
  result <- as.data.frame(result)
  colnames(result) <- colnames(data)
  return(result)
}

test <- datasetCleaning(data = core_energy_data,
                        dates = energy_weather_merged$datetime)
plot(test$pressure.Miami, type='l')
#rm(test)

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")


threshold_data <- makeExceedances(test, thresholds = getThresholds(test, 0.8))
apply(threshold_data, MARGIN = 2, 
      FUN = function(x){length(which(x > 0.0)) / length(x)})

s.clusters <- rep(5, length(test[1,]))
val.params <- findUnivariateParams(data = threshold_data, clusters_size = s.clusters)

val.params



horizon <- c(1,2,3,6,12,24)
s.sample <- 40000
cont_mat <- makeConditionalMatrices(data = test[,100:102],
                                    exceedeances = threshold_data[,100:102],
                                    q.s=getThresholds(test, 0.95)[100:102],
                                    horizon = horizon,
                                    params = val.params[100:102,],
                                    n_samples = s.sample,
                                    name = "conditional-mat-test")


cont_mat %>% print
names(cont_mat)
cont_mat[[1]]$quantiles.values


############################################################################################
################################## TRON PROBABILITIES ######################################
############################################################################################


### CREATION OF MATRICES

horizon <- c(1,2,3,6,12,24)
require(VineCopula)

list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")
list_of_list_horizons_vines <- list()

n_vars <- length(threshold_data[1,100:102])

fitExceedancesVines <- function(exceedances, list_of_matrix){
  #list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")
  list_of_list_horizons <- list_of_matrix
  list_of_list_horizons_vines <- list()
  
  n_vars <- length(exceedances[1,])
  
  for(h in horizon){
    list_of_vines_mat <- list()
    cat("Horizon: ", h, "\n")
    for(i in 1:n_vars){
      cat("--->", colnames(exceedances)[i], " ")
      list_of_vines_mat[[i]] <- RVineStructureSelect(
        data = list_of_list_horizons[[h]]$unif.values[[i]], familyset = c(3,4), type = 0,
        selectioncrit = "AIC", indeptest = TRUE, level = 0.05,
        trunclevel = NA, progress = FALSE, weights = NA, treecrit = "tau",
        se = FALSE, rotations = TRUE, method = "mle", cores = 7)
      cat("DONE in", "TIME", " \n")
    }
    list_of_list_horizons_vines[[h]] <- list_of_vines_mat
  }
  
  list.save(list_of_list_horizons_vines, file = "cond-mat-vines-12361224-v2.RData")
  return(list_of_list_horizons_vines)
}

fitExceedancesVines(threshold_data[,100:102], list_of_list_horizons)

list_of_list_horizons_vines_loaded <- list.load("cond-mat-vines-12361224-v2.RData")
list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")

computeTRONwithLists <- function(exceedances, list_vines, list_of_matrix){
  tron_probabilities <- list()
  set.seed(42)
  for(h in horizon){
    tron_proba_matrix <- matrix(0, nrow = length(exceedances[1,]), ncol = length(exceedances[1,]))
    colnames(tron_proba_matrix) <- colnames(exceedances)[1:(length(exceedances[1,]))]
    rownames(tron_proba_matrix) <- colnames(exceedances)[1:(length(exceedances[1,]))]
    tron_proba_matrix_sd <- tron_proba_matrix
    cat(paste("Horizon", h, "\n"))
    for(i in 1:n_vars){
      cat(paste("--> extreme in", colnames(tron_proba_matrix)[i]), "...")
      te.st <- RVineSim(RVM = list_vines[[h]][[i]], N = 100000)
      te.st <- te.st[,1:(length(te.st[1,])-1)]
      qq.values <- list_of_matrix[[h]]$quantiles.values[i,]
      qq.values <- c(qq.values)
      te.st <- t(apply(te.st, 1, function(x){x>qq.values}))
      tron_proba_matrix[i,] <- apply(te.st, MARGIN = 2, mean)
      tron_proba_matrix_sd[i,] <- (apply(te.st, MARGIN = 2, sd)/sqrt(length(te.st[,1])))
      cat("\t done\n")
    }
    tron_probabilities[[h]] <- list(mean=tron_proba_matrix, sd=tron_proba_matrix_sd)
  }
  tron_probabilities[[1]]$mean
  tron_probabilities[[1]]$sd
  
  list.save(tron_probabilities, file = "tron-cond-test.RData")
  return(tron_probabilities)
}

tron_temp <- computeTRONwithLists(exceedances = threshold_data[,100:102],
                                   list_vines = list_of_list_horizons_vines_loaded,
                                   list_of_matrix = list_of_list_horizons)

tron_temp[[12]]

computeTRON <- function(data, exceedances, q.s, horizon, params, n_samples, name){
  list_of_mat <- makeConditionalMatrices(data = data,
                                      exceedeances = exceedances,
                                      q.s=q.s,
                                      horizon = horizon,
                                      params = params,
                                      n_samples = n_samples,
                                      name = name)
  # list_of_mat <- list.load("conditional-mat-test.Rdata")
  list_vines <- fitExceedancesVines(exceedances = exceedances,
                                    list_of_matrix = list_of_mat)
  # list_vines <- list.load("cond-mat-vines-12361224-v2.RData")
  tron <- computeTRONwithLists(exceedances = exceedances, 
                               list_vines = list_vines,
                               list_of_matrix = list_of_mat)
  return(tron)
}


max_n_vars <- 110
res <- computeTRON(data = test[,100:max_n_vars],
                    exceedances = threshold_data[,100:max_n_vars],
                      q.s=getThresholds(test, 0.95)[100:max_n_vars],
                      horizon = horizon,
                      params = val.params[100:max_n_vars,],
                      n_samples = s.sample,
                      name = "conditional-mat-test")







# this is to observe the convergence speed
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
