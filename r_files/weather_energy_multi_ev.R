# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("multi_ev.R")
source("utils.R")

# loading data
merged_dataset_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/merged-datasets/"
setwd(merged_dataset_folder)
load("energy-weather-merge.Rda")

# ok <- ignoreColStartingWith(data = energy_weather_merged, tags = c("humidity",
#                                                                    "wind_direction"))
# dim(ok)
colnames(energy_weather_merged)

cols_to_ignore <- c("datetime", "index.x", "index.y")
types_to_ignore <- c("factor")
tags_to_ignore <- c(
         "wind_direction",
         "Beersheba",
         "Haifa",
         "Nahariyya",
         "Jerusalem",
         "Eilat",
         "Tel.Aviv")
dates <- strptime(energy_weather_merged$datetime, "%Y-%m-%d %H:%M:%S", tz = "GMT")
#strptime(energy_weather_merged$datetime, "%Y-%m-%d %H:%M:%S", tz = "GMT")

core_energy_data <- getCoreData(data = energy_weather_merged, 
                    ignore_tags = tags_to_ignore,
                    ignore_cols = cols_to_ignore,
                    ignore_data_type = types_to_ignore)

dim(energy_weather_merged)
dim(core_energy_data)

# dates is a vector of datetime
# data is a vector of data 

core_energy_data <- datasetCleaning(data = core_energy_data, dates = dates)


setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")

source("multi_ev.R")
horizon <- c(1,2,3,6,12,24)
s.sample <- 40000
val.params <- findUnivariateParams(data = core_energy_data[,100:105], clusters_size = rep(8, 6),
                                   thresholds = getThresholds(core_energy_data[,100:105], p.exceed = 0.8))


cont_mat <- makeConditionalMatrices(data = core_energy_data[,100:105],
                                    p.zeroes = 0.95,
                                    horizon = horizon,
                                    clusters_size = rep(5, 6),
                                    name = "conditional-mat-test")

cont_mat %>% print
names(cont_mat)
cont_mat[[1]]$quantiles.values


# Choosing the variables to include
core_energy_data$humidity.Vancouver %>% (function(x){head(x, 1000)}) %>% (function(x){plot(x, type = 'l')})
core_energy_data$PJME_hourly.PJME_MW %>% (function(x){head(x, 1000)}) %>% (function(x){plot(x, type = 'l')})
core_energy_data$ %>% (function(x){head(x, 5000)}) %>% (function(x){plot(x, type = 'l')})



############################################################################################
################################## TRON PROBABILITIES ######################################
############################################################################################


### CREATION OF MATRICES

horizon <- c(1,2,3,6,12,24)
require(VineCopula)

list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")
list_of_list_horizons_vines <- list()

n_vars <- length(threshold_data[1,100:105])

list_of_list_horizons_vines_loaded <- list.load("cond-mat-vines-12361224-v2.RData")
list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")

computeTRONwithLists <- function(data, horizons, list_vines, list_of_matrix, N=100000, save=F){
  tron_probabilities <- list()
  set.seed(42)
  n_vars <- length(data[1,])
  for(h in horizons){
    tron_proba_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
    colnames(tron_proba_matrix) <- colnames(data)
    rownames(tron_proba_matrix) <- colnames(data)
    tron_proba_matrix_sd <- tron_proba_matrix
    cat(paste("Horizon", h, "\n"))
    for(i in 1:n_vars){
      cat(paste("--> extreme in", colnames(tron_proba_matrix)[i]), "...")
      te.st <- RVineSim(RVM = list_vines[[h]][[i]], N = N)
      print(paste("min",min(te.st)))
      print(paste("max",min(te.st)))
      
      te.st <- te.st[,1:(length(te.st[1,])-1)]
      qq.values <- list_of_matrix[[h]]$quantiles.values[i,]
      print(qq.values)
      te.st <- t(apply(te.st, MARGIN = 1, FUN = function(x){x>qq.values}))
      print(apply(te.st, MARGIN = 2, mean))
      tron_proba_matrix[i,] <- t(apply(te.st, MARGIN = 2, mean))
      tron_proba_matrix_sd[i,] <- t(apply(te.st, MARGIN = 2, sd))/sqrt(length(te.st[,1]))
      cat("\t done\n")
    }
    tron_probabilities[[h]] <- list(mean=tron_proba_matrix, sd=tron_proba_matrix_sd)
  }
  tron_probabilities[[1]]$mean
  tron_probabilities[[1]]$sd
  if(save){
    list.save(tron_probabilities, file = "tron-cond-test.RData")
  }
  return(tron_probabilities)
}

indices_fit <- 37:72
list_of_list_horizons <- makeConditionalMatrices(data = test[,indices_fit],
                                    exceedeances = threshold_data[,indices_fit],
                                    q.s=getThresholds(test, 0.95)[indices_fit],
                                    horizon = horizon,
                                    params = val.params[indices_fit,],
                                    n_samples = s.sample,
                                    name = "conditional-mat-test")
list_vines <- fitExceedancesVines(threshold_data[,indices_fit], list_of_list_horizons)
tron_temp <- computeTRONwithLists(exceedances = threshold_data[,indices_fit],
                                   list_vines = list_of_list_horizons_vines_loaded,
                                   list_of_matrix = list_of_list_horizons)
plot(threshold_data[1:1000,40], type = 'l')


tron_temp[[24]]$mean

tron_temp[[12]]




cleaned_data <- core_energy_data
cleaned_data[,100:102] <- as.matrix(apply(cleaned_data[,100:102],
                                          MARGIN = 2, 
                                          FUN = function(x){deterministicCleaning(data=x, (1:length(cleaned_data[,1]))/24)}))

tmp <- ChoosingThresholds(data = cleaned_data[,100:102],
                          p.zeroes = c(0.95, 0.95, 0.95))
tmp
plot(tmp$temperature.Montreal$score$proba,
     round(tmp$temperature.Montreal$score$p.values, 3),
     lwd = 2, type="o", pch=23, lty=5)


source("multi_ev.R")
#'computeTRON allows to compute TRON probabilities very easily!
#'@param data clean dataset
#'@param p.zeroes a scalar or vector (as large as the number of columns of
#'  data). of probabilities to be an exceedance of zero (proba of NOT being an
#'  extreme).
#'@param horizons a integer or vector (of integers) of look-ahead horizons for
#'  extremes.
#'@param clusters a integer or vector (of integers) of clusters size in the
#'  autocorrelation sense. See \code{\link[ev.trawl]{GenerateParameters}}.
#'@param n_samples Number of samples to compute the TRON probabilites via
#'  Monte-Carlo.
#'@param name_matrices_file Default is NA. If NA, we use \code{makeFileName(NA,
#'  "_matrix", ".RData")}, otherwise we use \code{name_matrices_file.RData}.
#'@param name_vine_file Default is NA. If NA, we use \code{makeFileName(NA,
#'  "_matrix", ".RData")}, otherwise we use \code{name_vine_file.RData}.
#'@param name_tron_file Default is NA. If NA, we use \code{makeFileName(NA,
#'  "_matrix", ".RData")}, otherwise we use \code{name_tron_file.RData}.
#'@param save Logical (default is TRUE) to save matrices, vines and tron
#'  probabilities as RData files.
#'@return Returns a list of TRON probabilities with horizons as keys.
#'@seealso \code{\link[ev.trawl]{GenerateParameters}} for \code{clusters}.
computeTRON <- function(data, p.zeroes, horizons, clusters, n_samples,
                        name_matrices_file=NA, name_vine_file=NA, 
                        name_tron_file=NA, save=TRUE){
  name_matrices_file <- makeFileName(name_matrices_file, 
                                     tag = "_matrix",
                                     extension = ".RData")
  name_vine_file <- makeFileName(name_vine_file, 
                                     tag = "_vines",
                                     extension = ".RData")
  name_tron_file <- makeFileName(name_tron_file, 
                                     tag = "_tron",
                                     extension = ".RData")

  # Univariate parameters
  exceendances <- makeExceedances(data = data,
                                  thresholds = 
                                    getThresholds(data, p.exceed = p.zeroes),
                                  normalize = TRUE)
  
  # compute the matrices
  list_of_mat <- makeConditionalMatrices(data = data,
                                         p.zeroes = p.zeroes,
                                         horizons = horizons,
                                         clusters_size = clusters,
                                         n_samples = n_samples,
                                         name = name,
                                         save = F)
  if(save){
    rlist::list.save(list_of_mat, name_matrices_file) # save
  }
  
  # compute the vines
  list_vines <- fitExceedancesVines(horizons = horizons, 
                                    list_of_matrix = list_of_mat)
  if(save){
    rlist::list.save(list_vines, name_vine_file) #save
  }
  
  # compute TRON
  tron <- computeTRONwithLists(data = data,
                               horizons = horizons,
                               list_vines = list_vines,
                               list_of_matrix = list_of_mat)
  if(save){
    rlist::list.save(tron, name_tron_file) # save
  }
  return(tron)
}


cleaned_data <- core_energy_data
cleaned_data[,100:102] <- as.matrix(apply(cleaned_data[,100:102],
                                 MARGIN = 2, 
                                 FUN = function(x){deterministicCleaning(data=x, (1:length(cleaned_data[,1]))/24)}))

tron_test <- computeTRON(data = cleaned_data[,100:102],
                        p.zeroes = rep(0.95, 3),
                        horizons = c(1,24),
                        clusters = rep(5, 3),
                        n_samples = 100000)
tron_test[[1]]$mean
tron_test[[24]]$mean







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
    qq.values <- c(qq.values)
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
