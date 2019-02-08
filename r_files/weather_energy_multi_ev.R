# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("multi_ev.R")

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
         "wind_direction")

core_energy_data <- getCoreData(data = energy_weather_merged, 
                    ignore_tags = tags_to_ignore,
                    ignore_cols = cols_to_ignore,
                    ignore_data_type = types_to_ignore)

dim(energy_weather_merged)
dim(core_energy_data)

# dates is a vector of datetime
# data is a vector of data 

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("prep_univariate_latent_trawl_fit.R")
library("evir")

horizon <- c(1,2,3,6,12,24)
s.sample <- 40000
val.params <- findUnivariateParams(data = core_energy_data[,100:102], clusters_size = c(5,5,5))

cont_mat <- makeConditionalMatrices(data = core_energy_data[,100:102],
                                    p.zeroes = 0.95,
                                    horizon = horizon,
                                    params = val.params,
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

list_of_list_horizons_vines_loaded <- list.load("cond-mat-vines-12361224-v2.RData")
list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")

computeTRONwithLists <- function(exceedances, list_vines, list_of_matrix, N=100000){
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
      te.st <- RVineSim(RVM = list_vines[[h]][[i]], N = N)
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

makeRdmTimestamp <- function(){
  striped_time <- strptime(Sys.time(), 
                           format =  "%Y-%m-%d %H:%M:%S")
  elapsed <- as.integer(proc.time()[3])
  pasting <- sapply(X = c(year, month, day, hour, minute, second),
               FUN = function(x){x(striped_time)})
  pasting <- paste(pasting, collapse = "-")
  return(pasting %>% as.character)
}

makeRdmTimestamp()

makeFileName <- function(file_name, tag, extension){
  if(file_name %>% is.na){
    return(paste(makeRdmTimestamp(), "_matrix", ".RData", sep=""))
  }else{
    return(paste(file_name, ".RData", sep=""))
  }
}

makeFileName("ok", "_matrix", "RData")

#' computeTRON allows to compute TRON probabilities very easily!
computeTRON <- function(data, q.s, horizons, clusters, n_samples,
                        name_matrices_file=NA, name_vine_file=NA, name_tron_file=NA){
  name_matrices_file <- makeFileName(name_matrices_file, 
                                     tag = "_matrix",
                                     extension = ".RData")
  name_vine_file <- makeFileName(name_vine_file, 
                                     tag = "_vines",
                                     extension = ".RData")
  name_tron_file <- makeFileName(name_tron_file, 
                                     tag = "_tron",
                                     extension = ".RData")

  thresholds <- getThresholds(data = data, 
                              p.exceed = q.s)
  
  univ.params <- findUnivariateParams(data = data, 
                                      clusters_size = clusters)
  exceedances <- makeExceedances(data = data,
                                 thresholds = thresholds,
                                 normalize = TRUE)
  # compute the matrices
  list_of_mat <- makeConditionalMatrices(data = data,
                                      exceedeances = exceedances,
                                      q.s = q.s,
                                      horizon = horizons,
                                      params = univ.params,
                                      n_samples = n_samples,
                                      name = name,
                                      save = F)
  list_of_mat <- rlist::list.save(list_of_mat, name_matrices_file) # save
  
  # compute the vines
  list_vines <- fitExceedancesVines(exceedances = exceedances,
                                    list_of_matrix = list_of_mat)
  list_vines <- rlist::list.save(list_vines, name_vine_file) #save
  
  # compute TRON
  tron <- computeTRONwithLists(exceedances = exceedances, 
                               list_vines = list_vines,
                               list_of_matrix = list_of_mat)
  rlist::list.save(tron, name_tron_file) # save
  
  return(tron)
}

computeTRON(data = core_energy_data[,100:105],
            q.s = rep(0.95, 6),
            horizons = c(1,2,3),
            clusters = rep(5, 5),
            n_samples = 10000)









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
