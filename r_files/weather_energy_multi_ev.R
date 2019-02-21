# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source("utils.R")
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

# Now ew process the wind data
tags_to_collect <- GetCityNames(core_energy_data)
core_energy_data <- ConcatAndReplaceWind(core_energy_data, 
                                         tags=tags_to_collect)
dim(energy_weather_merged)
dim(core_energy_data)

# dates is a vector of datetime
# data is a vector of data 
clean_energy_data <- datasetSTLCleaning(data = core_energy_data, 
                                        frequency = 24,
                                        trend_window = 24*365/4,
                                        season_window = 24)
clean_energy_data <- datasetCleaning(data = clean_energy_data, dates = dates)

tags_west_coast <- c(
  "Vancouver" ,   "Portland"    , "Seattle"   ,   "Angeles"  ,    "Diego"    ,    "Vegas"    ,    "Phoenix",  "Francisco"   ,
  "Albuquerque"  ,"Denver"      , "Antonio")

tags_west_coast_light <- c(
  "Vancouver" ,   "Portland", "Seattle"   ,   "Angeles")

tags_east_coast <- c( "Dallas" ,      "Houston"    ,  "City"      ,   "Minneapolis" , "Louis"    ,   
"Chicago"  ,    "Nashville" ,   "Indianapolis", "Atlanta"  ,    "Detroit"    ,  "Jacksonville" ,"Charlotte"  ,  "Miami"   ,    
"Pittsburgh" ,  "Toronto"    ,  "Philadelphia" ,"York"       , "Montreal"   ,  "Boston" ) 

clean_east_data <- getCoreData(data = clean_energy_data, 
                               get_tags = tags_east_coast)
save(clean_east_data, file = "clean_east_data")
clean_west_data <- getCoreData(data = clean_energy_data, 
                               get_tags = tags_west_coast)
save(clean_west_data, file="clean_west_data")
clean_west_light_data <- getCoreData(data = clean_energy_data, 
                                     get_tags = tags_west_coast_light)


p.zeroes_guess <- 0.95
clusters_guess <- ChoosingClusters(clean_west_light_data, p.zeroes_guess)
horizons_guess <- c(1,2,3,4)


tron_west <- computeTRON(data = clean_west_light_data,
                         p.zeroes = p.zeroes_guess,
                         horizons = horizons_guess,
                         clusters = clusters_guess,
                         n_samples = 40000,
                         save = T,
                         sparse = F,
                         name_matrices_file = "matrix_west_light",
                         name_vine_file = "vine_west_light",
                         name_tron_file = "tron_west_light")
tron_west[[1]]$mean %>% (function(x){round(x,2)})
tron_west[[2]]$mean %>% (function(x){round(x,2)})
tron_west[[3]]$mean %>% (function(x){round(x,2)})

rlist::list.load("2019-2-20-11-47-20_params.RData")

vines_fitted <- rlist::list.load("vine_east_vines.RData")
install.packages('ggraph')
library('ggraph')
plot(vines_fitted[[1]][[1]])
plot(vines_fitted[[1]][[5]], var_names = 'legend')


print("Extreme in Pressure.Portland")
for(h in horizons_guess){
  data_to_print <- tron_west[[h]]$mean[3,]
  print(data_to_print %>% (function(x){round(x,2)}))
}

print("Extreme in temperature.Vancouver")
for(h in horizons_guess){
  data_to_print <- tron_west[[h]]$mean[5,]
  print(data_to_print %>% (function(x){round(x,2)}))
}

## Second set of proba

paste("Where do the extreme in", colnames(clean_west_light_data)[9], "come from?")
print(do.call(cbind, lapply(horizons_guess, function(h){tron_west[[h]]$mean[,9]}))) # Look at portland.sin

paste("Where do the extreme in", colnames(clean_west_light_data)[10], "come from?")
print(do.call(cbind, lapply(horizons_guess, function(h){tron_west[[h]]$mean[,10]})))

















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
                                    n_samples = 40000,
                                    name = "conditional-mat-test")

cont_mat %>% print
names(cont_mat)
cont_mat[[1]]$quantiles.values







# Choosing the variables to include
core_energy_data$humidity.Vancouver %>% (function(x){head(x, 1000)}) %>% (function(x){plot(x, type = 'l')})
core_energy_data$PJME_hourly.PJME_MW %>% (function(x){head(x, 1000)}) %>% (function(x){plot(x, type = 'l')})
core_energy_data$humidity.Portland %>% (function(x){head(x, 5000)}) %>% (function(x){plot(x, type = 'l')})



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





cleaned_data <- core_energy_data
cleaned_data[,100:102] <- as.matrix(apply(cleaned_data[,100:102],
                                 MARGIN = 2, 
                                 FUN = function(x){deterministicCleaning(data=x, (1:length(cleaned_data[,1]))/24)}))

tron_test <- computeTRON(data = cleaned_data[,100:102],
                        p.zeroes = rep(0.95, 3),
                        horizons = c(1,24),
                        clusters = rep(5, 3),
                        n_samples = 40000)
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
