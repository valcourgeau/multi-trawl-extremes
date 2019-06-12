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
core_energy_data_std <- apply(core_energy_data, MARGIN = 2,
                          FUN = function(x){(x)/sd(x)})
stl_clean_aep <- stl(core_energy_data_std[,88] %>% (function(x){ts(x, frequency = 24)}), s.window = 1, t.window = 24*365/4)
plot(stl_clean_aep)
par(mfrow=c(4,1),mar=c(4.1,5.1,0.1,0.1))
{
plot((1:length(stl_clean_aep$time.series[1:18000,1]))/24,
     ts(apply(stl_clean_aep$time.series[1:18000,] %>% as.matrix, 
                         FUN = sum, MARGIN = 1), frequency = 24) %>% as.vector, 
     main="", yaxt="n", ylab='Data', xaxt='n', xlab='Time (days)',
     type='l', cex.lab=1.8)
axis(2,cex.lab=3, cex.axis=1.7, ylab='Data')
axis(1,cex.lab=3, cex.axis=1.7, ylab='Time (days)')
for(v_x in seq(0,2500,by=365)){abline(v=v_x, lty=2, col='darkgrey')}

plot((1:length(stl_clean_aep$time.series[1:18000,1]))/24,
     stl_clean_aep$time.series[1:18000,1] %>% as.vector, 
     main="", yaxt="n", ylab='Seasonal', xaxt='n', xlab='Time (days)',
     type='l', cex.lab=1.8)
axis(2,cex.lab=3, cex.axis=1.7, ylab='Seasonal')
axis(1,cex.lab=3, cex.axis=1.7, ylab='Time (days)')
for(v_x in seq(0,2500,by=365)){abline(v=v_x, lty=2, col='darkgrey')}

plot((1:length(stl_clean_aep$time.series[1:18000,1]))/24,
     stl_clean_aep$time.series[1:18000,2] %>% as.vector, 
     main="", yaxt="n", ylab='Trend', xaxt='n', xlab='Time (days)',
     type='l', cex.lab=1.8)
axis(2,cex.lab=3, cex.axis=1.7, ylab='Trend')
axis(1,cex.lab=3, cex.axis=1.7, xlab='Time (days)')
for(v_x in seq(0,2500,by=365)){abline(v=v_x, lty=2, col='darkgrey')}

plot((1:length(stl_clean_aep$time.series[1:18000,1]))/24,
     stl_clean_aep$time.series[1:18000,3] %>% as.vector, 
     main="", yaxt="n", ylab='Remainder', xaxt='n', xlab='Time (days)',
     type='l', cex.lab=1.8)
axis(2,cex.lab=3, cex.axis=1.7, ylab='Remainder')
axis(1,cex.lab=3, cex.axis=1.7, xlab='Time (days)')
for(v_x in seq(0,2500,by=365)){abline(v=v_x, lty=2, col='darkgrey')}
}

clean_energy_data <- datasetSTLCleaning(data = core_energy_data, 
                                        frequency = 24,
                                        trend_window = round(24*365/4),
                                        season_window = 24)
clean_energy_data <- datasetCleaning(data = clean_energy_data, dates = dates)

tags_west_coast <- c(
  "Vancouver" ,   "Portland"    , "Seattle"   ,   "Angeles"  ,    "Diego"    ,    "Vegas"    ,    "Phoenix",  "Francisco"   ,
  "Albuquerque"  ,"Denver"      , "Antonio")

tags_west_coast_light <- c(
  "Vancouver" ,   "Portland", "PJMW", "PJME", "York")

tags_east_coast <- c(   "Minneapolis" , "Louis", "Nashville" ,   "Indianapolis", "Atlanta"  ,    "Detroit",   
"Pittsburgh" ,  "Toronto"    ,  "Philadelphia" ,"York"       ,  "Boston", "PJMW", "PJME", "Charlotte") 

tags_east_coast_light <- c(  "Detroit",  "Charlotte", "Toronto",
                      "York", "Vancouver"       ,  "Boston", "AEP", "DUQ") 

clean_east_data <- getCoreData(data = clean_energy_data, 
                               get_tags = tags_east_coast)
save(clean_east_data, file = "clean_east_data.Rda")
clean_west_data <- getCoreData(data = clean_energy_data, 
                               get_tags = tags_west_coast)
save(clean_west_data, file="clean_west_data.Rda")

clean_west_light_data <- getCoreData(data = clean_energy_data, 
                                     get_tags = tags_west_coast_light)
clean_east_light_data <- getCoreData(data = clean_energy_data, 
                                     get_tags = tags_east_coast_light)
# save(clean_east_light_data, file="clean_east_light_data.Rda")
load("clean_east_light_data.Rda")

exc_test <- makeExceedances(data = clean_east_light_data,
                            thresholds = getThresholds(clean_east_light_data, 0.96))

exc_test <- makeExceedances(data = clean_east_light_data,
                            thresholds = getThresholds(clean_east_light_data, 0.96))
set.seed(42)
clean_east_light_data_jittered <- apply(clean_east_light_data, function(x){x+rnorm(length(x), mean = 0, 0.1*sd(x))}, MARGIN = 2)
exc_test_jittered <- makeExceedances(data = clean_east_light_data_jittered,
                            thresholds = getThresholds(clean_east_light_data_jittered, 0.96))
acf(exc_test_jittered[,12])

p.zeroes_guess <- rep(0.96, ncol(clean_east_light_data))
clusters_guess <- ChoosingClusters(clean_east_light_data, p.zeroes_guess)
clusters_guess
horizons_guess <- c(1,2,3,6,12,24,48)
clusters_guess <- rep(4, ncol(clean_east_light_data))

setwd('~/GitHub/multi-trawl-extremes/r_files/')
source('draw_acf.R')
source('multi_ev.R')

for(i in 1:32){
  print(CustomMarginalMLE(exc_test_jittered[,i]))
}








# par(mfrow=c(8,8))
# for(i in 1:64){
#   plot(clean_east_data[1:5000,i], 
#        type='l',
#        ylab=colnames(clean_east_data)[i])
# }
# plot(, type='l')

# colnames(clean_east_data)

# ChoosingThresholds(clean_east_data, 0.96)

tron_east_light <- computeTRON(data = clean_east_light_data %>% as.data.frame() %>% as.matrix(),
                               p.zeroes = p.zeroes_guess,
                               horizons = horizons_guess,
                               clusters = clusters_guess,
                               conditional_on = 1:ncol(clean_east_light_data),
                               n_samples = 40000,
                               save = T,
                               sparse = F,
                               name_matrices_file = paste("matrix_east_light_1_to_72_4th"),
                               name_vine_file = paste("vine_east_light_1_to_72_4th"),
                               name_tron_file = paste("tron_east_ligh_1_to_72_4th"))

for(h in horizons_guess){
  print(t(tron_east_light[[h]]$mean[1:5,19]) %>% (function(x){round(x,2)}))
}

# create the final table

CreateFinalTable <- function(tron_data, horizons, pick = 19){
  tron_to_save <- matrix('o', ncol = length(vars_names), nrow = length(horizons_guess))
  i <- 1
  vars_names <- colnames(tron_east_light[[horizons[1]]]$mean)
  reordered_vars <- c(19:20,1:18, which(grepl('sin', vars_names)), which(grepl('cos', vars_names)))
  for(h in horizons){
    tron_mat <- tron_data[[h]]$mean[pick,]
    tron_mat_sd <- tron_data[[h]]$sd[pick,]
    
    
    
    # reordering
    
    tron_to_save[i,] <- apply(rbind(tron_mat[reordered_vars], tron_mat_sd[reordered_vars]), 
                              FUN = function(x){
                                paste(round(x[1], 2),' (', round(x[2] / sqrt(200), 3), ')', sep = '0')
                                      },
                              MARGIN = 2)
    i <- i + 1
  }
  colnames(tron_to_save) <- vars_names[reordered_vars]
  rownames(tron_to_save) <- horizons
  return(tron_to_save)
}

final_table <- CreateFinalTable(tron_data = tron_east_light,
                 horizons = horizons_guess,
                 pick = 19) %>% t

write.csv(final_table, '~/GitHub/multi-trawl-extremes/results/weather-energy-tron-v3.csv')


indices_show <- c(19, 20, 1:18, 21:32)
mat <- matrix(0, ncol=length(indices_show), nrow=length(horizons_guess))
i <- 1
for(h in horizons_guess){
  mat[i,] <- (t(tron_east_light[[h]]$mean[indices_show, 19]) %>% (function(x){round(x,2)}))
  i <- i + 1
}

rownames(mat) <- horizons_guess

mat
colnames(mat) <- (rownames(tron_east_light[[h]]$mean)[indices_show])

colnames(mat) <- c("AEP", "DUQ", 
                   "Humidity VC", "Humidity DT", "Humidity CH", "Humidity TO", "Humidity NY", "Humidity BO", 
                   "Pressure CV", "Pressure DT", "Pressure CH", "Pressure TO", "Pressure NY", "Pressure BO",
                   "Temperature VC"," Temperature DT"," Temperature CH", "Temperature TO", "Temperature NY", "Temperature BO",
                   "Wind Speed VC cos"," Wind Speed VC sin", "Wind Speed DT cos"," Wind Speed DT sin", "Wind Speed CH cos", "Wind Speed CH sin",
                   "Wind Speed TO cos", "Wind Speed TO sin", "Wind Speed NY cos", "Wind Speed NY sin", "Wind Speed BO cos", "Wind Speed BO sin")

write.csv(t(mat), "prob_1_to_72_2nd.csv")
t(mat)[21:32,]
read.csv("prob_1_to_72_t.csv")[21:32,]
read.csv("prob_1_to_72_t.csv")[21:32,-1] - t(mat)[21:32,]


df_mat <- data.frame(mat)
colnames(df_mat) <- rownames(tron_east_light[[h]]$mean)[indices_show]

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
color_to_use <- getPalette(length(indices_show))
plot(horizons_guess[1:6], df_mat[1:6,2], type = 'o', ylim = c(0,0.35), col=color_to_use[2], lwd=2, xaxt="n", yaxt="n")
for(i in 3:32){
  lines(horizons_guess[1:6], df_mat[1:6,i], col = color_to_use[i], lwd = 2, type='o')
}
axis(side=1,at=horizons_guess,tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0))
axis(side=2,at=seq(0,0.3,by=0.05),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0))

tron_east_light[[1]]$mean %>% (function(x){round(x,2)})
tron_east_light[[2]]$mean %>% (function(x){round(x,2)})
tron_east_light[[3]]$mean %>% (function(x){round(x,2)})

tron_east_light[[1]]$mean[19:20,] %>% (function(x){round(x,2)})


tron_west[[1]]$mean %>% (function(x){round(x,2)})
tron_west[[2]]$mean %>% (function(x){round(x,2)})
tron_west[[3]]$mean %>% (function(x){round(x,2)})

east_v <- rlist::list.load("vine_east_light_1_to_72_vines.RData")
east_m <- rlist::list.load("matrix_east_light_1_to_72_matrix.RData")
for(h in horizons_guess){
  print(VineCopula::RVineGofTest(data = east_m[[h]]$unif.values[[19]],
                                 RVM = east_v[[h]][[19]],
                                 method = "White"))
}



rlist::list.load("2019-2-21-17-8-35_params.RData")

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


cor_matrix <- matrix(0, nrow = 6, ncol = 6)
cor_matrix[1,2] <- -0.33
cor_matrix[1,3] <- -0.31
cor_matrix[1,4] <- -0.27
cor_matrix[1,5] <- -0.65
cor_matrix[1,6] <- -0.59


cor_matrix[2,3] <- 0.35
cor_matrix[2,4] <- 0.28
cor_matrix[2,5] <- 0.35
cor_matrix[2,6] <- 0.33


cor_matrix[3,4] <- 0.38
cor_matrix[3,5] <- 0.49
cor_matrix[3,6] <- 0.40

cor_matrix[4,5] <- 0.43
cor_matrix[4,6] <- 0.43

cor_matrix[5,6] <- 0.76
cor_matrix <- cor_matrix + t(cor_matrix)
diag(cor_matrix) <- 1
colnames(cor_matrix) <- c("O3",
                          "CO",
                          "SO2",
                          "PM10",
                          "NO",
                          "NO2")
rownames(cor_matrix) <- c("O3",
                          "CO",
                          "SO2",
                          "PM10",
                          "NO",
                          "NO2")
corrplot::corrplot(cor_matrix,
                   tl.col = "black",
                   method = "color",
                   order = "hclust",
                   number.cex = 1.7,
                   addCoef.col = "white",
                   addrect = 2)






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
