# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)

# loading data
merged_dataset_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/merged-datasets/"
setwd(merged_dataset_folder)
load("energy-weather-merge.Rda")


ignoreColStartingWith <- function(data, tags, return.filters=F){
  # return.filters 
  filtered_by_tag <- sapply(tags, function(x){grepl(x, colnames(data))})
  filtered_by_tag <- apply(filtered_by_tag, 
                           MARGIN = 1, 
                           FUN = function(x){Reduce('|',x)}) # merging the filters
  if(return.filters){
    return(!filtered_by_tag)
  }
  return(data[,!filtered_by_tag]) 
}
ok <- ignoreColStartingWith(data = energy_weather_merged, tags = c("humidity",
                                                                   "wind_direction"))
dim(ok)
colnames(energy_weather_merged)

getCoreData <- function(data, ignore_tags=c(), ignore_cols=c(), ignore_data_type=c()){
  return(data[, ignoreColStartingWith(data, ignore_tags, return.filters = T) &
                (!colnames(data) %in% ignore_cols) & 
                !(sapply(data, class) %in% ignore_data_type)])
}
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
require(stats)
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
s.clusters <- c(10, 15, 11, 11, 13, 13)

findUnivariateParams <- function(data, clusters_size){
  n_vars <- length(data[1,])
  val_params <- matrix(0, nrow = n_vars, ncol = 4)
  for(i_agent in 1:n_vars){
    print(i_agent)
    tryCatch(
      val_params[i_agent,] <- generate_parameters(data[,i_agent],
                                                cluster.size = clusters_size[i_agent])
    , error = function(e) {
      val_params[i_agent,] <- 0
    })
  }
  return(val_params)
}

dim(test)
dim(energy_weather_merged)

getThresholds <- function(data, p.exceed){
  if(p.exceed < 0 | p.exceed > 1)
    stop('p.exceed should be between 0 and 1.')
  
  # deals with data as a vector
  if(is.vector(data)){
    return(quantile(data, p.exceed[1]))
  }
  
  if(length(p.exceed) == 1){
    return(sapply(data, function(x){quantile(x, p.exceed)}))
  }else{
    if(length(data[1,]) != p.exceed){
      return(sapply(data, function(x){quantile(x, p.exceed)}))
    }else{
      stop('p.exceed should either be a scalar or as wide as data.')
    }
  }
}

makeExceedances <- function(data, thresholds, normalize=TRUE){
  if(length(data[1,]) != length(thresholds)){
    stop('thresholds and data have non-comforting sizes. 
         Tip: Use getThresholds to get the right size.')
  }
  
  if(!is.vector(data)){
    n <- length(data[,1])
    rep_thres <- t(matrix(rep(thresholds, n), ncol=n))
    epd <- (data - rep_thres) * (data > rep_thres)
  }else{
    epd <- (data - thresholds) * (data > thresholds)
    
  }
  if(normalize){
    epd <- apply(X = as.matrix(epd), MARGIN = 2, 
                 FUN = function(x){return(x/sd(x))})
    
  }  
  epd <- data.frame(epd)
  colnames(epd) <- colnames(data)
  return(epd)
}

threshold_data <- makeExceedances(test, thresholds = getThresholds(test, 0.8))
getThresholds(test, 0.8)[4]
length(which(test$humidity.Seattle > 100))

threshold_data$humidity.Vancouver
apply(threshold_data, MARGIN = 2, 
      FUN = function(x){length(which(x > 0.0)) / length(x)})

#plot(threshold_data$pressure.Houston)
#plot(energy_weather_merged$pressure.Boston)

s.clusters <- rep(5, length(test[1,]))
val.params <- findUnivariateParams(data = threshold_data, clusters_size = s.clusters)

val.params


plgpd_unif_at_zero <- function(x, p.zero, alpha, beta, kappa){
  if(p.zero < 0 | p.zero > 1) stop("p.zero should be between 0 and 1.")
  if(x == 0)
    return(runif(n = 1, min = 0, max = p.zero))
  else{
    return(p.zero + (1-p.zero)*(1-max(0, (1+sign(alpha)*x/(beta+kappa))^{-alpha})))
  }
}

plgpd_unif_at_zero.row <- function(xs, p.zeroes, params.mat){
  # params.mat contains alpha beta rho kappa
  res <- rep(0, length(xs))
  for(i in 1:length(xs)){
    res[i] <- plgpd_unif_at_zero(x = xs[i],
                    p.zero = p.zeroes[i],
                    alpha = params.mat[i,1],
                    beta = params.mat[i,2],
                    kappa = params.mat[i,4])
  }
  
  return(res)
}

computePZero <- function(params){
  # Computes probability of having zero given univariate extreme value model
  return(1-(1+params[,4]/params[,2])^{-abs(params[,1])})
}

makeConditionalMatrices <- function(data, exceedeances, q.s, 
                                    horizon, save=T, params, n_samples = length(data[,1]),
                                    name="conditional-matrices"){
  
  p.zeroes <- computePZero(params)
  # epd_cdf <- apply(X = exceedeances, MARGIN = 1, 
  #                  FUN = function(x){
  #                    return(plgpd.row(xs = x, p.zeroes = p.zeroes, params.mat = params))
  #                    }
  #                  )
  # epd_cdf <- t(epd_cdf)

  exceedeances_cdf_ecdf <- exceedeances
  for(i in 1:length(exceedeances[1,])){
    # This creates ordered uniform samples in the same order 
    # as in data.
    exceedeances_cdf_ecdf[which(exceedeances[,i]==0), i] <- 
      ecdf(data[which(exceedeances[,i]==0), i])(data[which(exceedeances[,i]==0), i]) * p.zeroes[i]
  }

  list_of_list_horizons <- list()
  n_vars <- length(data[1,])
  for(h in horizon){
    list_of_matrices_conditional <- list()
    # creates a square matrix with all the vars in
    quantile.update.values <- matrix(0, 
                                     nrow = length(exceedeances[1,]), 
                                     ncol = length(exceedeances[1,]))
    colnames(quantile.update.values) <- colnames(exceedeances)
    rownames(quantile.update.values) <- colnames(exceedeances)
    
    for(i in 1:n_vars){
      # creates a temporary matrix with cols equal to number of nvars + 1
      # and rows such that the i-th component is an extreme
      mat_temp <- matrix(0,
                         nrow = length(which(exceedeances[1:(s.sample-h), i] > 0)),
                         ncol = n_vars+1)
      # filtering the i-th eCDF column with extremes 
      temp <- exceedeances_cdf_ecdf[which(exceedeances[1:(s.sample-h), i] > 0), i]
      # addting those values in the nvars + 1 column
      mat_temp[,n_vars+1] <- ecdf(temp)(temp)
      
      for(j in 1:n_vars){
        # filtering data of j-th vars when i-th is extreme h timesteps before
        data_j <- exceedeances_cdf_ecdf[which(exceedeances[1:(s.sample-h), i] > 0)+h, j]
        # computing the probability that j-th was an extreme as well 
        # h timesteps after i-th was an extreme
        quantile.update.values[i, j] <- mean(data_j <= q.s[j])
        data_j <- ecdf(data_j)(data_j)
        # saving the unif values of j-th var used here.
        mat_temp[,j] <- data_j 
      }
      
      colnames(mat_temp) <- c(colnames(exceedeances), colnames(exceedeances)[i])
      list_of_matrices_conditional[[i]] <- mat_temp
    }

    
    list_of_list_horizons[[h]] <- list(unif.values=list_of_matrices_conditional,
                                       quantiles.values=quantile.update.values)
  }
  
  if(save){
    file_name <- paste(name, ".RData", sep="")
    cat(paste("Backing up the conditional matrices as", file_name, "..."))
    list.save(list_of_list_horizons,
              file=file_name)
    cat("done\n")
  }
  return(list_of_list_horizons)
}

horizon <- c(1,2,3,6,12,24)
s.sample <- 40000
cont_mat <- makeConditionalMatrices(data = test[,100:118],
                                    exceedeances = threshold_data[,100:118],
                                    q.s=getThresholds(test, 0.95)[100:118],
                                    horizon = horizon,
                                    params = val.params[100:118,],
                                    n_samples = s.sample,
                                    name = "conditional-mat-test")


cont_mat %>% print
names(cont_mat)
cont_mat[[1]]$quantiles.values
