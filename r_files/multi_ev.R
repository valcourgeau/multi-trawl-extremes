# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)
require(stats)
require(evir)
source("prep_univariate_latent_trawl_fit.R")

#' @examples 
#' ok <- ignoreColStartingWith(data = energy_weather_merged, tags = c("humidity",
#'                                                                  "wind_direction"))
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


#' @examples cols_to_ignore <- c("datetime", "index.x", "index.y")
#' types_to_ignore <- c("factor")
#' tags_to_ignore <- c("humidity",
#'                     "wind_direction")
#' 
#' core_energy_data <- getCoreData(data = energy_weather_merged, 
#'                                 ignore_tags = tags_to_ignore,
#'                                 ignore_cols = cols_to_ignore,
#'                                 ignore_data_type = types_to_ignore)
getCoreData <- function(data, ignore_tags=c(), ignore_cols=c(), ignore_data_type=c()){
  return(data[, ignoreColStartingWith(data, ignore_tags, return.filters = T) &
                (!colnames(data) %in% ignore_cols) & 
                !(sapply(data, class) %in% ignore_data_type)])
}



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

#' @examples 
#' test <- datasetCleaning(data = core_energy_data,
#'                         dates = energy_weather_merged$datetime)
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

#' @examples getThresholds(core_energy_data[,1:3], c(0.1,0.5,0.8))
getThresholds <- function(data, p.exceed){
  if(any(p.exceed < 0) | any(p.exceed > 1)){
    stop('p.exceed should be between 0 and 1.')
  }
  
  # deals with data as a vector
  if(is.vector(data)){
    return(quantile(data, p.exceed[1]))
  }
  
  if(length(p.exceed) == 1){
    return(sapply(data, function(x){quantile(x, p.exceed)}))
  }else{
    if(length(data[1,]) == length(p.exceed)){
      return(sapply(1:length(data[1,]), function(x){quantile(data[,x], p.exceed[x])}))
    }else{
      stop('p.exceed should either be a scalar or as wide as data.')
    }
  }
}
getThresholds(core_energy_data[,1:3], c(0.1,0.5,0.8))


makeThresholdsMatrix <- function(data, thresholds){
  n <- length(data[,1])
  rep_thres <- t(matrix(rep(thresholds, n), ncol=n))
  
  return(rep_thres)
}

#' @param data clean dataset;
#' @param thresholds d dimensional vector of threshold values;
#' @param normalize Logical. Normalise extremes to have sd = 1;
#' @examples
#' threshold_data <- makeExceedances(test, thresholds = getThresholds(test, 0.8))
makeExceedances <- function(data, thresholds, normalize=TRUE){
  if(!is.vector(thresholds)){
    stop('Thresholds should be a vector of extreme threshold value.')
  }
  if(!is.vector(data)){
    if(length(data[1,]) != length(thresholds)){
      stop('thresholds and data have non-comforting width. 
         Tip: Use getThresholds to get the right size.')
    }
    rep_thres <- makeThresholdsMatrix(data = data, thresholds = thresholds)
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

#' @param data clean dataset
#' @param q
makeConditionalMatrices <- function(data, q.s, 
                                    horizon, save=T, params, n_samples = length(data[,1]),
                                    name="conditional-matrices"){
  
  p.zeroes <- computePZero(params)
  thres <- getThresholds(data = data, p.exceed = p.zeroes)
  exceedeances <- makeExceedances(data = data, thresholds = q.s, normalize = T)
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
                         nrow = length(which(exceedeances[1:(n_samples-h), i] > 0)),
                         ncol = n_vars+1)
      # filtering the i-th eCDF column with extremes 
      temp <- exceedeances_cdf_ecdf[which(exceedeances[1:(n_samples-h), i] > 0), i]
      # addting those values in the nvars + 1 column
      mat_temp[,n_vars+1] <- ecdf(temp)(temp)
      
      for(j in 1:n_vars){
        # filtering data of j-th vars when i-th is extreme h timesteps before
        data_j <- exceedeances_cdf_ecdf[which(exceedeances[1:(n_samples-h), i] > 0)+h, j]
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

#' @examples fitExceedancesVines(threshold_data[,100:102], list_of_list_horizons)
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

