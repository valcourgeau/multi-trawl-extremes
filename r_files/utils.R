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



#' This function replicates the vector to the length of data 
#' by vertically stacking it.
#' @param data cleaned dataset
#' @param vector_to_rep vector that will be replicated
#' @examples 
#' data <- matrix(runif(9), ncol=3)
#' vector_to_rep <- c(1,2,3)
#' makeMatrix(data = data, vector_to_rep = vector_to_rep)
makeMatrix <- function(data, vector_to_rep){
  n <- length(data[,1])
  rep_thres <- t(matrix(rep(vector_to_rep, n), ncol=n))
  
  return(rep_thres)
}

#' Creates a string from current system clock.
#' @return a string in the format "%Y-%m-%d-%H-%M-%S"
#' @example makeRdmTimestamp()
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

#' Allows user to create a filename using a main filename,
#' a tag and an extension.
#' @param file_name Main name of the file. If not given, replaced 
#'                  by a timestamp (see makeRdmTimestamp).
#' @param tag Additional tag to append if file_name is not given.
#' @param extension File extension (e.g. "RData).
#' @return A string in the format file_name + tag + extension.
#' @examples 
#' makeFileName(file_name="ok", tag="tag", extension=".RData") #returns "file.RData".
#' makeFileName(tag="tag", extension=".csv") #returns makeRdmTimestamp + "tag.csv".
makeFileName <- function(file_name=NA, tag, extension){
  if(file_name %>% is.na){
    return(paste(makeRdmTimestamp(), tag, ".RData", sep=""))
  }else{
    return(paste(file_name, ".RData", sep=""))
  }
}

makeFileName("ok", "_matrix", "RData")


#' A wrapper for the automated threshold selection from Bader et al. (2018).
#'
#' @param data dataset (vector or matrix)
#' @param p.zeroes scalar or vector (to the length of \code{data}) of
#'   probability of NOT having an extreme value (or threshold probability)
#' @return List with column names as keys and threshold selection tests results
#'   as values.
#' @seealso Bader, Brian; Yan, Jun; Zhang, Xuebin. Automated threshold selection
#'   for extreme value analysis via ordered goodness-of-fit tests with
#'   adjustment for false discovery rate. Ann. Appl. Stat. 12 (2018), no. 1,
#'   310--329. doi:10.1214/17-AOAS1092.
#'   https://projecteuclid.org/euclid.aoas/1520564474
ChoosingThresholds <- function(data, p.zeroes){
  if(!is.vector(data)){
    ncol <- length(data[1,])
  }
  if(length(p.zeroes) == 1){
    p.zeroes <- rep(p.zeroes, ncol)
  }
  col_names <- colnames(data)
  if(length(p.zeroes))
  
  tests_results <- list() 
  
  for(i in 1:ncol){
    thres_test <- threshold_test(data[,i], p.zero = p.zeroes[i])
    tests_results[[col_names[i]]] <- thres_test
  }
  
  return(tests_results)
}