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
