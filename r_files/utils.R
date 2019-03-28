
#' @param data dataset
#' @param tags Column tags (e.g. partial name) to remove.
#' @param return.filters Logical flag (default FALSE). Whether to return the
#'   filtered indices (TRUE) or filtered data (FALSE).
#' @examples
#' ok <- ignoreColStartingWith(data = energy_weather_merged, tags = c("humidity",
#'                                                                  "wind_direction"))
ignoreColStartingWith <- function(data, tags, return.filters=F){
  # return.filters 
  filtered_by_tag <- sapply(tags, function(x){grepl(x, colnames(data))}) # find the right columns
  filtered_by_tag <- apply(filtered_by_tag, 
                           MARGIN = 1, 
                           FUN = function(x){Reduce('|',x)}) # merging the filters
  if(return.filters){
    return(!filtered_by_tag)
  }else{
    return(data[,!filtered_by_tag]) 
  }
}

#' Returns a filter for columns of \code{data} which contains \code{tags}.
#' @param data dataset
#' @param tags Column tags (e.g. partial name) to remove.
#' @param return.filters Logical flag (default FALSE). Whether to return the
#'   filtered indices (TRUE) or filtered data (FALSE).
#' @examples
#' ok <- GetColStartingWith(data = energy_weather_merged, tags = c("humidity",
#'                                                                  "wind_direction"))
GetColStartingWith <-  function(data, tags, return.filters=F){
  filtered <- !ignoreColStartingWith(data = data, tags = tags,
                                     return.filters = T)
  if(return.filters){
    return(filtered)
  }else{
    return(data[, filtered])
  }
}


#'
#' @param data dataset
#' @param ignore_tags Vector of colname tags to remove/ignore.
#' @param get_tags Vector of colname tags (incomplete name) to keep.
#' @param ignore_cols Vector of EXACT colnames to remove/ignore.
#' @param ignore_data_type Vector of datatypes to remove/ignore, e.g. factor or
#'   numeric.
#' @note Ignoring is preferred over keeping columns, conflicts between
#'   \code{get_tags} and \core{ignore_tags} should be checked accordingly.
#' @examples cols_to_ignore <- c("datetime", "index.x", "index.y")
#' types_to_ignore <- c("factor")
#' tags_to_ignore <- c("humidity",
#'                     "wind_direction")
#'
#' core_energy_data <- getCoreData(data = energy_weather_merged,
#'                                 ignore_tags = tags_to_ignore,
#'                                 ignore_cols = cols_to_ignore,
#'                                 ignore_data_type = types_to_ignore)
getCoreData <- function(data, ignore_tags=c(), get_tags=c(), ignore_cols=c(), ignore_data_type=c()){
  flag <- rep(T, length(colnames(data)))
  if(length(ignore_tags) > 0){
    flag <- ignoreColStartingWith(data, ignore_tags, return.filters = T)
  }
  if(length(get_tags) > 0){
    flag <- flag & GetColStartingWith(data, get_tags, return.filters = T)
  }
  return(data[, flag &
                (!colnames(data) %in% ignore_cols) & 
                !(sapply(data, class) %in% ignore_data_type)])
}

# cols_to_ignore <- c("datetime", "index.x", "index.y")
# types_to_ignore <- c("factor")
# tags_to_ignore <- c("humidity",
#                     "wind_direction")
# 
# core_energy_data <- getCoreData(data = energy_weather_merged,
#                                 ignore_tags = c(),
#                                 get_tags = tags_to_ignore,
#                                 ignore_cols = cols_to_ignore,
#                                 ignore_data_type = types_to_ignore)

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

DeterministicCleaningSequential <- function(dates, data, 
                                     wday = F, month = F,
                                     quarter = F,
                                     hour = F, seasons = NA,
                                     p.adjust.method="bonferroni"){
  #dates <- strptime(energy_weather_merged$datetime, "%Y-%m-%d %H:%M:%S", tz = "GMT")
  if(length(dates) != length(data)){
    stop('Data and dates should have same length.')
  }
  
  fitting_mat <- as.matrix(rep(1, length(data))) # creating skeleton of matrix
  
  if(!any(is.na(seasons)) & is.vector(seasons)){
    temp <- do.call(cbind,
                    lapply(seasons, 
                           FUN = function(x){
                             cbind(cos(2*pi*(1:length(dates))/x),
                                   sin(2*pi*(1:length(dates))/x))
                           }))
    fitting_mat <- cbind(fitting_mat, temp)
  }
  
  if(quarter){
    temp <- do.call(cbind,
                    lapply(1:3, 
                           FUN = function(i){quarter(dates) == i}))
    fitting_mat <- cbind(fitting_mat, temp)
  }
  if(wday){
    temp <- do.call(cbind,
                    lapply(1:6, 
                           FUN = function(i){wday(dates) == i}))
    fitting_mat <- cbind(fitting_mat, temp)
  }
  if(month){
    temp <- do.call(cbind,
                    lapply(1:12, 
                           FUN = function(i){month(dates) == i}))
    fitting_mat <- cbind(fitting_mat, temp)
  }
  if(hour){
    temp <- do.call(cbind,
                    lapply(1:23, 
                           FUN = function(i){hour(dates) == i}))
    fitting_mat <- cbind(fitting_mat, temp)
  }
  
  fit <- lm(data ~ fitting_mat) # we use fitting_mat[,-1] to remove skeleton
  
  # adjusting the method to take in account multiple p.values testing
  p.v.adjusted <- summary(fit)$coefficients[,4]
  p.v.adjusted <- p.adjust(p.v.adjusted, method = p.adjust.method)
  # TODO test that!
  # print(p.v.adjusted)
  
  fitting_indices <- which(p.v.adjusted < 0.05)
  print(length(fitting_indices))
  # if(1 %in% fitting_indices){
  #   fitting_indices <- fitting_indices[-1]
  # }
  #fitting_mat <- fitting_mat[,fitting_indices-1]
  
  fitting_mat <- fitting_mat[,fitting_indices]
  return(lm(data ~ fitting_mat)$residuals)
}


#' @param data cleaned dataset
#' @param frequency frequency of sampling of underlying timeseries per unit of
#'   time (hourly data, frequency = 24).
#' @param trend_window window to look for trends in STL algorithm.
#' @param season_window  window to look for seasonality in STL algorithm.
#' @seealso \code{\link[stats]{stl}}.
#' @examples
#' datasetSTLCleaning(data = core_energy_data,
#'                           frequency=24,
#'                           trend_window=24*365/4, # quartly trend changes
#'                           season_window=24) # daily seasonality
datasetSTLCleaning <- function(data, frequency, trend_window, season_window){
  result <- apply(data, MARGIN=2,
                  FUN = function(x){
                    return(stl(ts(x, frequency = frequency), 
                               t.window = trend_window, 
                               s.window = season_window)$time.series[,3] %>% as.vector)
                  })
  result <- data.frame(result)
  colnames(result) <- colnames(data)
  return(result)
}

#' @param data cleaned dataset
#' @param dates vector of numerical dates (timestamps)
#' @param frequenc frequency of sampling of underlying timeseries per unit of
#'   time (hourly data, frequency = 24).
#' @param 
#' @param sequential Logical flag (default = TRUE). Whether to use the
#'   sequential creation of fitting matrix with wday, month, quarter, hour,
#'   seasons of 12, 24, 48 and 24*365.
#' @examples
#' test <- datasetCleaning(data = core_energy_data,
#'                         dates = energy_weather_merged$datetime)
datasetCleaning <- function(data, dates, sequential=T){
  
  
  #we apply linear modelling cleaning on all columns
  result <- apply(data, MARGIN = 2, 
                  FUN = function(x){
                    if(sequential){
                      DeterministicCleaningSequential(data = x, dates = dates,
                                                      wday = T, month = T, quarter = T,
                                                      hour = T, seasons = c(12, 24, 48, 24*365))
                    }else{
                      deterministicCleaning(dates = dates, data = data)
                    }
                  }
  )

  result <- data.frame(result)
  colnames(result) <- colnames(data)
  return(result)
}

#' Returns the projection of \code{wind} using \code{direction} (in degrees) and
#' potentially rotated using \code{offset}.
#'
#' @param wind vector of wind data.
#' @param direction vector of wind direction (in degrees from 0 to 360).
#' @param offset offset angle to rotate
#' @return a matrix with two columns with \code{cos} \code{sin} projections
#'   respectively.
#' @examples
#' WindSpeedProjection(rnorm(mean=100, n=100), direction=rnorm(n=100)*10+150, name="wind_speed")
WindSpeedProjection <- function(wind, direction, name=NA, offset=0){
  wind <- abs(wind) # taking absolute value just in case
  direction <- ((direction+offset) %% 360)/360
  result <- cbind(wind * cos(direction), wind * sin(direction))
  result <- data.frame(result)
  if(!is.na(name)){
    colnames(result) <- c(paste(name, ".cos", sep = ""),
                          paste(name, ".sin", sep = ""))
  }
  return(result)
}


#' @examples 
#' ExtractProcessWindData(energy_weather_merged, "New.York") %>% head
ExtractProcessWindData <- function(data, tag){
  speed_filter <- GetColStartingWith(data = data, 
                                     tags = "wind_speed.", 
                                     return.filters = T)

  direction_filter <- GetColStartingWith(data=data, 
                                         tags = "wind_direction.", 
                                         return.filters = T)

  data_tagged <- GetColStartingWith(data = data, 
                                    tags = tag, 
                                    return.filters = F)
  cols_tagged <- GetColStartingWith(data = data, 
                                    tags = tag, 
                                    return.filters = T)
  direction_filter <- direction_filter & cols_tagged
  speed_filter <- speed_filter & cols_tagged
  
  data_direction <- data[, direction_filter] %>% as.vector
  data_speed <- data[, speed_filter] %>% as.vector
  name_col <- colnames(data_tagged)[GetColStartingWith(data = data_tagged, tags = "wind_speed", return.filters = T)]
  
  return(WindSpeedProjection(wind = data_speed, direction = data_direction, name = name_col))
}
#test 
# sqrt(rowSums(asd(energy_weather_merged, "New.York")^2))
# energy_weather_merged$wind_speed.New.York


#' @examples 
#' ProcessWindProjections(data = energy_weather_merged, c("New.York", "Detroit")) %>% head
ProcessWindProjections <- function(data, tags){
  return(
    do.call(cbind, 
            lapply(tags, function(x){
                                ExtractProcessWindData(data = data,
                                                       tag = x)})
            )
        )
}

#' @examples
#' ConcatAndReplaceWind(energy_weather_merged, c("New.York", "Houston")) %>% colnames
ConcatAndReplaceWind <- function(data, tags){
  concat_vals <- ProcessWindProjections(data = data, tags = tags) %>% abs
  removing_tags_index <- ignoreColStartingWith(data, c("wind_direction", "wind_speed"), return.filters = T)
  #print(removing_tags_index)
  concat_vals <- cbind(data[,removing_tags_index], concat_vals)
  return(concat_vals)
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

#' Returns the list of cities in the column names.
#' @param column_names collection of strings
#' @return List of cities
GetCityNames <- function(data){
  sol <- colnames(data)[!grepl(pattern = "_", colnames(data))]
  return(unique(sub('.*\\.', '', sol)))
}

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
    return(paste(makeRdmTimestamp(), tag, extension, sep=""))
  }else{
    return(paste(file_name, tag, extension, sep=""))
  }
}

makeFileName("ok", "_matrix", ".RData")


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


#' A function that returns a first guess for cluster size based on monotonicity of the autocorrelation function and
#' @param data dataset (vector or matrix)
#' @param p.zeroes scalar or vector (to the length of \code{data}) of
#'   probability of NOT having an extreme value (or threshold probability)
#' @return vector of cluster size respectively for each column.
ChoosingClusters <- function(data, p.zeroes){
  lag_max <- min(150, length(data[,1])/100)
  acf_min <- 0.2
  exc_data <- makeExceedances(data = data, 
                              thresholds = getThresholds(data,p.zeroes))
  clusters_tmp <- apply(exc_data, MARGIN = 2,
                    FUN = function(x){
                      acf_data <- acf(x, lag.max = lag_max, plot = F)$acf
                      first_increase <-  which.max(diff(acf_data) > 0)
                      if(first_increase == 1){
                        first_increase <- lag_max
                      }
                      return(min(first_increase, which.max(acf_data < acf_min)))
                    })
  return(clusters_tmp)
}

