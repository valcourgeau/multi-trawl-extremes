library(evir)
library(forecast)
library(lubridate)
library(imputeTS)

merged_dataset_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/merged-datasets/"
energy_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/hourly-energy-consumption/"
weather_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/historical-hourly-weather-data/"


reading_stack_csv <- function(){
  temp <- list.files(pattern="*.csv")
  myfiles <- lapply(temp, read.csv)
  names(myfiles) <- sapply(temp, function(x){sub('\\.csv$', '', x)})
  
  return(myfiles)
} 

concat_csv <- function(by.col="datetime"){
  myfiles <- reading_stack_csv()
  data <- data.frame()
  n_files <- length(myfiles)

  for(i in 1:n_files){
    names(myfiles[[i]]) <- c(names(myfiles[[1]])[1],
                             paste(names(myfiles)[i], ".", 
                                   names(myfiles[[i]])[-1], sep = ""))
    if(i == 1){
      data <- myfiles[[i]]
    }else{
      data <- merge(x=data, y=myfiles[[i]], by = by.col)
    }
  }
  return(data)
}

create_date_index <- function(data, freq=24){
  # freq = 24 means hourly data
  n <- length(data[,1])
  data['index'] <- (0:(n-1)) / freq
  return(data)
}

filter_na <- function(data, exclude.cols=NA, mode=NA){
  if(is.na(mode)){
    return(data[rowSums(is.na(data)) == 0,])
  }else{
    if(mode == "interp"){
      data[,-which(colnames(data) %in% exclude.cols)] <- 
        na.interpolation(data[,-which(colnames(data) %in% exclude.cols)])
    }else{
      stop('mode should be NA or "inter".')
    }
  }
  return(data)
}

prep_data <- function(by.col="datetime"){
  cat("concatenation...")
  ptm <- proc.time()
  data <- concat_csv(by.col = by.col)
  cat(paste("done in", round((proc.time()-ptm)[3], digits = 2), "s.\n"))
  
  cat("time indexing data...")
  ptm <- proc.time()
  data <- create_date_index(data)
  cat(paste("done in", round((proc.time()-ptm)[3], digits = 2), "s.\n"))
  
  cat("NA handling...")
  ptm <- proc.time()
  data <- filter_na(data, exclude.cols = c(by.col, "index"), mode = "interp")
  cat(paste("done in", round((proc.time()-ptm)[3], digits = 2), "s.\n"))
  cat("done\n")
  
  return(data)
}

# weather data
setwd(weather_folder)
hourly_weather_merged <- prep_data()
save(hourly_weather_merged, file=paste(merged_dataset_folder, "weather-data.Rda", sep = ""))

setwd(energy_folder)
hourly_energy_merged <- prep_data("Datetime")
save(hourly_energy_merged, file= paste(merged_dataset_folder, "energy-data.Rda", sep = ""))

names(hourly_energy_merged) <- c(names(hourly_weather_merged)[1], names(hourly_energy_merged)[-1])
names(hourly_energy_merged)
energy_weather_merged <- merge(x = hourly_weather_merged, y = hourly_energy_merged, by="datetime")
save(energy_weather_merged, file=paste(merged_dataset_folder, "energy-weather-merge.Rda", sep=""))
#load(file = "energy-weather-merge.Rda")
