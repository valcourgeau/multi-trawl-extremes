setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/")
library("zoo")

prep.data.daily.kings <- function(rawdata, levels_l){
  pd2 <- data.frame(vapply(levels_l, function(x){rawdata$Value[rawdata$Species == x]}, 
                           FUN.VALUE = as.double(1:(floor(length(rawdata[,1])/length(levels_l))))))
  days_data <- rawdata$ReadingDateTime[rawdata$Species == levels_l[1]]
  for(i in 1:length(levels_l)){
    cat(levels_l[i], length(which(is.na(pd2[,i])))/length(pd2[,1]), "\n")
    pd2[,i] <- zoo::na.aggregate(object = zoo(pd2[,i], order.by = days_data), by=lubridate::wday)
  }
  cat("OVERALL", length(which(rowSums(is.na(pd2))>0))/length(pd2[,1]), "\n")
  rawdata <- cbind(days_data, pd2)
  cat("OVERALL", length(which(rowSums(is.na(rawdata))>0))/length(rawdata[,1]), "\n")
  
  rawdata <- cbind(1:length(rawdata[,1]), rawdata)
  rawdata[,2] <- format(as.POSIXct(rawdata[,2],format='%d/%m/%Y %H:%M'),format='%d/%m/%Y')
  colnames(rawdata) <- c("index", "date", levels_l)
  rawdata <- rawdata[,]
  return(rawdata)
}

prep.save.daily.kings <- function(rawdata, levels_l, overwrite=F, name="data.csv"){
  preppy <- prep.data.daily.kings(rawdata, levels_l)
  if(overwrite | !file.exists(name)){
    write.csv(preppy, file = name, row.names=F)
  }
  #return(preppy)
}

prep.save.daily.kings(read.csv("daily_cromwell.csv"), levels_l = c("CO", "NO", "NO2", "PM10", "SO2"),
                overwrite = T, name = "daily_cromwell_2000_2017.csv")
prep.save.daily.kings(read.csv("daily_bloomsbury.csv"), levels_l = c("O3", "CO", "NO", "NO2", "PM10", "SO2"),
                overwrite = T, name = "daily_bloomsbury_2000_2017.csv")

prep.data.hourly.kings <- function(rawdata, levels_l){
  pd2 <- data.frame(vapply(levels_l, function(x){rawdata$Value[rawdata$Species == x]}, 
                           FUN.VALUE = as.double(1:(floor(length(rawdata[,1])/length(levels_l))))))
  days_data <- rawdata$ReadingDateTime[rawdata$Species == levels_l[1]]
  for(i in 1:length(levels_l)){
    cat(levels_l[i], length(which(is.na(pd2[,i])))/length(pd2[,1]), "\n")
    pd2[,i] <- zoo::na.aggregate(object = zoo(pd2[,i], order.by = days_data), by=lubridate::wday)
  }
  cat("OVERALL", length(which(rowSums(is.na(pd2))>0))/length(pd2[,1]), "\n")
  rawdata <- cbind(days_data, pd2)
  cat("OVERALL", length(which(rowSums(is.na(rawdata))>0))/length(rawdata[,1]), "\n")
  
  rawdata <- cbind(1:length(rawdata[,1]), rawdata)
  rawdata <- cbind(1:length(rawdata[,1]), rawdata)
  rawdata[,2] <- format(as.POSIXct(days_data,format='%d/%m/%Y %H:%M'), format='%d/%m/%Y')
  rawdata[,3] <- format(as.POSIXct(days_data,format='%d/%m/%Y %H:%M'), format='%H:%M')
  colnames(rawdata) <- c("index", "date", "time", levels_l)
  rawdata <- rawdata[which(rowSums(is.na(rawdata))==0),]
  #rawdata <- rawdata[,]
  return(rawdata)
}

prep.save.hourly.kings <- function(rawdata, levels_l, overwrite=F, name="data.csv"){
  preppy <- prep.data.hourly.kings(rawdata, levels_l)
  if(overwrite | !file.exists(name)){
    write.csv(preppy, file = name, row.names=F)
  }
  #return(preppy)
}

prep.save.hourly.kings(read.csv("raw_hourly_bloomsbury_2000_2017.csv"), levels_l = c("O3", "CO", "NO", "NO2", "PM10", "SO2"),
                      overwrite = T, name = "hourly_bloomsbury_2000_2017.csv")
