setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")
col_names <- c("Date", "Time",
               "O3", "Status O3",
               "NO", "Status NO",
               "NO2", "Status NO2",
               "NO2x", "Status NO2x",
               "SO2", "Status SO2",
               "SO2x", "Status SO2x",
               "CO", "Status CO",
               "PM10", "Status PM10",
               "PM2.5", "Status PM2.5")
col_classes <- c("Date", "Date", 
                 "numeric", "character", 
                 "numeric", "character", 
                 "numeric", "character", 
                 "numeric", "character", 
                 "numeric", "character", 
                 "numeric", "character", 
                 "numeric", "character",
                 "numeric", "character")

bl <- read.csv("bloomsbury_1994_1998.csv", sep = ",",
               as.is = T)
bl <- bl[1:10000,]
# list of packages
library(spatstat)
library(evir)
library(gPdtest)
library(DEoptim)
library(fExtremes)