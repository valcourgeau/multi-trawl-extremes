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

setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
# list of packages
library(spatstat)
library(evir)
library(gPdtest)
library(DEoptim)
library(fExtremes)

# source files
source('pairwise_latent_trawl.R')

# Function
{
  create_ig_trawl <- function(dataset, tail_ind, qs, rho, kappa)
  {
    d <- length(dataset[1,]) - 1
    params <- rep(0,4*d)
    
    for(i in 1:d){
      fit <- evir::gpd(dataset[,i+1][tail_ind[i][[1]]], method = "ml", threshold = qs[i])
      params[(4*(i-1)+1):(4*(i-1)+2)] <- fit$par.ests
      params[(4*(i-1)+3):(4*(i-1)+4)] <- c(rho, kappa)
    }
    
    return(params)
  }
  
  loglik_pl <- function(times, values, deltas, logscale=T){
    d <- length(values[1,])
    return(function(params){
              return(pl_final(times = times,
                            values = values,
                            deltas = deltas,
                            params = params,
                            logscale = logscale))
            })
  }
}

# Prepare data
{
  colnames(bl) <- col_names
  bl <- subset(bl, select=c("Date", "O3", "NO2x", "SO2x", "NO", "CO"))
  valid_data_indices <- which(is.na(bl)) %% length(bl$O3)
  bl <- bl[-valid_data_indices,]
  bl <- bl[-39413,]
  SAME_QUANTILE <- 0.7
  SAME_QUANTILE_G <- 0.7
  
  bl_q <-  c(
    quantile(bl$O3, SAME_QUANTILE, na.rm = T)[[1]],
    quantile(bl$NO2x, SAME_QUANTILE, na.rm = T)[[1]],
    quantile(bl$SO2x, SAME_QUANTILE, na.rm = T)[[1]],
    quantile(bl$NO, SAME_QUANTILE, na.rm = T)[[1]],
    quantile(bl$CO, SAME_QUANTILE, na.rm = T)[[1]]
  )
  
  bl_pdf <- c(
    approxfun(density(bl$O3, na.rm = T)),
    approxfun(density(bl$NO2x, na.rm = T)),
    approxfun(density(bl$SO2x, na.rm = T)),
    approxfun(density(bl$NO, na.rm = T)),
    approxfun(density(bl$CO, na.rm = T))
  )
  
  bl_cdf <- c(
    CDF(density(bl$O3, na.rm = T)),
    CDF(density(bl$NO2x, na.rm = T)),
    CDF(density(bl$SO2x, na.rm = T)),
    CDF(density(bl$NO, na.rm = T)),
    CDF(density(bl$CO, na.rm = T))
  )
  
  bl_tail <- list(
    O3=which(bl$O3 > bl_q[1]),
    NO2=which(bl$NO2x > bl_q[2]),
    SO2=which(bl$SO2x > bl_q[3]),
    NO=which(bl$NO > bl_q[4]),
    CO=which(bl$CO > bl_q[5])
  )
  
  bl_values <- bl[,-1]
  bl_times <- matrix(rep(seq(0,1,length.out = length(bl[,1])), length(bl_values[1,])),
                     ncol=length(bl_values[1,]))
  
  # Fitting marginals
  bl_rho <- 0.3
  bl_kappa <- 10
  bl_deltas <- rep(2, 5)
  params_init <- create_ig_trawl(bl, bl_tail, bl_q, bl_rho, bl_kappa)
  
  
  # Create Gumbel marginals
  gb <- create_gumbel(dataset = bl, quantiles = bl_q, cdfs = bl_cdf, params = params_init)
  
  gb_q <-  c(
    quantile(gb$O3, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$NO2x, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$SO2x, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$NO, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$CO, SAME_QUANTILE_G, na.rm = T)[[1]]
  )
  
  gb_tail <- list(
    O3=which(gb$O3 > gb_q[1]),
    NO2=which(gb$NO2x > gb_q[2]),
    SO2=which(gb$SO2x > gb_q[3]),
    NO=which(gb$NO > gb_q[4]),
    CO=which(gb$CO > gb_q[5])
  )
}

params_init <- create_ig_trawl(bl, bl_tail, bl_q, bl_rho, bl_kappa)
fn_to_optim <- loglik_pl(times = bl_times,
                         values = bl_values,
                         deltas = bl_deltas)
fn_to_optim(params_init)
optim(par = params_init, fn = fn_to_optim)
