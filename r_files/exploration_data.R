setwd("C:/Users/Valentin/Documents/MRes/Data")
source("data_extraction.R")
library(stats)
library(zoo)

plot_acf_pacf <- function(data_to_plot, na_action=na.pass, lag.max=40)
{
  dev.off(dev.list()["RStudioGD"]) # clean the plotting device
  par(mfrow=c(2,5), mai=c(.6,0.3,0.4,0.2)) # set number of plots and margins
  acf(data_to_plot$O3, na.action = na_action, main = "ACF O3", lag.max = lag.max)
  acf(data_to_plot$NO2x, na.action = na_action, main = "ACF NO2",lag.max =  lag.max)
  acf(data_to_plot$SO2x, na.action = na_action, main = "ACF SO2x",lag.max =  lag.max)
  acf(data_to_plot$NO, na.action = na_action, main = "ACF NO",lag.max =  lag.max)
  acf(data_to_plot$PM10, na.action = na_action, main = "ACF PM10",lag.max =  lag.max)
  
  pacf(data_to_plot$O3, na.action = na_action, main = "PACF O3",lag.max =  lag.max)
  pacf(data_to_plot$NO2x, na.action = na_action, main = "PACF NO2",lag.max =  lag.max)
  pacf(data_to_plot$SO2x, na.action = na_action, main = "PACF SO2x",lag.max =  lag.max)
  pacf(data_to_plot$NO, na.action = na_action, main = "PACF NO",lag.max =  lag.max)
  pacf(data_to_plot$PM10, na.action = na_action, main = "PACF PM10", lag.max = lag.max)
}

plot_acf_pacf_exc <- function(data_to_plot, na_action=na.pass, lag.max=40, marg.threshold = rep(0.0, times=5))
{
  dev.off(dev.list()["RStudioGD"]) # clean the plotting device
  par(mfrow=c(2,5), mai=c(.6,0.3,0.4,0.2)) # set number of plots and margins
  
  O3_transformed <-    pmax(data_to_plot$O3 - marg.threshold[1], 0.0)
  NO2_transformed <-   pmax(data_to_plot$NO2x - marg.threshold[2], 0.0)
  SO2x_transformed <-  pmax(data_to_plot$SO2x - marg.threshold[3], 0.0)
  NO_transformed <-    pmax(data_to_plot$NO - marg.threshold[4], 0.0)
  PM10_transformed <-  pmax(data_to_plot$PM10 - marg.threshold[5], 0.0)
  
  acf(O3_transformed, na.action = na_action, main = "ACF O3", lag.max = lag.max)
  acf(NO2_transformed, na.action = na_action, main = "ACF NO2",lag.max =  lag.max)
  acf(SO2x_transformed, na.action = na_action, main = "ACF SO2x",lag.max =  lag.max)
  acf(NO_transformed, na.action = na_action, main = "ACF NO",lag.max =  lag.max)
  acf(PM10_transformed, na.action = na_action, main = "ACF PM10",lag.max =  lag.max)
  
  pacf(O3_transformed, na.action = na_action, main = "PACF O3",lag.max =  lag.max)
  pacf(NO2_transformed, na.action = na_action, main = "PACF NO2",lag.max =  lag.max)
  pacf(SO2x_transformed, na.action = na_action, main = "PACF SO2x",lag.max =  lag.max)
  pacf(NO_transformed, na.action = na_action, main = "PACF NO",lag.max =  lag.max)
  pacf(PM10_transformed, na.action = na_action, main = "PACF PM10", lag.max = lag.max)
}


plot_histograms <- function(data_to_plot, marg.threshold=rep(0.0, times=5), na_action=na.pass, brk_hist = 40,
                            brk_exc = 10, exc_without_zero=F)
{
  dev.off(dev.list()["RStudioGD"]) # clean the plotting device
  par(mfrow=c(2,5), mai=c(.6,0.5,0.4,0.2)) # set number of plots and margins
  breaks_each = brk_hist
  hist(data_to_plot$O3, main = "Hist O3", breaks = breaks_each, probability = T)
  hist(data_to_plot$NO2, main = "Hist NO2", breaks = breaks_each, probability = T)
  hist(data_to_plot$SO2x, main = "Hist SO2x", breaks = breaks_each, probability = T)
  hist(data_to_plot$NO, main = "Hist NO", breaks = breaks_each, probability = T)
  hist(data_to_plot$PM10, main = "Hist PM10", breaks = breaks_each, probability = T)
  
  O3_transformed <-    pmax(data_to_plot$O3 - marg.threshold[1], 0.0)
  NO2_transformed <-   pmax(data_to_plot$NO2x - marg.threshold[2], 0.0)
  SO2x_transformed <-  pmax(data_to_plot$SO2x - marg.threshold[3], 0.0)
  NO_transformed <-    pmax(data_to_plot$NO - marg.threshold[4], 0.0)
  PM10_transformed <-  pmax(data_to_plot$PM10 - marg.threshold[5], 0.0)
  
  if(exc_without_zero){
    O3_transformed <- O3_transformed[O3_transformed > 0.0]
    NO2_transformed <- NO2_transformed[NO2_transformed > 0.0]
    SO2x_transformed <- SO2x_transformed[SO2x_transformed > 0.0]
    NO_transformed <- NO_transformed[NO_transformed > 0.0]
    PM10_transformed <- PM10_transformed[PM10_transformed > 0.0]
  }
  
  breaks_each = brk_exc
  hist(O3_transformed, main = "Hist exceedances O3", breaks = breaks_each, probability = T)
  hist(NO2_transformed, main = "Hist exceedances NO2", breaks = breaks_each, probability = T)
  hist(SO2x_transformed, main = "Hist exceedances SO2x", breaks = breaks_each, probability = T)
  hist(NO_transformed, main = "Hist exceedances NO", breaks = breaks_each, probability = T)
  hist(PM10_transformed, main = "Hist exceedances PM10", breaks = breaks_each, probability = T)
  
}

empirical_mean_excess <- function(univ_data, plot=F)
{
  # univ_data: 1-d array of data
  # returns: 2 * n matrix (x_k, \hat{e}(x_k))
  univ_data <- univ_data[!is.na(univ_data)]
  n <- length(univ_data)
  results <- matrix(0.0, 2, n)
  lin.reg.order <- n
  lin.reg.sd <- Inf
  for (i in 1:n)
   {
     u <- univ_data[i]
     # u is the current threshold
     univ_data_above_threshold <- univ_data[which(univ_data > u)]
     results[1,i] <- u
     if(length(univ_data_above_threshold) == 0)
     {
       results[2,i] <- 0
     } else
     {
       results[2,i] <- sum(univ_data_above_threshold - u)/ length(univ_data_above_threshold) 
     }
  }
  
  if(plot)
   {
     # under the assumption of iid samples
     dev.off(dev.list()["RStudioGD"]) # clean the plotting device
     par(mfrow=c(1,1), mai=c(.9,0.8,0.5,0.1)) # set number of plots and margins
     x.vals <- results[1,]
     y.vals <- results[2,]
     plot(x.vals, y.vals, main = "Mean excess", 
          xlab = "Value", ylab = "Mean excess function")
    }
  
   return(results)
}

empirical_median_excess <- function(univ_data, plot=F)
{
  # univ_data: 1-d array of data
  # returns: 2 * n matrix (x_k, \hat{e}(x_k))
  univ_data <- univ_data[!is.na(univ_data)]
  n <- length(univ_data)
  results <- matrix(0.0, 2, n)
  lin.reg.order <- n
  lin.reg.sd <- Inf
  for (i in 1:n)
  {
    u <- univ_data[i]
    # u is the current threshold
    univ_data_above_threshold <- univ_data[which(univ_data > u)]
    results[1,i] <- u
    if(length(univ_data_above_threshold) == 0)
    {
      results[2,i] <- 0
    } else
    {
      results[2,i] <- median(univ_data_above_threshold-u, na.rm = T)[[1]]
    }
  }
  
  if(plot)
  {
    # under the assumption of iid samples
    dev.off(dev.list()["RStudioGD"]) # clean the plotting device
    par(mfrow=c(1,1), mai=c(.9,0.8,0.5,0.1)) # set number of plots and margins
    x.vals <- results[1,]
    y.vals <- results[2,]
    plot(x.vals, y.vals, main = "Median excess", 
         xlab = "Value", ylab = "Median excess function")
  }
  
  return(results)
}


empirical_threshold_estimate <- function(data_matrix, plot=F)
{
  # data_matrix: matrix resulting from empirical mean/median excess function
  # returns: threshold s.t. best linear approx with max point
  
  data_to_use <- data_matrix[,order(data_matrix[1,])]
  n <- length(data_to_use[,1][!is.na(data_to_use[1,])])
  x.vals <- data_to_use[1,]
  y.vals <- data_to_use[2,]
  
  lin.sd <- Inf
  lin.order <- 1
  
  for(i in 1:(n-10))
  {
    model <- lm(y.vals[i:n] ~ x.vals[i:n])
    if(lin.sd > (sd(model$residuals))/i^3)
    {
      lin.sd <- sd(model$residuals)/i^3
      lin.order <- i
    }
  }
  
  the_value <- x.vals[lin.order]
  the_model <- lm(y.vals[which(x.vals >=the_value)] ~ x.vals[which(x.vals >=the_value)])
  plot(x.vals, y.vals, main = paste("Mean excess with best-fitted straight line (", the_value, ")"))
  abline(v = x.vals[lin.order], lty = 3)
  lines(x.vals[which(x.vals >=the_value)], the_model$fitted.values, lwd = 2, lty = 5, col = "red")
  print(lin.order)
  return(the_value)
}

plot_stl <- function(univ_data, name="data", frequ=365, na.action=na.locf)
{
  univ.ts <- ts(univ_data, names = name, frequency = frequ)
  stl_test <- stl(univ.ts, na.action = na.locf, s.window = "periodic")
  plot(stl_test, main = paste("SLT for ",name))
  
  return(stl_test)
}
