setwd("C:/Users/Valentin/Documents/MRes/Data")

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

data_bloomsbury <- read.csv("bloomsbury_1994_1998.csv", sep = ",",
                            as.is = T)
colnames(data_bloomsbury) <- col_names

hat_pdf <- function(x, xi, beta, q, pdf_f, ecdf_f){
  x <- as.numeric(x)

    if(x > q){
      return((1.0 - ecdf_f(q)) * dgpd(x, xi = xi, beta = beta)[1])
    } else{
      return(pdf_f(x))
    }
}

vec_hat_pdf <- function(x_vec, xi, beta, q, pdf_f, ecdf_f, without_zero=F){
  data_to_use <- x_vec[which(x_vec > q)]
  res <- rep(0, length(data_to_use))
  for(i in 1:length(res)){
    res[i] <- hat_pdf(data_to_use[i], xi, beta, q, pdf_f, ecdf_f)
  }
  if(without_zero){
    return(res[res>0])  
  } else {
    return(res)
  }
}

vec_full_hat_pdf <- function(x_vec, xi, beta, q, pdf_f, ecdf_f, without_zero=F){
  res <- rep(0, length(x_vec))
  for(i in 1:length(res)){
    if(is.na(x_vec[i])){
      res[i] <- NA
    } else {
      res[i] <- hat_pdf(x_vec[i], xi, beta, q, pdf_f, ecdf_f)
    }
  }
  if(without_zero){
    return(res[res>0])  
  } else {
    return(res)
  }
}

library(lubridate)
month_data <- month(as.POSIXlt(data_bloomsbury$Date, format="%Y-%m-%d"))
winter_points <- which(month_data >= 2 & month_data >= 11)
early_summer_points <- which(month_data >= 4 & month_data <= 7)

bloom_winter <- data_bloomsbury[winter_points,]
bloom_summer <- data_bloomsbury[early_summer_points,]
library('fExtremes')
par(mfrow=c(1,1), mai=c(.9,0.8,0.3,0.1))
mePlot(bloom_winter$NO2x)
plot(bloom_winter$O3, bloom_winter$NO2x)
plot(bloom_summer$O3, bloom_summer$NO2x)

q_winter_O3 <- quantile(bloom_winter$O3, 0.95, na.rm = T)[[1]]
q_winter_NO2 <- quantile(bloom_winter$NO2x, 0.95, na.rm = T)[[1]]
q_winter_SO2 <- quantile(bloom_winter$SO2x, 0.95, na.rm = T)[[1]]
q_winter_NO <- quantile(bloom_winter$NO, 0.95, na.rm = T)[[1]]
q_winter_CO <- quantile(bloom_winter$CO, 0.95, na.rm = T)[[1]]

ecdf_winter_O3 <- ecdf(bloom_winter$O3)
ecdf_winter_NO2 <- ecdf(bloom_winter$NO2x)
ecdf_winter_SO2 <- ecdf(bloom_winter$SO2x)
ecdf_winter_NO <- ecdf(bloom_winter$NO)
ecdf_winter_CO <- ecdf(bloom_winter$CO)

pdf_winter_O3 <- approxfun(density(bloom_winter$O3, na.rm = T))
pdf_winter_NO2 <- approxfun(density(bloom_winter$NO2x, na.rm = T))
pdf_winter_SO2 <- approxfun(density(bloom_winter$SO2x, na.rm = T))
pdf_winter_NO <- approxfun(density(bloom_winter$NO, na.rm = T))
pdf_winter_CO <- approxfun(density(bloom_winter$CO, na.rm = T))

# params 2 * 5 matrix [xi, beta; xi, beta; xi, beta; xi, beta; xi, beta] 
# ORDER: O3, NO2, SO2, NO, CO

qs_winter <- c(q_winter_O3, q_winter_NO2, q_winter_SO2, q_winter_NO, q_winter_CO)
cdf_winter <- c(ecdf_winter_O3, ecdf_winter_NO2, ecdf_winter_SO2, ecdf_winter_NO, ecdf_winter_CO)
pdf_winter <- c(pdf_winter_O3, pdf_winter_NO2, pdf_winter_SO2, pdf_winter_NO, pdf_winter_CO)

library('gPdtest')
xi_beta_gpd <- function(data, quantile)
{
  gpd.fit(data[data > q])
}


params <- rep(0, 10)

### Initial guess
# O3
mePlot(bloom_winter$O3[bloom_winter$O3 > q_winter_O3 & !is.na(bloom_winter$O3)])
gpd.test(bloom_winter$O3[bloom_winter$O3 > q_winter_O3 & !is.na(bloom_winter$O3)])
fit <- gpd.fit(bloom_winter$O3[bloom_winter$O3 > q_winter_O3 & !is.na(bloom_winter$O3)], method = "amle")
params[1:2] <- as.numeric(fit)

# NO2
mePlot(bloom_winter$NO2x[bloom_winter$NO2x > q_winter_NO2 & !is.na(bloom_winter$NO2x)])
gpd.test(bloom_winter$NO2x[bloom_winter$NO2x > q_winter_NO2 & !is.na(bloom_winter$NO2x)])
fit <- gpd.fit(bloom_winter$NO2x[bloom_winter$NO2x > q_winter_NO2 & !is.na(bloom_winter$NO2x)], method = "amle")
params[3:4] <- as.numeric(fit)

# SO2
mePlot(bloom_winter$SO2x[bloom_winter$SO2x > q_winter_SO2 & !is.na(bloom_winter$SO2x)])
gpd.test(bloom_winter$SO2x[bloom_winter$SO2x > q_winter_SO2 & !is.na(bloom_winter$SO2x)])
fit <- gpd.fit(bloom_winter$SO2x[bloom_winter$SO2x > q_winter_SO2 & !is.na(bloom_winter$SO2x)], method = "amle")
params[5:6] <- as.numeric(fit)

# NO
mePlot(bloom_winter$NO[bloom_winter$NO > q_winter_NO & !is.na(bloom_winter$NO)])
gpd.test(bloom_winter$NO[bloom_winter$NO > q_winter_NO & !is.na(bloom_winter$NO)])
fit <- gpd.fit(bloom_winter$NO[bloom_winter$NO > q_winter_NO & !is.na(bloom_winter$NO)], method = "amle")
params[7:8] <- as.numeric(fit)

# CO
mePlot(bloom_winter$CO[bloom_winter$CO > q_winter_CO & !is.na(bloom_winter$CO)])
gpd.test(bloom_winter$CO[bloom_winter$CO > q_winter_CO & !is.na(bloom_winter$CO)])
fit <- gpd.fit(bloom_winter$CO[bloom_winter$CO > q_winter_CO & !is.na(bloom_winter$CO)], method = "amle")
params[9:10] <- as.numeric(fit)





full_hat_pdf <- function(params, qs, pdf_fs, ecdf_fs)
{
  # O3
  temp <- vec_hat_pdf(bloom_winter$O3, params[1], params[2], 
                      qs[1], pdf_fs[1][[1]], ecdf_fs[1][[1]])
  full_logLik <- sum(log(temp[!is.na(temp)]))
  
  #NO2x
  temp <- vec_hat_pdf(bloom_winter$NO2x, params[3], params[4], 
                      qs[2], pdf_fs[2][[1]], ecdf_fs[2][[1]])
  full_logLik <- full_logLik + sum(log(temp[!is.na(temp)]))
  
  # SO2x
  temp <- vec_hat_pdf(bloom_winter$SO2x, params[5], params[6], 
                      qs[3], pdf_fs[3][[1]], ecdf_fs[3][[1]])
  full_logLik <- full_logLik + sum(log(temp[!is.na(temp)]))
  
  # NO
  temp <- vec_hat_pdf(bloom_winter$NO, params[7], params[8], 
                      qs[4], pdf_fs[4][[1]], ecdf_fs[4][[1]])
  full_logLik <- full_logLik + sum(log(temp[!is.na(temp)]))
  
  # CO
  temp <- vec_hat_pdf(bloom_winter$CO, params[9], params[10], 
                      qs[5], pdf_fs[5][[1]], ecdf_fs[5][[1]])
  full_logLik <- full_logLik + sum(log(temp[!is.na(temp)]))
  
  return(full_logLik)
}

gradient_hat_pdf <- function(params, qs, pdf_fs, ecdf_fs)
{
  full_gradient <- rep(0,10)
  # O3
  temp <- vec_hat_pdf(bloom_winter$O3, params[1], params[2], 
                      qs[1], pdf_fs[1][[1]], ecdf_fs[1][[1]])
  
  
  full_gradient[2] <- (-1 + (params[1]* (params[1]+1)*(bloom_winter$O3 - qs[1]))/params[1]) / params[1]
  
  #NO2x
  temp <- vec_hat_pdf(bloom_winter$NO2x, params[3], params[4], 
                      qs[2], pdf_fs[2][[1]], ecdf_fs[2][[1]])
  full_logLik <- full_logLik + sum(1/(temp[!is.na(temp)]))
  
  # SO2x
  temp <- vec_hat_pdf(bloom_winter$SO2x, params[5], params[6], 
                      qs[3], pdf_fs[3][[1]], ecdf_fs[3][[1]])
  full_logLik <- full_logLik + sum(1/(temp[!is.na(temp)]))
  
  # NO
  temp <- vec_hat_pdf(bloom_winter$NO, params[7], params[8], 
                      qs[4], pdf_fs[4][[1]], ecdf_fs[4][[1]])
  full_logLik <- full_logLik + sum(1/(temp[!is.na(temp)]))
  
  # CO
  temp <- vec_hat_pdf(bloom_winter$CO, params[9], params[10], 
                      qs[5], pdf_fs[5][[1]], ecdf_fs[5][[1]])
  full_logLik <- full_logLik + sum(1/(temp[!is.na(temp)]))
  
  return(full_logLik)
}

fct_full_hat_pdf <- function(params){
  return(full_hat_pdf(params, qs_winter, pdf_fs = pdf_winter, ecdf_fs = ecdf_winter))
}

gr_full_hat_pdf <- function(params){
  return(gradient_hat_pdf(params, qs_winter, pdf_fs = pdf_winter, ecdf_fs = ecdf_winter))
}


? optim
library('optimr')
library('evir')
opm(params, fct_full_hat_pdf, method = "BFGS")

params
qqparetoPlot(bloom_winter$O3[bloom_winter$O3 > q_winter_O3 & !is.na(bloom_winter$O3)])
hist(rgpd(1000,xi=params[1], beta= params[2]),breaks = 50)
hist(bloom_winter$O3[bloom_winter$O3 > q_winter_O3 & !is.na(bloom_winter$O3)]-40, breaks = 20)

hist(rgpd(1000,xi=params[9], beta= params[10]),breaks = 50)
hist(bloom_winter$CO[bloom_winter$CO > q_winter_CO & !is.na(bloom_winter$CO)], breaks = 20)

bloom_winter$O3.gumbel <- -log(-log(vec_full_hat_pdf(bloom_winter$O3,
                                                xi = params[1],
                                                beta = params[2],
                                                q = q_winter_O3,
                                                pdf_f = pdf_winter_O3,
                                                ecdf_f = ecdf_winter_O3)))
bloom_winter$NO2x.gumbel <- -log(-log(vec_full_hat_pdf(bloom_winter$NO2x,
                                                     xi = params[3],
                                                     beta = params[4],
                                                     q = q_winter_NO2,
                                                     pdf_f = pdf_winter_NO2,
                                                     ecdf_f = ecdf_winter_NO2)))
bloom_winter$SO2x.gumbel <- -log(-log(vec_full_hat_pdf(bloom_winter$SO2x,
                                                     xi = params[5],
                                                     beta = params[6],
                                                     q = q_winter_SO2,
                                                     pdf_f = pdf_winter_SO2,
                                                     ecdf_f = ecdf_winter_SO2)))
bloom_winter$NO.gumbel <- -log(-log(vec_full_hat_pdf(bloom_winter$NO,
                                                     xi = params[7],
                                                     beta = params[8],
                                                     q = q_winter_NO,
                                                     pdf_f = pdf_winter_NO,
                                                     ecdf_f = ecdf_winter_NO)))
bloom_winter$CO.gumbel <- -log(-log(vec_full_hat_pdf(bloom_winter$CO,
                                                     xi = params[9],
                                                     beta = params[10],
                                                     q = q_winter_CO,
                                                     pdf_f = pdf_winter_CO,
                                                     ecdf_f = ecdf_winter_CO)))

q_gumbel_winter_O3 <- quantile(bloom_winter$O3.gumbel, 0.95, na.rm = T)[[1]]
q_gumbel_winter_NO2 <- quantile(bloom_winter$NO2x.gumbel, 0.95, na.rm = T)[[1]]
q_gumbel_winter_SO2 <- quantile(bloom_winter$SO2x.gumbel, 0.95, na.rm = T)[[1]]
q_gumbel_winter_NO <- quantile(bloom_winter$NO.gumbel, 0.95, na.rm = T)[[1]]
q_gumbel_winter_CO <- quantile(bloom_winter$CO.gumbel, 0.95, na.rm = T)[[1]]

q_gumbel_full <- c(q_gumbel_winter_O3,
                   q_gumbel_winter_NO2,
                   q_gumbel_winter_SO2,
                   q_gumbel_winter_NO,
                   q_gumbel_winter_CO)

ecdf_gumbel_winter_O3 <- ecdf(bloom_winter$O3.gumbel)
ecdf_gumbel_winter_NO2 <- ecdf(bloom_winter$NO2x.gumbel)
ecdf_gumbel_winter_SO2 <- ecdf(bloom_winter$SO2x.gumbel)
ecdf_gumbel_winter_NO <- ecdf(bloom_winter$NO.gumbel)
ecdf_gumbel_winter_CO <- ecdf(bloom_winter$CO.gumbel)

ecdfs_winter_gumbel <- c(ecdf_gumbel_winter_O3,
                  ecdf_gumbel_winter_NO2,
                  ecdf_gumbel_winter_SO2,
                  ecdf_gumbel_winter_NO,
                  ecdf_gumbel_winter_CO)

pdf_gumbel_winter_O3 <- approxfun(density(bloom_winter$O3.gumbel, na.rm = T))
pdf_gumbel_winter_NO2 <- approxfun(density(bloom_winter$NO2x.gumbel, na.rm = T))
pdf_gumbel_winter_SO2 <- approxfun(density(bloom_winter$SO2x.gumbel, na.rm = T))
pdf_gumbel_winter_NO <- approxfun(density(bloom_winter$NO.gumbel, na.rm = T))
pdf_gumbel_winter_CO <- approxfun(density(bloom_winter$CO.gumbel, na.rm = T))

pdfs_gumbel <- c(pdf_gumbel_winter_O3,
                 pdf_gumbel_winter_NO2,
                 pdf_gumbel_winter_SO2,
                 pdf_gumbel_winter_NO,
                 pdf_gumbel_winter_CO)

make_indices <- function(dataset, quantiles, ecdf_fs){
  C_indexes <- which(dataset$O3.gumbel > quantiles[1] 
                     | dataset$NO2x.gumbel > quantiles[2]
                     | dataset$SO2x.gumbel > quantiles[3]
                     | dataset$NO.gumbel > quantiles[4]
                     | dataset$CO > quantiles[5])
  
  C1_indexes <- C_indexes[which(ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes]) 
                                > ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes]) 
                                & 
                                  ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes])
                                > ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes])
                                & 
                                  ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes])
                                > ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                &
                                  ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes])
                                > ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes]))]
  
  
  C2_indexes <- C_indexes[which(ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes]) 
                                > ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes]) 
                                & 
                                  ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes]) 
                                > ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes])
                                & 
                                  ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes]) 
                                > ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                &
                                  ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes]) 
                                > ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes]))]
  
  C3_indexes <- C_indexes[which(ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes]) 
                                > ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes]) 
                                & 
                                  ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes]) 
                                > ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes])
                                & 
                                  ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes]) 
                                > ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                &
                                  ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes]) 
                                > ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes]))]
  
  C4_indexes <- C_indexes[which(ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                > ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes]) 
                                & 
                                  ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                > ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes])
                                & 
                                  ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                > ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes]) 
                                &
                                  ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes])
                                > ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes]))]
  
  C5_indexes <- C_indexes[which(ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes])
                                > ecdf_fs[1][[1]](dataset$O3.gumbel[C_indexes]) 
                                & 
                                  ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes])
                                > ecdf_fs[2][[1]](dataset$NO2x.gumbel[C_indexes])
                                & 
                                  ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes])
                                > ecdf_fs[3][[1]](dataset$SO2x.gumbel[C_indexes]) 
                                &
                                  ecdf_fs[5][[1]](dataset$CO.gumbel[C_indexes])
                                > ecdf_fs[4][[1]](dataset$NO.gumbel[C_indexes]))]
  return(list("C1"=C1_indexes,"C2"=C2_indexes,"C3"=C3_indexes,"C4"=C4_indexes,"C5"=C5_indexes))
}


indices_full <- make_indices(bloom_winter, q_gumbel_full, ecdf_fs = ecdfs_gumbel)
indices_full["C1"]


single_cond_log_lik <- function(params, dataset, d, vec_indices)
{
  temp <- 0
  n <- length(dataset$O3.gumbel)
  a <- params[1:d]
  b <- params[(d+1):(2*d)]
  mu <-  params[(2*d+1):(3*d)]
  sigma.o <- params[(3*d+1):(4*d)]
  
  
  # O3
  mu_vec_1 <- a[2] * dataset$NO2x.gumbel[vec_indices[1+1]] + mu[2] * (dataset$NO2x.gumbel[vec_indices[1+1]])^b[2]
  mu_vec_2 <- a[3] * dataset$SO2x.gumbel[vec_indices[1+1]] + mu[3] * (dataset$SO2x.gumbel[vec_indices[1+1]])^b[3]
  mu_vec_3 <- a[4] * dataset$NO.gumbel[vec_indices[1+1]] + mu[4] * (dataset$NO.gumbel[vec_indices[1+1]])^b[4]
  mu_vec_4 <- a[5] * dataset$CO.gumbel[vec_indices[1+1]] + mu[5] * (dataset$CO.gumbel[vec_indices[1+1]])^b[5]
  sigma_vec_1 <- sigma.o[2] * (dataset$NO2x.gumbel[vec_indices[1+1]])^b[2]
  sigma_vec_2 <- sigma.o[3] * (dataset$SO2x.gumbel[vec_indices[1+1]])^b[3]
  sigma_vec_3 <- sigma.o[4] * (dataset$NO.gumbel[vec_indices[1+1]])^b[4]
  sigma_vec_4 <- sigma.o[5] * (dataset$CO.gumbel[vec_indices[1+1]])^b[5]
  
  temp <- temp-sum(log(sigma_vec_1) + 0.5 * ((dataset$NO2x.gumbel[vec_indices[1+2]] - mu_vec_1)/sigma_vec_1)^2)
  temp <- temp-sum(log(sigma_vec_2) + 0.5 * ((dataset$SO2x.gumbel[vec_indices[1+2]] - mu_vec_2)/sigma_vec_2)^2)
  temp <- temp-sum(log(sigma_vec_3) + 0.5 * ((dataset$NO.gumbel[vec_indices[1+2]] - mu_vec_3)/sigma_vec_3)^2)
  temp <- temp-sum(log(sigma_vec_4) + 0.5 * ((dataset$CO.gumbel[vec_indices[1+2]] - mu_vec_4)/sigma_vec_4)^2)
  
  
  # NO2x
  mu_vec_1 <- a[1] * dataset$O3.gumbel[vec_indices[1+2]] + mu[1] * (dataset$O3.gumbel[vec_indices[1+2]])^b[1]
  mu_vec_2 <- a[3] * dataset$SO2x.gumbel[vec_indices[1+2]] + mu[3] * (dataset$SO2x.gumbel[vec_indices[1+2]])^b[3]
  mu_vec_3 <- a[4] * dataset$NO.gumbel[vec_indices[1+2]] + mu[4] * (dataset$NO.gumbel[vec_indices[1+2]])^b[4]
  mu_vec_4 <- a[5] * dataset$CO.gumbel[vec_indices[1+2]] + mu[5] * (dataset$CO.gumbel[vec_indices[1+2]])^b[5]
  sigma_vec_1 <- sigma.o[1] * (dataset$O3.gumbel[vec_indices[1+2]])^b[1]
  sigma_vec_2 <- sigma.o[3] * (dataset$SO2x.gumbel[vec_indices[1+2]])^b[3]
  sigma_vec_3 <- sigma.o[4] * (dataset$NO.gumbel[vec_indices[1+2]])^b[4]
  sigma_vec_4 <- sigma.o[5] * (dataset$CO.gumbel[vec_indices[1+2]])^b[5]
  
  temp <- temp-sum(log(sigma_vec_1) + 0.5 * ((dataset$O3.gumbel[vec_indices[1+2]] - mu_vec_1)/sigma_vec_1)^2)
  temp <- temp-sum(log(sigma_vec_2) + 0.5 * ((dataset$SO2x.gumbel[vec_indices[1+2]] - mu_vec_2)/sigma_vec_2)^2)
  temp <- temp-sum(log(sigma_vec_3) + 0.5 * ((dataset$NO.gumbel[vec_indices[1+2]] - mu_vec_3)/sigma_vec_3)^2)
  temp <- temp-sum(log(sigma_vec_4) + 0.5 * ((dataset$CO.gumbel[vec_indices[1+2]] - mu_vec_4)/sigma_vec_4)^2)
  
  # SO2x
  mu_vec_1 <- a[1] * dataset$O3.gumbel[vec_indices[1+3]] + mu[1] * (dataset$O3.gumbel[vec_indices[1+3]])^b[1]
  mu_vec_2 <- a[2] * dataset$NO2x.gumbel[vec_indices[1+3]] + mu[2] * (dataset$NO2x.gumbel[vec_indices[1+3]])^b[2]
  mu_vec_3 <- a[4] * dataset$NO.gumbel[vec_indices[1+3]] + mu[4] * (dataset$NO.gumbel[vec_indices[1+3]])^b[4]
  mu_vec_4 <- a[5] * dataset$CO.gumbel[vec_indices[1+3]] + mu[5] * (dataset$CO.gumbel[vec_indices[1+3]])^b[5]
  sigma_vec_1 <- sigma.o[1] * (dataset$O3.gumbel[vec_indices[1+3]])^b[1]
  sigma_vec_2 <- sigma.o[2] * (dataset$NO2x.gumbel[vec_indices[1+3]])^b[2]
  sigma_vec_3 <- sigma.o[4] * (dataset$NO.gumbel[vec_indices[1+3]])^b[4]
  sigma_vec_4 <- sigma.o[5] * (dataset$CO.gumbel[vec_indices[1+3]])^b[5]
  
  temp <- temp-sum(log(sigma_vec_1) + 0.5 * ((dataset$O3.gumbel[vec_indices[1+3]] - mu_vec_1)/sigma_vec_1)^2)
  temp <- temp-sum(log(sigma_vec_2) + 0.5 * ((dataset$NO2x.gumbel[vec_indices[1+3]] - mu_vec_2)/sigma_vec_2)^2)
  temp <- temp-sum(log(sigma_vec_3) + 0.5 * ((dataset$NO.gumbel[vec_indices[1+3]] - mu_vec_3)/sigma_vec_3)^2)
  temp <- temp-sum(log(sigma_vec_4) + 0.5 * ((dataset$CO.gumbel[vec_indices[1+3]] - mu_vec_4)/sigma_vec_4)^2)
  
  
  
  # NO
  mu_vec_1 <- a[1] * dataset$O3.gumbel[vec_indices[1+4]] + mu[1] * (dataset$O3.gumbel[vec_indices[1+4]])^b[1]
  mu_vec_2 <- a[2] * dataset$NO2x.gumbel[vec_indices[1+4]] + mu[2] * (dataset$NO2x.gumbel[vec_indices[1+4]])^b[2]
  mu_vec_3 <- a[3] * dataset$SO2x.gumbel[vec_indices[1+4]] + mu[3] * (dataset$SO2x.gumbel[vec_indices[1+4]])^b[3]
  mu_vec_4 <- a[5] * dataset$CO.gumbel[vec_indices[1+4]] + mu[5] * (dataset$CO.gumbel[vec_indices[1+4]])^b[5]
  sigma_vec_1 <- sigma.o[1] * (dataset$O3.gumbel[vec_indices[1+4]])^b[1]
  sigma_vec_2 <- sigma.o[2] * (dataset$NO2x.gumbel[vec_indices[1+4]])^b[2]
  sigma_vec_3 <- sigma.o[4] * (dataset$SO2x.gumbel[vec_indices[1+4]])^b[3]
  sigma_vec_4 <- sigma.o[5] * (dataset$CO.gumbel[vec_indices[1+4]])^b[5]
  
  temp <- temp-sum(log(sigma_vec_1) + 0.5 * ((dataset$O3.gumbel[vec_indices[1+4]] - mu_vec_1)/sigma_vec_1)^2)
  temp <- temp-sum(log(sigma_vec_2) + 0.5 * ((dataset$NO2x.gumbel[vec_indices[1+4]] - mu_vec_2)/sigma_vec_2)^2)
  temp <- temp-sum(log(sigma_vec_3) + 0.5 * ((dataset$SO2x.gumbel[vec_indices[1+4]] - mu_vec_3)/sigma_vec_3)^2)
  temp <- temp-sum(log(sigma_vec_4) + 0.5 * ((dataset$CO.gumbel[vec_indices[1+4]] - mu_vec_4)/sigma_vec_4)^2)
  
  
  # CO
  mu_vec_1 <- a[1] * dataset$O3.gumbel[vec_indices[1+5]] + mu[1] * (dataset$O3.gumbel[vec_indices[1+5]])^b[1]
  mu_vec_2 <- a[2] * dataset$NO2x.gumbel[vec_indices[1+5]] + mu[2] * (dataset$NO2x.gumbel[vec_indices[1+5]])^b[2]
  mu_vec_3 <- a[3] * dataset$SO2x.gumbel[vec_indices[1+5]] + mu[3] * (dataset$SO2x.gumbel[vec_indices[1+5]])^b[3]
  mu_vec_4 <- a[4] * dataset$NO.gumbel[vec_indices[1+5]] + mu[4] * (dataset$NO.gumbel[vec_indices[1+5])^b[4]
  sigma_vec_1 <- sigma.o[1] * (dataset$O3.gumbel[vec_indices[1+5]])^b[1]
  sigma_vec_2 <- sigma.o[2] * (dataset$NO2x.gumbel[vec_indices[1+5]])^b[2]
  sigma_vec_3 <- sigma.o[4] * (dataset$SO2x.gumbel[vec_indices[1+5]])^b[3]
  sigma_vec_4 <- sigma.o[4] * (dataset$NO.gumbel[vec_indices[1+5]])^b[4]
  
  temp <- temp-sum(log(sigma_vec_1) + 0.5 * ((dataset$O3.gumbel[vec_indices[1+5]] - mu_vec_1)/sigma_vec_1)^2)
  temp <- temp-sum(log(sigma_vec_2) + 0.5 * ((dataset$NO2x.gumbel[vec_indices[1+5]] - mu_vec_2)/sigma_vec_2)^2)
  temp <- temp-sum(log(sigma_vec_3) + 0.5 * ((dataset$SO2x.gumbel[vec_indices[1+5]] - mu_vec_3)/sigma_vec_3)^2)
  temp <- temp-sum(log(sigma_vec_4) + 0.5 * ((dataset$NO.gumbel[vec_indices[1+5]] - mu_vec_4)/sigma_vec_4)^2)
  
  return(temp)
}
