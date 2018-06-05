setwd("~/GitHub/multi-trawl-extremes/r_files/")
source('latent_novel_simulation.R')
source('pairwise_latent_trawl.R')
source("imp_latent_noven_final.R")
source("latent_novel_simulation.R")
source("imp_warren.R")


# with taylor expansion
require(hypergeo)
alpha <- 4
beta <- 4
rho <- 0.2
kappa <- 3.11

zeta <- function(alpha, beta, kappa){
  res.zeta <- (alpha-1)*beta^alpha / ((beta+kappa)^{alpha-1} * alpha)
  res.zeta <- res.zeta * hypergeo(A = alpha-1, B = alpha-1, C = alpha, z = - beta/(beta+kappa))
  return(Re(res.zeta))
}

get_t1 <- function(alpha, beta, kappa, rho){
  # TODO adapt to other trawl functions
  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta)^{1-alpha} + (2*alpha-3)/((alpha-2)*(alpha-1)))
  d_times <- beta * log(1+kappa/beta)/(alpha*(alpha-1)-1)*(alpha*beta/(alpha-2)*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + alpha/(alpha-2)*zeta(alpha-1,beta,kappa) + zeta(alpha+1,beta,kappa))
  return(-rho * (d_plus - d_times))
}

get_estimate_rho <- function(alpha, beta, kappa, index, data){
  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta)^{1-alpha} + (2*alpha-3)/((alpha-2)*(alpha-1)))
  d_times <- beta * log(1+kappa/beta)/(alpha*(alpha-1)-1)*(alpha*beta/(alpha-2)*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + alpha/(alpha-2)*zeta(alpha-1,beta,kappa) + zeta(alpha+1,beta,kappa))

  d_times <- 2*beta/((alpha-2)*(alpha-1))*((1+2*kappa/beta)^{2-alpha}*log(1+kappa/beta)+zeta(alpha = alpha-1, beta = beta, kappa = kappa))
  return(var(data) / (index * (d_plus-d_times)))
}

alpha <- 4
beta <- 4
rho <- 0.2
kappa <- 3.11

# testing procedure for rho
latent_1 <- read.csv("latent_4_4_3.113118_0.2_1000_v3.csv", sep = " ")
test_rho <- rep(0,500)
for(i in 1:500){
  test_rho[i] <- get_estimate_rho(alpha = alpha, beta = beta, kappa = kappa, 5, latent_1[,i])
}

# barplot(vapply(1:40, function(i){return(mean(as.numeric(latent_1[1,]) * as.numeric(latent_1[i,])) - mean(as.numeric(latent_1[i,]))^2)}, 1.0))
# mean()

par(mfrow=c(2,2), mar=c(4.02,4,1,2.42))
hist(test_rho, breaks = 30, main = "", probability = T, xlab = expression(rho))
abline(v=(rho), col = "red", lwd = 2)
hist(log(test_rho), breaks = 30, main = "", probability = T, xlab = expression(log(rho)))
abline(v=log(rho), col = "red", lwd = 2)


# NO TRF
## CASE 1
# General parameters
alpha_true <- 4
beta_true <- 4
rho_true <- 0.2
kappa_true <- 3.113118

latent_1 <- read.csv("latent_4_4_3.113118_0.2_1000_v3.csv", sep = " ")
#latent_1 <- latent_1[,1:30]
mom_l1 <- initial_guess_trawl(latent_1)
index<-2
n_timestamps <- 4400 - 1 
n_sims <- 1000
mat_res <- matrix(0, nrow = n_sims, ncol = 4)
mat_res_o <- read.table("mat_res_4_4_id_1_372", header = F)
mat_res[1:372,] <- as.matrix(mat_res_o)
dim(mat_res)

for(index in 373:n_sims){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(latent_1[1:n_timestamps,index][latent_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(latent_1[1:n_timestamps,index] > 0))/length(latent_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*params_to_work_with[1]
  params_to_work_with[3] <- if(mom_l1$rho[index] != Inf) mom_l1$rho[index] else mom_l1$mean_rho
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = latent_1[1:n_timestamps,index],
                                          delta = 4,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = F)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(2.5,0.5,exp(-2.5),exp(0.2))
  upper_b <- c(10,10,exp(-0.01),exp(3))
  
  params_to_work_with
  lower_b
  upper_b
  
  system.time(res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, factr=1e12),
               method = "L-BFGS-B", lower = lower_b, upper=upper_b))
  print(index)
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}

#write.table(case_44, row.names = F, col.names = F, file = "mat_res_4_4_id_1_500.csv", sep = ",")

results <- mat_res[1:500,]
boxplot(results[,4])


## CASE 2
# General parameters
alpha_true <- 4
beta_true <- 1
rho_true <- 0.05
kappa_true <- 0.75

latent_1 <- read.csv("latent_9_1_0.75_0.05_500.csv", sep = " ")
# latent_1 <- latent_1[,1:2]
mom_l1 <- initial_guess_trawl(latent_1)
index<-1
n_timestamps <- 4400 - 1 
n_sims <- 500
mat_res <- matrix(0, nrow = n_sims, ncol = 4)

for(index in 1:n_sims){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(latent_1[1:n_timestamps,index][latent_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(latent_1[1:n_timestamps,index] > 0))/length(latent_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*params_to_work_with[1]
  params_to_work_with[3] <- if(mom_l1$rho[index] != Inf) mom_l1$rho[index] else mom_l1$mean_rho
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  #print(params_to_work_with)
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = latent_1[1:n_timestamps,index],
                                          delta = 6,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = F)
  # lower_b <- 0.8*params_to_work_with
  # lower_b[3] <- -2.3
  # upper_b <- 1.3*params_to_work_with
  # upper_b[3] <- -0.001
  # 
  lower_b <- c(1,
               0.1,
               0.01,
               0.2
               )
  upper_b <- c(17,
               15,
               0.8,
               1.0
               )
  
  params_to_work_with
  lower_b
  upper_b
  
  res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=3, factr=1e12),
               method = "L-BFGS-B", lower = lower_b, upper=upper_b)
  print(index)
  print(params_to_work_with)
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}

## CASE 3
# General parameters
alpha_true <- -4
beta_true <- 4
rho_true <- 0.2
kappa_true <- 9

latent_1 <- read.csv("latent_minus_4_4_9_0.2_500.csv", sep = ",")
# latent_1 <- latent_1[,1:2]
mom_l1 <- initial_guess_trawl(latent_1)
index<-3
n_timestamps <- 4400 - 1 
n_sims <- 500
mat_res <- matrix(0, nrow = n_sims, ncol = 4)

test_rho <- rep(0,500)
for(i in 1:500){
  test_rho[i] <- get_estimate_rho(alpha = 4, beta = 12, kappa = kappa_true, 7,
                                  trf_inv_g(z = latent_1[,i], alpha = alpha_true, beta = beta_true, 
                                            kappa = kappa_true, offset_scale = 12+kappa_true, offset_shape = 4))
}
boxplot(test_rho)
summary(test_rho)
evir::gpd(trf_inv_g(z = latent_1[,i], alpha = alpha_true, beta = beta_true, 
                    kappa = kappa_true, offset_scale = 3+kappa_true, offset_shape = 3), threshold = 0)$par.ests

for(index in 443:n_sims){ # TODO change 43
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(latent_1[1:n_timestamps,index][latent_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(latent_1[1:n_timestamps,index] > 0))/length(latent_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- fit_marginal[2]*abs(params_to_work_with[1])
  params_to_work_with[4] <- params_to_work_with[2]*(p_nz^{1/params_to_work_with[1]}-1)/ (p_nz^{1/params_to_work_with[1]})
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  
  
  beta_trf <- params_to_work_with[4]/(p_nz^{-1/4}-1)
  params_to_work_with[3] <- get_estimate_rho(alpha = 4, beta = beta_trf, kappa = params_to_work_with[4], 5, 
                                             trf_inv_g(z = latent_1[,i], alpha = params_to_work_with[1], beta = params_to_work_with[2],
                                                       kappa = params_to_work_with[4], offset_scale = beta_trf+params_to_work_with[4], offset_shape = 4))
  params_to_work_with[3] <- min(params_to_work_with[3], 1.0)
  if(is.na(params_to_work_with[3])){
    params_to_work_with[3] <- runif(min = 0.1,0.5, n = 1)
  }
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = latent_1[1:n_timestamps,index],
                                          delta = 4,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = T)
  # lower_b <- 0.8*params_to_work_with
  lower_b[3] <- -2.3
  upper_b <- 1.3*params_to_work_with
  upper_b[3] <- -0.001

  lower_b <- c(-8,
               1,
               0.01,
               1
  )
  upper_b <- c(-1,
               15,
               0.6,
               20
  )
  
  params_to_work_with
  lower_b
  upper_b
  
  #print(params_to_work_with)
  res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=1, factr=1.5e11),
               method = "L-BFGS-B", lower = lower_b, upper=upper_b)
  print(index)
  print(params_to_work_with)
  print(res$par)
  mat_res[index,] <- res$par
  #print(mat_res[index,])
}

#write.table(mat_res, row.names = F, col.names = F, file = "mat_res_minus_4_4_0.2_9_500.csv", sep = ",")

## CASE 4
# General parameters
alpha_true <- -4
beta_true <- 1
rho_true <- 0.1
kappa_true <- 15

latent_1 <- read.csv("latent_minus_4_1_15_0.1_500.csv", sep = ",")
# latent_1 <- latent_1[,1:2]
mom_l1 <- initial_guess_trawl(latent_1)
index <- 3
n_timestamps <- 4400 - 1 
n_sims <- 500
mat_res <- matrix(0, nrow = n_sims, ncol = 4)

for(index in 1:n_sims){
  params_to_work_with <- rep(0, 4)
  
  # Find Kappa
  fit_marginal <-  fExtremes::gpdFit(latent_1[1:n_timestamps,index][latent_1[1:n_timestamps,index] > 0], u= 0)@fit$fit$par
  p_nz <- length(which(latent_1[1:n_timestamps,index] > 0))/length(latent_1[1:n_timestamps,index])
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[4] <- fit_marginal[2]*abs(params_to_work_with[1])
  params_to_work_with[2] <- (1-p_nz^(-1/params_to_work_with[1]))/(p_nz^(-1/params_to_work_with[1]))
  params_to_work_with[4] <- params_to_work_with[4] - params_to_work_with[2]
  
  beta_trf <- params_to_work_with[4]/(p_nz^{-1/4}-1)
  params_to_work_with[3] <- get_estimate_rho(alpha = 4, beta = beta_trf, kappa = params_to_work_with[4], 13, 
                                             trf_inv_g(z = latent_1[,i], alpha = params_to_work_with[1], beta = params_to_work_with[2],
                                                       kappa = params_to_work_with[4], offset_scale = beta_trf+params_to_work_with[4], offset_shape = 4))
  params_to_work_with[3] <- min(params_to_work_with[3], 1.0)
  if(is.na(params_to_work_with[3])){
    params_to_work_with[3] <- runif(min = 0.1,0.5, n = 1)
  }
  mat_res[index,] <- params_to_work_with
  params_to_work_with
  
  ##TODO CHECK NAMES OF VARIABLES !!!!
  fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:n_timestamps,
                                          values = latent_1[1:n_timestamps,index],
                                          delta = 6,
                                          lambda = 0.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = T)
  # lower_b <- 0.8*params_to_work_with
  lower_b[3] <- -2.3
  upper_b <- 1.3*params_to_work_with
  upper_b[3] <- -0.001

  lower_b <- c(-10,
               0.05,
               0.01,
               5
  )
  upper_b <- c(-1,
               5,
               0.6,
               25
  )
  
  # params_to_work_with
  # lower_b
  # upper_b
  # 
  res <- optim(fn = fn_to_optim, par = params_to_work_with[c(1,2,3,4)], control = list(trace=1, factr=1.5e11),
               method = "L-BFGS-B", lower = lower_b, upper=upper_b)
  print(index)
  print(params_to_work_with)
  print(res$par)
  mat_res[index,] <- res$par
}

#write.table(mat_res, row.names = F, col.names = F, file = "mat_res_minus_4_1_0.1_15_500.csv", sep = ",")


#library(DEoptim)
#DEoptim::DEoptim(fn_to_optim, lower = lower_b, upper_b, control = DEoptim.control(trace=T, itermax=200))

fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:4399,
                                        values = latent_1[1:4399,index],
                                        delta = 4,
                                        lambda = 0.0,
                                        model_vars_names = univ_model_vars_names,
                                        fixed_names = c("kappa", "rho"),
                                        fixed_params = c(3.11, 0.2),#c(params_to_work_with[c(2,3,4)]),
                                        logscale = T,
                                        transformation = F)
#plot(seq(0.1,2,length.out = 20), vapply(seq(0.1,2,length.out = 20), fn_to_optim, 1))
n_points_grid <- 30
alpha_s <- seq(0.5,6,length.out = n_points_grid)
beta_s <- seq(0.5,6,length.out = n_points_grid)
vis_mat <- matrix(0, ncol=n_points_grid, nrow=n_points_grid)

for(index_a in 1:n_points_grid){
  for(index_b in 1:n_points_grid){
    vis_mat[index_a,index_b] <- fn_to_optim(c(alpha_s[index_a], beta_s[index_b]))
  }
  print(index_a)
}
par(mfrow=c(1,2))
library(plot3D)
persp3D(alpha_s, beta_s, vis_mat, theta = 30, phi = 30, xlab= expression(alpha), ylab="beta", zlab="(minus) log-PL")
library(ggplot2)
contour(x = alpha_s, y = beta_s, z = vis_mat, xlab="Alpha", ylab="Beta", nlevels = 40, zlim=c(0,110000), col = colorRamp(c("red", "green"))(90))
lines(x=alpha_s/6*6, y =alpha_s/8*exp(alpha_s/8)*6-0.7, col = "red", lwd = 2)
par(mfrow=c(1,1))



### TRF Plot

fn_to_optim <- loglikelihood_pl_univ_ic(times = 1:4399,
                                        values = latent_1[1:4399,index],
                                        delta = 4,
                                        lambda = 0.0,
                                        model_vars_names = univ_model_vars_names,
                                        fixed_names = c("kappa", "rho"),
                                        fixed_params = c(9, 0.2),#c(params_to_work_with[c(2,3,4)]),
                                        logscale = T,
                                        transformation = T)
#plot(seq(0.1,2,length.out = 20), vapply(seq(0.1,2,length.out = 20), fn_to_optim, 1))
n_points_grid <- 30
alpha_s <- seq(-10,-1,length.out = n_points_grid)
beta_s <- seq(1,10,length.out = n_points_grid)
vis_mat <- matrix(0, ncol=n_points_grid, nrow=n_points_grid)

for(index_a in 1:n_points_grid){
  for(index_b in 1:n_points_grid){
    vis_mat[index_a,index_b] <- fn_to_optim(c(alpha_s[index_a], beta_s[index_b]))
  }
  print(index_a)
}

par(mfrow=c(1,2))
library(plot3D)
persp3D(alpha_s, beta_s, vis_mat, theta = -10, phi = 10, xlab= expression(alpha), ylab="beta", zlab="(minus) log-PL")
library(ggplot2)
contour(x = alpha_s, y = beta_s, z = vis_mat, xlab="Alpha", ylab="Beta", nlevels = 50, zlim=c(12000,15000), col = colorRamp(c("red", "green"))(90))
lines(x=alpha_s/6*6, y = abs(alpha_s*2)*exp(alpha_s/16)*2-10, col = "red", lwd = 2)
par(mfrow=c(1,1))












quire(grDevices) # for colours
x <- -6:16
op <- par(mfrow = c(2, 2))
contour(outer(x, x), method = "edge", vfont = c("sans serif", "plain"))
z <- outer(x, sqrt(abs(x)), FUN = "/")
image(x, x, z)
contour(x, x, z, col = "pink", add = TRUE, method = "edge",
        vfont = c("sans serif", "plain"))
contour(x, x, z, ylim = c(1, 6), method = "simple", labcex = 1,
        xlab = quote(x[1]), ylab = quote(x[2]))
contour(x, x, z, ylim = c(-6, 6), nlev = 20, lty = 2, method = "simple",
        main = "20 levels; \"simple\" labelling method")
par(op)

## Persian Rug Art:
x <- y <- seq(-4*pi, 4*pi, len = 27)
r <- sqrt(outer(x^2, y^2, "+"))
opar <- par(mfrow = c(2, 2), mar = rep(0, 4))
for(f in pi^(0:3))
  contour(cos(r^2)*exp(-r/f),
          drawlabels = FALSE, axes = FALSE, frame = TRUE)

rx <- range(x <- 10*1:nrow(volcano))
ry <- range(y <- 10*1:ncol(volcano))
ry <- ry + c(-1, 1) * (diff(rx) - diff(ry))/2
tcol <- terrain.colors(12)
par(opar); opar <- par(pty = "s", bg = "lightcyan")
plot(x = 0, y = 0, type = "n", xlim = rx, ylim = ry, xlab = "", ylab = "")
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = tcol[8], border = "red")
contour(x, y, volcano, col = tcol[2], lty = "solid", add = TRUE,
        vfont = c("sans serif", "plain"))
title("A Topographic Map of Maunga Whau", font = 4)
abline(h = 200*0:4, v = 200*0:4, col = "lightgray", lty = 2, lwd = 0.1)

## contourLines produces the same contour lines as contour
plot(x = 0, y = 0, type = "n", xlim = rx, ylim = ry, xlab = "", ylab = "")
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = tcol[8], border = "red")
contour(x, y, volcano, col = tcol[1], lty = "solid", add = TRUE,
        vfont = c("sans serif", "plain"))
line.list <- contourLines(x, y, volcano)
invisible(lapply(line.list, lines, lwd=3, col=adjustcolor(2, .3)))
par(opar)


## NO TRF boxplots
case_44 <- read.table("mat_res_4_4_id_1_500.csv", sep = ",")
case_91 <- read.table("mat_res_9_1_0.05_0.75_500.csv", sep = ",")

par(mfrow=c(2,4), mar=c(1.02,2.82,2.82,3.42))
name_d <- c(expression(alpha), expression(beta), expression(rho), expression(kappa))
true_44 <- c(4,4,0.2,3.11)
true_91 <- c(4,1,0.05,0.75)
for(idex in 1:4){
  boxplot(case_44[,idex], main = (name_d[idex]))
  abline(h = true_44[idex], col = "red")
}
for(idex in 1:4){
  if(idex == 3){
    boxplot(case_91[,idex]/4, main = (name_d[idex]))
  }else{
    boxplot(case_91[,idex], main = (name_d[idex]))
  }
  
  abline(h = true_91[idex], col = "red")
}
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))

par(mfrow=c(1,1), mar=c(4.1,4.5,1.1,2.1))
summary(lm(case_44[,2] + case_44[,4] ~ case_44[,1]))
summary(lm(case_91[,2] + case_91[,4]~ case_91[,1]))
plot((case_44[,2]+case_44[,2])/case_44[,1], xlab="Simulation", ylab=expression((beta+kappa)/alpha), ylim=c(-0.5,4), pch = 6)
points((case_91[,2]+case_91[,4])/case_91[,1], xlab=expression(alpha), ylab=expression(beta), pch = 0)
abline(h=0.48 - 0.22/2, col = "blue", lwd = 3, lty = 2)
abline(h=2.12-1.26/8, col = "darkgreen", lwd = 3, lty = 2)
legend(420, 4.08, c("Case 1", "Case 2"), col = c("black", "black", 6), lty = c(NA, NA, 1), lwd = c(1,1), pch = c(6, 0, 4),
       merge = TRUE, bg = "gray95")
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))



## TRF boxplots
case_44 <- read.table("mat_res_minus_4_4_0.2_9_500.csv", sep = ",")
case_41 <- read.table("mat_res_minus_4_1_0.1_15_500.csv", sep = ",")

par(mfrow=c(2,4), mar=c(1.02,2.82,2.82,3.42))
name_d <- c(expression(alpha), expression(beta), expression(rho), expression(kappa))
true_44 <- c(-4,4,0.2,9)
true_41 <- c(-4,1,0.1,15)
for(idex in 1:4){
  if(idex == 3){
    boxplot(case_44[,idex], main = (name_d[idex]), ylim= c(0.1,0.4))
  }else{
    boxplot(case_44[,idex], main = (name_d[idex]))
  }
  abline(h = true_44[idex], col = "red")
}
for(idex in 1:4){
  if(idex == 3){
    boxplot(case_41[,idex], main = (name_d[idex]), ylim= c(0.01,0.3))
  }else{
    boxplot(case_41[,idex], main = (name_d[idex]))
  }
  
  abline(h = true_41[idex], col = "red")
}
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))

par(mfrow=c(1,1), mar=c(4.1,4.5,1.1,2.1))
summary(lm(case_44[,2] + case_44[,4] ~ case_44[,1]))
summary(lm(case_41[,2] + case_41[,4]~ case_41[,1]))
plot((case_44[,2]+case_44[,4])/case_44[,1], xlab="Simulation", ylab=expression((beta+kappa)/alpha), ylim=c(-5.5,-1.5), pch = 6)
points((case_41[,2]+case_41[,4])/case_41[,1], xlab=expression(alpha), ylab=expression(beta), pch = 0, lwd=2)
abline(h=-1.63-4.49/5.8, col = "blue", lwd = 3, lty = 2)
abline(h=- 2.56 - 3.35/3, col = "darkgreen", lwd = 3, lty = 2)
legend(440, -1.5, c("Case 3", "Case 4"), col = c("black", "black", 6), lty = c(NA, NA, 1), lwd = c(1,2), pch = c(6, 0, 4),
       merge = TRUE, bg = "gray95")
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))


# estimating the ACF

alpha <- 4
beta <- 1
rho <- 0.05
kappa <- 0.75

latent_1 <- read.csv("latent_9_1_0.75_0.05_500.csv", sep = " ")
test_rho <- rep(0,500)
for(i in 1:500){
  test_rho[i] <- get_estimate_rho(alpha = alpha, beta = beta, kappa = kappa, 15, latent_1[,i])
}

hist(test_rho, breaks = 30, main = "", probability = T, xlab = expression(rho))
abline(v=(rho), col = "red", lwd = 2)
hist(log(test_rho), breaks = 30, main = "", probability = T, xlab = expression(log(rho)))
abline(v=log(rho), col = "red", lwd = 2)

par(mfrow=c(1,1))

trawl_list_grp <- collection_trawl(times = 1:50, params = list(rho=rho), type = "exp", prim = F)
trawl_list_prim_grp <- collection_trawl(times = 1:50, params = list(rho=rho), type = "exp", prim = T)
tp1 <- 0.0
tp2 <- 0.0
tp3 <- 0.0
for(i in 1:200){
  gen_exc_grp <- rltrawl(alpha = alpha,
                         beta = beta,
                         kappa = kappa,
                         times = 1:50,
                         trawl_fs = trawl_list_grp,
                         trawl_fs_prim = trawl_list_prim_grp,
                         n = 1,
                         transformation = F)
  tp1 <- tp1 + exp(-kappa*(gen_exc_grp[10]+gen_exc_grp[11]))/(gen_exc_grp[10]*gen_exc_grp[11])
  tp2 <- tp2 + exp(-kappa*(gen_exc_grp[10]+gen_exc_grp[12]))/(gen_exc_grp[10]*gen_exc_grp[12])
  tp3 <- tp3 + exp(-kappa*(gen_exc_grp[10]+gen_exc_grp[13]))/(gen_exc_grp[10]*gen_exc_grp[13])
}
tp1/200
tp2/200
tp3/200

rho_univ <- optim(par=params_to_work_with[4],
                 fn = fn_to_optim,
                 control = list(trace=4, maxit=10, pgtol=1e-3, parscale=rep(0.5, 1)),
                 method = "L-BFGS-B",
                 lower=-3,
                 upper=1)
rho_univ$par


gen_exc <- rlexceed(alpha = alpha,
                    beta = beta,
                    kappa = kappa,
                    times = times,
                    trawl_fs = trawl_list,
                    trawl_fs_prim = trawl_list_prim,
                    n = 1,
                    transformation = F)

exceed_p
compute_mean_std <- function(values_array){
  # values_array contains the time series with first axis as time and second as # of time series
  n_dim <- length(values_array[1,])
  n_values <- length(values_array[,1])
  results <- rep(0, n_dim)
  for(index in 1:n_dim){
    results[index] <- length(which(values_array[,index] > 0))/ n_values
  }
  
  return(list(mean=mean(results), sd=sd(results), prob=results))
}

compute_mean_std(values_array = exceed_p)
alpha_values <- seq(alpha/5, alpha*3, length.out = 100)
beta_values <- seq(beta/5, beta*3, length.out = 100)
kappa_values <- seq(kappa/5, kappa*3, length.out = 100)
paras_mom <- mom_gpd(exceed_p)

fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c(),
                                                           fixed_params = c(),
                                                           params = c(exp(x), beta+kappa),
                                                           model_vars_names = c("alpha", "beta"),
                                                           logscale = T,
                                                           transformation = F,
                                                           n_moments = 4))}
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("alpha", "beta"),
                                                           fixed_params = c(alpha, beta),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha","beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res <- optim(par = c(0.3*0.8), fn = fn_to_optim, method = "BFGS")
exp(res$par)
alpha
beta

# TRF

## alpha
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("beta", "kappa"),
                                                           fixed_params = c(beta, kappa),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha", "beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res_alpha <- optim(par = log(paras_mom$alpha[10]), fn = fn_to_optim, method = "BFGS")
exp(res_alpha$par)
alpha

plot(alpha_values, vapply(alpha_values,
                          function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                                      fixed_names = c("beta", "kappa"),
                                                                      fixed_params = c(beta, kappa),
                                                                      params = c(x),
                                                                      model_vars_names = c("alpha", "beta", "kappa"),
                                                                      logscale = T,
                                                                      transformation = T,
                                                                      n_moments = 4))},
                          1.0), type = "l", ylab="log-likelihood")
abline(v=alpha, col="red")
abline(v=exp(res_alpha$par[1]), col = "blue")

## beta
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("alpha", "kappa"),
                                                           fixed_params = c(alpha, kappa),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha", "beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res_beta <- optim(par = log(paras_mom$beta[10]), fn = fn_to_optim, method = "BFGS")
exp(res_beta$par)
beta
plot(beta_values, vapply(beta_values,
                          function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                                      fixed_names = c("alpha", "kappa"),
                                                                      fixed_params = c(alpha, kappa),
                                                                      params = c(x),
                                                                      model_vars_names = c("alpha", "beta", "kappa"),
                                                                      logscale = T,
                                                                      transformation = T,
                                                                      n_moments = 4))},
                          1.0), type = "l", ylab="log-likelihood")
abline(v=beta, col="red")
abline(v=exp(res_beta$par[1]), col = "blue")

## kappa
fn_to_optim <- function(x){return(-marginal_gpd_likelihood(values = trawl_p[,10][trawl_p[,10] > 0.0],
                                                           fixed_names = c("alpha", "beta"),
                                                           fixed_params = c(alpha, beta),
                                                           params = c(exp(x[1])),
                                                           model_vars_names = c("alpha", "beta", "kappa"),
                                                           logscale = T,
                                                           transformation = T,
                                                           n_moments = 4))}
res_kappa <- optim(par = log(0.8*0.3), fn = fn_to_optim, method = "BFGS")
exp(res_kappa$par)
kappa
plot(kappa_values, vapply(kappa_values,
                         function(x){return(-marginal_gpd_likelihood(values = trf_inv_g(z=trawl_p[,10][trawl_p[,10] > 0.0], alpha = alpha,
                                                                                        beta = beta, kappa = 0.25/0.75, 
                                                                                        offset_scale = trf_find_offset_scale(alpha = alpha, beta = beta, 
                                                                                              kappa = 0.25/0.75, offset_shape = 4)+x,
                                                                                        offset_shape = 4),
                                                                     fixed_names = c("alpha"),
                                                                     fixed_params = c(4),
                                                                     params = c(trf_find_offset_scale(alpha = alpha, beta = beta, 
                                                                                                      kappa = 0.25/0.75, offset_shape = 4)+x),
                                                                     model_vars_names = c("alpha", "beta"),
                                                                     logscale = T,
                                                                     transformation = F,
                                                                     n_moments = 4))},
                         1.0), type = "l", ylab="log-likelihood")
abline(v=kappa, col="red")
abline(v=exp(res_kappa$par[1]), col = "blue")

## ALL
ididi <- 100
fn_to_optim <- function(x){return(-marginal_simple_lik(values = latent_1[,ididi][latent_1[,ididi] > 0.0], params = x))}
optim(fn_to_optim, par = c(5,3), method = "L-BFGS-B", lower = c(0.1,0.1), upper = c(10,10))
gpd_fit(latent_1[,ididi][latent_1[,ididi] > 0.0], initial_guess =c(8,2))


# Initial guess
## MoM on GPD(alpha, beta+kappa)
## Use this and proba > 0 to get kappa: kappa ~ (beta+kappa)*(1-p^{1/alpha})
## Now we have alpha, beta and kappa individually


setwd("~/GitHub/multi-trawl-extremes/data/simulations/")
e_notrf_4000 <- read.csv(file = "exceed_4000a4b2k0.32485.csv")
dim(e_notrf_4000)
delta <- 4

ig_t_notrf_4000 <- initial_guess_trawl(e_notrf_4000)
params_to_work_with <- c((ig_t_notrf_4000$alpha[1]),
                         (ig_t_notrf_4000$beta[1]),
                         log(ig_t_notrf_4000$rho[1]),
                         log(ig_t_notrf_4000$kappa[1]))

fn_to_optim <- loglikelihood_pl_univ_ic(times = (1:4000),
                                        values = e_notrf_4000[1:4000,1],
                                        delta = delta,
                                        lambda = 0.0,
                                        model_vars_names = univ_model_vars_names,
                                        fixed_names = c("alpha"),
                                        fixed_params = c(4.09),
                                        logscale = T,
                                        transformation = F)

system.time(fn_to_optim(params_to_work_with))

optim(params_to_work_with[c(1)], fn_to_optim, method = "BFGS", 
      hessian = F, control = list(trace=5))$par

optim(params_to_work_with[c(2,3,4)], fn_to_optim, method = "L-BFGS-B", 
      hessian = F, control = list(trace=5), lower = c(0.5*2.5,-3,-3), upper = c(1.5*2.5,2,2))

# diagnostic plots
trial_alpha <- seq(0.01, 1, length.out = 15)
system.time(plot(trial_alpha, vapply(trial_alpha, fn_to_optim, 1)))

xs <- seq(0.1, 2, length.out = 50)
plot(xs, evir::dgpd(xs, xi = 1/5, beta = (0.1)/5))
lines(xs, evir::dgpd(xs, xi = 1/alpha, beta = (beta+kappa)/alpha))


# Warren process
lw_sim <- matrix(0, nrow = length(times), ncol = 10)
for(index in 1:10){
 lw_sim[,index] <- rwprocess(alpha = alpha,
                             beta = beta,
                             kappa = kappa,
                             rho = rho,
                             timesteps = length(times),
                             n=1)
}

lw_mom <- mom_gpd(lw_sim)
acf(lw_sim[,5])
