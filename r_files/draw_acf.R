
acf_trawl <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50){
  seq_kappa <- seq(kappa, end_seq, by = delta)
  
  b_h_minus_0 <- - alpha  * (1-exp(-rho*h))
  b_0_minus_h <- - alpha  * (1-exp(-rho*h))
  b_0_h <- - alpha * exp(-rho*h)
  
  res <- 0
  first_mom <- 0
  res_0 <- 0
  beta <- beta
  
  # first_mom <- 0
 
  for(x in seq_kappa){
    x <- x + delta / 2
    # first_mom <- first_mom + (1+x/beta)^{-alpha}*(1+y/beta)^{-alpha}
    for(y in seq_kappa){
      y <- y + delta / 2
      res <- res + (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
      res_0 <- res_0 + (1+(x+y)/beta)^{-alpha}
    }
  }
  
  res <- res * delta^2
  res_0 <- res_0 * delta^2
  # first_mom <- first_mom * delta
  first_mom_sq <-  ((1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha - 1))^2
  # print(res/res_0)
  # print(first_mom^2)
  # res <- (res - first_mom^2)^1 # first_mom_sq
  # res_0 <- (res_0 - first_mom^2)^1 # first_mom_sq
  return((res-first_mom_sq)/(res_0-first_mom_sq))
}

d_plus <- function(alpha, beta, kappa){
  val <- beta^2 / ((alpha - 2) * (alpha - 1))
  val <- val*(log(1+2*kappa/beta) + (2*alpha-3)/((alpha-2)*(alpha-1)))
  val <- val*(1+2*kappa/beta)^{2-alpha}
  return(val)
}

d_times <- function(alpha, beta, kappa){
  val <- 2*beta / ((alpha - 2) * (alpha - 1))
  val <- val*(beta*(1+2*kappa/beta)^{2-alpha}*log(1+kappa/beta) + 
                ev.trawl::Zeta(kappa = kappa, beta = beta, alpha = alpha-1))
  return(val)
}

d_plus_times <- function(alpha, beta, kappa){
  val <- beta/((alpha-2)*(alpha-1))*log(1+2*kappa/beta)*(beta*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + ev.trawl::Zeta(alpha=alpha-1, beta=beta, kappa=kappa))
  val_2 <-beta/(alpha*(alpha-1))*log(1+kappa/beta)*ev.trawl::Zeta(alpha = alpha+1, beta=beta, kappa=kappa)
  val_3 <- beta^2/(alpha^3*(alpha-1))*(1+kappa/beta)^{-2*alpha}
  return(val + val_2 + val_3)
}

d_2_plus <-function(alpha, beta, kappa){
  val <- beta^2 / ((alpha - 2) * (alpha - 1))
  val <- val*(log(1+2*kappa/beta)^2 + 2/(alpha-2)*log(1+2*kappa/beta) + 2/(alpha-2)^2)
  val <- val*(1+2*kappa/beta)^{2-alpha}
  
  val_2 <- 2*beta^2/((alpha - 2) * (alpha - 1)^2)
  val_2 <- val_2*(log(1+2*kappa/beta) + 1/(alpha-2) + 1/(alpha-1))*(1+2*kappa/beta)^{2-alpha}
  return(val + val_2)
}

d_2_times <- function(alpha, beta, kappa){
  val <- 2*beta^2/((alpha-1)*(alpha-2))*((alpha-1)*log(+kappa/beta)^2*(1+2*kappa/beta)^{2-alpha} + log(1+kappa/beta)*(1+2*kappa/beta)^{1-alpha} + ev.trawl::Zeta(alpha=alpha,beta=beta,kappa=kappa))
  val <- val + 2*d_plus_times(alpha = alpha, beta = beta, kappa = kappa)
  return(val)
}

d_2_plus_times <- function(alpha, beta, kappa){
  val <- 4*beta/(alpha-1)^2 * ((log(1+2*kappa/beta) + log(1+kappa/beta) + 1/(alpha-1) ) * (1+2*kappa/beta)^{1-alpha} + ev.trawl::Zeta(alpha=alpha,beta=beta,kappa=kappa))
  val <- val + 2*beta^2/((alpha-1)*(alpha-2))*(log(1+2*kappa/beta)^2 + 2/(alpha-2)*log(1+2*kappa/beta)+2/(alpha-2)^2)*(1+2*kappa/beta)^{2-alpha}
}

var_latent_trawl <- function(alpha, beta, kappa){
  val <- 1/(alpha-2)*(1+2*kappa/beta)^{2-alpha} - 1/(alpha-1)*(1+kappa/beta)^{2-2*alpha}
  val <- val * beta^2/(alpha-1)
  return(val)
}
  
acf_trawl_approx <- function(h, alpha, beta, kappa, rho){
  acf_vals <- vapply(h, 
                     function(x){
                        var_value <- var_latent_trawl(alpha,beta,kappa)
                        first_order <- rho * alpha * (d_plus(alpha,beta,kappa)-d_times(alpha,beta,kappa))*x/sqrt(var_value)
                        second_order_1 <-  alpha*(d_2_plus(alpha,beta,kappa)-2*d_2_plus_times(alpha,beta,kappa)+d_2_times(alpha,beta,kappa))
                        second_order <- rho^2*alpha * (d_plus(alpha,beta,kappa)-d_times(alpha,beta,kappa) + second_order_1)*x^2
                        return(1+first_order)#+second_order)
  }, 1)
 return(acf_vals)
}

al <- 1/val_p3[1,1]
be <- val_p3[1,2]/val_p3[1,1] - val_p3[1,4]
ka <- val_p3[1,4]

plot(1:10/200, vapply(1:10/200, function(h){
  acf_trawl(h, alpha = al, beta = be, kappa = ka, 
            rho = 0.1, delta = 0.5)}, 1), lty=2, type='l')
vapply(1:10/200, function(h){
  acf_trawl(h, alpha = al, beta = be, kappa = ka, 
            rho = 0.1, delta = 0.2)}, 1)
lines(1:10/200, acf_trawl_approx(1:10/200, al,be,ka, 0.1))
acf_trawl_approx(1:10/200, al,be,ka, 0.1)

library('evir')
val_copy <- val_p3
par(mfrow=c(6,4), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:6){
  alpha_tmp <- 1/val_p3[i,1]
  beta_tmp <- val_p3[i,2]/val_p3[i,1] - val_p3[i,4]
  kappa_tmp <- val_p3[i,4]
  # rho_tmp <- val_p3[i,3]
  
  depth <- 20
  kk <- acf(epd[,i], lag.max = depth, plot=F)
  n_trials <- 50
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  rho_tab <- seq(log(1e-4), log(20), length.out = n_trials) %>% exp
  index <- 1
  acf(epd[,i], lag.max = depth, ylab=paste('ACF', colnames(epd)[i]))
  for(rho_iter in rho_tab){
    acf_vals <- vapply(c(0.05, 1:depth), function(h){
      acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                rho = rho_iter, delta = 1)}, 1)
    mae_tab[index] <- sum(abs(kk$acf - acf_vals)) 
    mse_tab[index] <- sum((kk$acf - acf_vals)^2) #+ 
      # sum(abs(kk$acf - vapply(c(0.05, 1:depth), function(h){
      #             acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
      #                       rho = rho_iter)}, 1)))
    if(index %% 2 == 0){
      acf_vals %>% (function(x){lines(c(0.05, 1:depth), x, col = 'red')})
    }
    index <- index + 1
  }
  
  
  #plot(rho_tab %>% log, mse_tab)
  plot(rho_tab %>% log, mse_tab, type='b', 
       xlab=paste('log-', expression(rho), 'green=MAE, orange=MSE'), 
       ylab='L1-Error with ACF', col = 'orange', lty = 2)
  lines(rho_tab %>% log, mae_tab, type='b', col = 'green', lty = 2)
  rho_tmp <- rho_tab[which.min(mse_tab)]
  print(rho_tmp)
  print(rho_tab[which.min(mae_tab)])
  val_copy[i,3] <- rho_tmp %>% log 
  # # ACF plots with red line given by model
  # acf(epd[,i], lag.max = depth)
  # vapply(c(0.05, 1:depth), function(h){
  #   acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
  #             rho = rho_tmp)}, 1) %>% (function(x){lines(c(0.05, 1:depth), x, col = 'red')})
  # 
  # Histogram
  hist(epd[,i][epd[,i] >0 & epd[,i] < 20], 
       breaks = 100, probability = T, main = '', 
       ylab=paste('Density', colnames(epd)[i]),
       xlab=paste(colnames(epd)[i], 'values'))
  lines(0:2000/20, dgpd(0:2000/20, shape = val_p3[i,1], scale=val_p3[i,2]), col = 'red')
  
  # QQ-plot
  data <- epd[,i][epd[,i] > 0]
  xi <- val_p3[i,1]
  if (xi == 0) {
    y <- qexp(ppoints(data))
  }
  if (xi != 0) {
    y <- evir::qgpd(ppoints(data), xi = xi, beta = val_p3[i,2])
  }
  plot(sort(data), y, xlab = "", ylab = "")
  #title(xlab = paste(colnames(epd)[i], "(ordered)"), ylab = paste('GPD quantile', expression(xi), '=',round(xi,3)))
  abline(lsfit(sort(data), y))
  title(xlab=paste(colnames(epd)[i], '(ordered)'), ylab=paste('GPD Quantile'))
  abline(v = quantile(epd[,i][epd[,i] >0], 0.95), col = 'blue', lwd = 1, lty = 2)
}


par(mfrow=c(6,4), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:6){
  alpha_tmp <- 1/val_p3[i,1]
  beta_tmp <- val_p3[i,2]/val_p3[i,1] - val_p3[i,4]
  kappa_tmp <- val_p3[i,4]
  rho_tmp <- val_p3[i,3]
  
  depth <- 20
  kk <- acf(epd[,i], lag.max = depth, plot=F)

  acf(epd[,i], lag.max = depth, ylab=paste('ACF', colnames(epd)[i]), cex.lab=1.5, cex.axis=1.5)
  lines(c(0.05, 1:depth), vapply(c(0.05, 1:depth), function(h){
    acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
              rho = rho_tmp, delta = 1)}, 1), col='red', lty=2, lwd=2)
  
  depth <- 20
  kk <- acf(epd[,i], lag.max = depth, plot=F)
  n_trials <- 50
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  rho_tab <- seq(log(1e-4), log(20), length.out = n_trials) %>% exp
  index <- 1
  # acf(epd[,i], lag.max = depth, ylab=paste('ACF', colnames(epd)[i]))
  for(rho_iter in rho_tab){
    acf_vals <- vapply(c(0.05, 1:depth), function(h){
      acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                rho = rho_iter, delta = 1)}, 1)
    mae_tab[index] <- sum(abs(kk$acf - acf_vals)) 
    mse_tab[index] <- sum((kk$acf - acf_vals)^2) #+ 
    # sum(abs(kk$acf - vapply(c(0.05, 1:depth), function(h){
    #             acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
    #                       rho = rho_iter)}, 1)))
    # if(index %% 2 == 0){
    #   acf_vals %>% (function(x){lines(c(0.05, 1:depth), x, col = 'red')})
    # }
    index <- index + 1
  }
  
  
  #plot(rho_tab %>% log, mse_tab)
  plot(rho_tab %>% log, mse_tab, type='o', 
       xlab= paste('log-rho'),#expression(rho), 
       ylab='Error with ACF', col = 'orange', lty = 2, cex=1.2, pch = 3, cex.lab=1.5, cex.axis=1.5)
  lines(rho_tab %>% log, mae_tab, type='o', col = 'green', lty = 2, cex=1.2, pch = 9)
  legend("topright", 
         legend = c("MAE", "MSE"), 
         col = c('green', 
                 'orange'), 
         pch = c(3,9), 
         bty = "n", 
         pt.cex = 2, 
         cex = 1.2, 
         text.col = "black", 
         horiz = F) 
         #inset = c(0.1, 0.1))
  abline(v=log(rho_tmp), col ='red', lty=2, lwd=2)
  
  #rho_tmp <- rho_tab[which.min(mse_tab)]
 
  # # ACF plots with red line given by model
  # acf(epd[,i], lag.max = depth)
  # vapply(c(0.05, 1:depth), function(h){
  #   acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
  #             rho = rho_tmp)}, 1) %>% (function(x){lines(c(0.05, 1:depth), x, col = 'red')})
  # 
  
  
  # Histogram
  hist(epd[,i][epd[,i] >0 & epd[,i] < 20], 
       breaks = 100, probability = T, main = '', 
       ylab=paste('Density', colnames(epd)[i]),
       xlab=paste(colnames(epd)[i], 'values'),cex.lab=1.4, cex.axis=1.4)
  lines(0:2000/20, eva::dgpd(0:2000/20, shape = val_p3[i,1], scale=val_p3[i,2]), col = 'red', lty=2, lwd=2)
  
  # QQ-plot
  data <- epd[,i][epd[,i] > 0]
  xi <- val_p3[i,1]
  if (xi == 0) {
    y <- qexp(ppoints(data))
  }
  if (xi != 0) {
    y <- evir::qgpd(ppoints(data), xi = xi, beta = val_p3[i,2])
  }
  plot(sort(data), y, xlab = "", ylab = "", cex.lab=1.5, cex.axis=1.5)
  #title(xlab = paste(colnames(epd)[i], "(ordered)"), ylab = paste('GPD quantile', expression(xi), '=',round(xi,3)))
  abline(lsfit(sort(data), y))
  title(xlab=paste(colnames(epd)[i], '(ordered)'), ylab=paste('GPD Quantile'), cex.lab=1.4)
  abline(v = quantile(epd[,i][epd[,i] >0], 0.95), col = 'blue', lwd = 1, lty = 4)
}

print(val_copy[,3] %>% exp)


kk <- acf(epd[,i], lag.max = 20, plot=F)
n_trials <- 50
mse_tab <- rep(0, n_trials)
rho_tab <- seq(log(1e-6), log(20), length.out = n_trials) %>% exp
index <- 1
for(rho_iter in rho_tab){
  mse_tab[index] <- sum((kk$acf - vapply(c(0.05, 1:20), function(h){
    acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
              rho = rho_iter)}, 1))^2)
  index <- index + 1
}
plot(rho_tab %>% log, mse_tab)
rho_tab[which.min(mse_tab)]


k <- 0

acf(epd[,1], lag.max = 40)
for(rho in seq(log(0.0067), log(20), length.out = 10) %>% exp){
  a0 <- acf_trawl(0.05, alpha = 1/0.2, beta = 5, kappa = 4, rho = rho)
  a_more <- vapply(1:40, function(h){
    acf_trawl(h, alpha = 1/0.15, beta = 5, kappa = 4, rho = rho)}
    , 1)
  a_more <- c(a0, a_more)
  a_more <-  a_more
  lines(0:40, a_more / a_more[1], col = rgb(green = 0,red = (k+1) / 11 , blue = 0, alpha=0.8))
  k <- k + 1
}

a_more
a0

start_params <- c(0.146896, 1.639081, -4, 3.824522) 

(1+3/(2.14/0.24 - 3))^{-1/0.23}

fn_to_optim <- function(x){
  # x_list <- list(alpha=x[1],
  #                beta=(x[2]),
  #                rho=x[3],
  #                kappa=x[4])
  x_list <- list(alpha=1/x[1],
                 beta=x[2]/abs(x[1]) - x[4],
                 rho=exp(x[3]),
                 kappa=x[4])
  
  if(x_list[["beta"]] <= 0){
    return(1e8)
  }
  
  
  return(-log(abs(x[1])^{-3}) - (x[3]) - 
           ParamsListFullPL(
             times = 1:5000/5000, 
             values = epd[1:5000,1], 
             delta = 5,
             params = x_list, logscale = T, 
             transformation = F))
}

par(mfrow=c(6,1), mar=c(4.1,4.1,0.5,0.5))
for(i in 1:6){
  start_params <- val_copy[i,]
  
  fn_to_optim <- function(x){
    # x_list <- list(alpha=x[1],
    #                beta=(x[2]),
    #                rho=x[3],
    #                kappa=x[4])
    x_list <- list(alpha=1/x[1],
                   beta=x[2]/abs(x[1]) - x[4],
                   rho=(x[3]),
                   kappa=x[4])
    
    if(x_list[["beta"]] <= 0){
      return(1e8)
    }
    
    
    return(-log(abs(x[1])^{-3}) - 
             ev.trawl::ParamsListFullPL(
               times = 1:length(epd[,i]), 
               values = epd[1:8000,i], 
               delta = s.clusters[i],
               params = x_list, 
               logscale = T, 
               transformation = T))
  }
  
  rho_vals <- seq(-2, 1, length.out = 20)
  pl_vals <- vapply(rho_vals, function(x){
    start_params[3] <- exp(x)
    return(fn_to_optim(start_params))
  }, 1)
  plot(rho_vals, pl_vals, xlab= expression(rho), ylab=paste('- log-PL', colnames(epd)[i]))
  print(rho_vals[which.min(pl_vals)])
  
  # plot(seq(-5, 3, length.out = 20) , vapply(seq(-10, 3, length.out = 20), function(x){
  #   start_params[3] <- (x)
  #   return(fn_to_optim(start_params))
  # }, 1), xlab= expression(rho), ylab=' - log-PL')
}



