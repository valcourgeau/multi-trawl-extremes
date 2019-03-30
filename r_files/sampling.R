# Implementation of the edging up algorithm
library("VineCopula")
library("rlist")

get_name <-  function(j, tree, RVM) {
  M <- RVM$Matrix
  d <- nrow(M)
  # variable names
  print(cat("Copula:", M[c(j, (d - tree + 1):d), j], "\n"))
  nams <- RVM$names[M[c(j, (d - tree + 1):d), j]]
  # conditioned set
  bef <- paste(nams[2],
               nams[1],
               sep = ",",
               collapse = "")
  # conditioning set
  aft <- if (length(nams) > 2) {
    gsub(" ",  ",", do.call(paste, as.list(nams[3:length(nams)])))
  }  else ""
  # paste together
  sep <- if (length(nams) > 2) " ; " else ""
  paste(bef, aft, sep = sep, collapse = "")
}

get_copula_vars <-  function(j, tree, RVM) {
  M <- RVM$Matrix
  d <- nrow(M)
  # variable names
  vars <- M[c(j, (d - tree + 1):d), j]
  if(tree ==1){
    return(list(LHS=vars[1:2], RHS=NA))
  }else{
    return(list(LHS=vars[1:2], RHS=vars[3:length(vars)]))
  }
}

get_copula <-  function(j, tree, RVM) {
  M <- RVM$Matrix
  d <- nrow(M)
  # variable names
  cop_fam <- RVM$family[d-tree+1, j]
  cop_par1 <- RVM$par[d-tree+1, j]
  cop_par2 <- RVM$par2[d-tree+1, j]
  return(VineCopula::BiCop(family=cop_fam,
                    par = cop_par1,
                    par2 = cop_par2))
}

RVineCondSimSingle <- function(target, sample_target, RVM){
  d <- nrow(RVM$Matrix)
  sample_vine <- rep(NA, d)
  sample_vine[target] <- sample_target
  for(tree in 1:d){
    for(cop_number in 1:(d-tree)){
      already_sampled <- which(!is.na(sample_vine)) # updating already_sampled
      vars_lhs_rhs <- get_copula_vars(cop_number, tree, RVM)
      if(tree > 1){
        if(all(vars_lhs_rhs$RHS %in% already_sampled)){
          already_left <- which(vars_lhs_rhs$LHS %in% already_sampled)
          if(length(already_left) == 1){
            # Here we know if have sampled the conditional vars + ONE and only ONE of LHS so sample the remaining one.
            
            to_sample <- vars_lhs_rhs$LHS[-already_left]
            assertthat::are_equal(length(to_sample), 1)
            assertthat::are_equal(length(already_left), 1)
            sample_vine[to_sample] <- VineCopula::BiCopCondSim(N = 1,
                                                               cond.val = sample_vine[already_left],
                                                               cond.var = 1,
                                                               obj = get_copula(cop_number, tree, RVM))
          }
        }
      }else{
        # first tree level here
        already_left <- which(vars_lhs_rhs$LHS %in% already_sampled)
        if(length(already_left) == 1){
          # Here we know if have sampled the conditional vars + ONE and only ONE of LHS so sample the remaining one.
          to_sample <- vars_lhs_rhs$LHS[-already_left]
          assertthat::are_equal(length(to_sample), 1)
          assertthat::are_equal(length(already_left), 1)
          sample_vine[to_sample] <- VineCopula::BiCopCondSim(N = 1,
                                                             cond.val = sample_vine[vars_lhs_rhs$LHS[already_left]],
                                                             cond.var = already_left, # WARNING might need to be replaced by 1 instead
                                                             obj = get_copula(cop_number, tree, RVM))
        }
      }
    }
  }
  return(sample_vine)
}

RVineCondSim <- function(target, samples, RVM, nsub=1){
  N <- length(samples)
  d <- nrow(RVM$Matrix)
  # if(nsub == 1){
  #   sample_mat <- matrix(NA, nrow=N, ncol=)
  #   sample_mat[,target] <- samples
  # }else{
    sample_mat <- array(NA, dim = c(nsub, d, N))
    sample_mat[,target,] <- samples
  #}
  
  pb <- txtProgressBar(min = 0, max = nsub, initial = 0, char = "=",
                 width = 20, title, label, style = 3, file = "")
  for(sub in 1:nsub){
      sub_ind <- 0
      sample_mat[sub,,] <- vapply(samples, 
                                  function(x){
                                    sol <- RVineCondSimSingle(target = target, sample_target = x, RVM = RVM)
                                    sub_ind <- sub_ind + 1
                                    setTxtProgressBar(pb, sub + sub_ind/nsub)
                                    return(sol)
                                  }, FUN.VALUE = rep(0,d))
      
  }
  
  close(pb)
  if(nsub == 1){
    return(sample_mat[1,,])
  }else{
    return(sample_mat)
  }
}

GetConditionalMatrix <- function(data, h, conditional_on){
  n <- nrow(data)
  col_names <- colnames(data)
  data <- cbind(data[h:n,], data[1:(n-h+1),conditional_on])
  colnames(data) <- c(col_names, paste(col_names[conditional_on], "_ex", sep = ""))
  return(data)
}



GetPredictionResults <- function(simulations, conditional_data, conditional_data_pred, extreme_index_pred, extreme_quantiles, conditional_extreme_quantiles){
  prediction_results <- list()
  nvars <- ncol(conditional_data)-1
  var_names <- colnames(conditional_data)[1:nvars]
  
  for(var in 1:nvars){
    counter <- 0
    counter_f <- 0
    counter_t <- 0
    n_examples <- dim(simulations)[3]
    
    n_extreme <- 0 
    n_no_extreme <- 0
    
    for(k in 1:n_examples){
      binary_response <- simulations[,var,k] > conditional_extreme_quantiles[var]
      nsub <- length(simulations[,var,k])
      proportion_extreme <- length(which(binary_response)) / nsub
      
      mean_extreme_response <- mean(binary_response) # same as proportion_extreme
      sd_extreme_response <- sd(binary_response)
      
      lower_bound <- mean_extreme_response - 1.96*sd_extreme_response/sqrt(nsub)
      upper_bound <- mean_extreme_response + 1.96*sd_extreme_response/sqrt(nsub)
      
      prediction <- lower_bound > 1-conditional_extreme_quantiles[var]
      pred_upper <- upper_bound < 1-conditional_extreme_quantiles[var]
      
      actual <- conditional_data_pred[extreme_index[k], var] > extreme_quantiles[var]
      
      if(actual){
        n_extreme <- n_extreme + 1
      }else{
        n_no_extreme <- n_no_extreme + 1
      }
      
      # if-else structure for prediction, sensitivity and specificity.
      if(prediction == actual){
        if(actual == F){
          counter_f <- counter_f + 1
        }
        else{
          if(actual == T){
            counter_t <- counter_t + 1
          }
        }
        counter <- counter + 1
      }else{
        if(pred_upper == actual){
          if(actual == F){
            counter_f <- counter_f + 1
          }
          else{
            if(actual == T){
              counter_t <- counter_t + 1
            }
          }
          counter <- counter + 1
        }
      }
    }
    
    speci <- round(counter_f / n_no_extreme, 2)
    if(n_no_extreme == 0){
      speci <- 1
    }
    sensi <- round(counter_t / n_extreme, 2)
    if(n_extreme == 0){
      sensi <- 1
    }
    
    cat('Extreme in', colnames(data_for_pred)[var], '\t h =', h_selected, ' \n')
    cat('Accuracy \t \t \t', counter / n_examples, ' \n')
    cat('Specificity \t\t\t', speci, '\t found ', counter_f, '\t total', n_no_extreme,  '\n')
    cat('Sensitivity\t\t\t', sensi, '\t found ', counter_t, '\t total', n_extreme,  '\n')
    
    prediction_results[[var_names[var]]] <- list(accuracy=counter / n_examples,
                                                 specificity=speci,
                                                 sensitivity=sensi,
                                                 n_extreme_found=counter_t,
                                                 n_not_extreme_found=counter_f,
                                                 n_extreme=n_extreme,
                                                 n_not_extreme=n_no_extreme)
  }
  
  return(prediction_results)
}



setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
daily_bloom <- read.csv("/../Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/hourly_bloomsbury_2000_2017.csv")

library(evir)
library(forecast)
library(lubridate)
stlpd <- daily_bloom
n_vars <- length(stlpd[1,]) - 3
for(i_agent in 4:(n_vars+3)){
  fitting_matrix <- cbind(cos(2*pi*1:length(stlpd$index)/200),
                          sin(2*pi*1:length(stlpd$index)/200),
                          cos(2*pi*1:length(stlpd$index)/14),
                          sin(2*pi*1:length(stlpd$index)/14),
                          #as.numeric(isWeekend(stlpd$date)==T))
                          vapply(1:3, FUN = function(i){quarter(stlpd$date) == i}, FUN.VALUE = quarter(as.Date(stlpd$date))),
                          #vapply(1:12, FUN = function(i){month(as.Date(stlpd$date)) == i}, FUN.VALUE = month(as.Date(stlpd$date))),
                          vapply(1:6, FUN = function(i){wday(as.Date(stlpd$date)) == i}, FUN.VALUE = wday(as.Date(stlpd$date))),
                          vapply(1:23, FUN = function(i){hour(as.Date(lubridate::hm(stlpd$time))) == i}, FUN.VALUE = wday(as.Date(stlpd$date))))
  
  fit <- lm(stlpd[,i_agent] ~ fitting_matrix)
  summary(fit)
  fitting_indices <- which(summary(fit)$coefficients[,4] < 0.05)
  if(1 %in% fitting_indices){
    fitting_indices <- fitting_indices[-1]
  }
  fitting_matrix <- fitting_matrix[,fitting_indices-1]
  print(fitting_indices)
  stlpd[,i_agent] <- lm(stlpd[,i_agent] ~ fitting_matrix)$residuals
}
#daily_bloom <- daily_bloom[,-c(1:3)]
daily_bloom <- stlpd[,-c(1:3)]

train_index <- 1:as.integer(0.8*nrow(daily_bloom))
daily_train <- daily_bloom[train_index,]
daily_test <- daily_bloom[-train_index,]


tron_air_pollution <- computeTRON(data = daily_train,
                          p.zeroes = rep(0.95, 6),
                          horizons = c(1,2,3,4,5,6,12,24),
                          clusters = c(10,15,11,11,13,13),
                          n_samples = 100000,
                          name_matrices_file = "matrix_air_pollution_rerun",
                          name_vine_file = "vines_air_pollution_rerun",
                          name_tron_file = "tron_air_pollution_rerun",
                          save = T,
                          sparse = F)



vines_air_pollution <- rlist::list.load("vines_air_pollution_rerun_vines.RData")

horizons_vines <- c(1,2,3,4,5,6,12,24)

mat_res <- matrix(0, ncol=ncol(daily_bloom), nrow=length(horizons_vines))
k <- 1
for(h in horizons_vines){
  mat_res[k,] <- tron_air_pollution[[h]]$mean[2,]
  k <- k + 1
}
mat_res

library(VineCopula)
par(mfrow=c(2,3), mar=c(4.1,4.1,0.5,0.5))
vines_of_interest <-vines_air_pollution[[12]][1][[1]]
vines_of_interest$names <- c(colnames(daily_train), 'O3 ex')
plot(vines_of_interest,  type=1, 
     interactive=F, label.bg="white", label.col="black", label.cex = 1.2, edge.lwd=1.3,
     edge.labels=c("family-par"), edge.label.cex=1.5, edge.label.col="blue")

contour(vines_of_interest)


vines_1_72 <- rlist::list.load("hourly-bloomsbury-vines-12361224-v2.RData")

h_selected <- 12
#v_temp <- vines_1_72[[h_selected]][[1]]

#v_temp <- RVineMatrixNormalize(v_temp) # vars in order from d to 1
v_temp # we condition on the first variable 

alpha <- 4.1
beta <- 4.12
kappa <- 4.5
rho <- 0.04


conditional_on <- 2
v_temp <- vines_air_pollution[[h_selected]][[conditional_on]]
data_for_fitting <- daily_train
data_for_pred <- daily_test

q_extremes <- apply(data_for_fitting, FUN=function(x){quantile(x, 0.95)}, MARGIN = 2)
q_extremes_plus_one <- c(q_extremes, q_extremes[conditional_on])

cond_data <- GetConditionalMatrix(data_for_fitting, h=h_selected, conditional_on = conditional_on)
cond_data_test <- GetConditionalMatrix(data_for_pred, h=h_selected, conditional_on = conditional_on)

extreme_index <- which(cond_data[,ncol(cond_data)] > q_extremes[conditional_on]) # finds rows on which last col is extreme
extreme_index_test <- which(cond_data_test[,ncol(cond_data_test)] > q_extremes[conditional_on]) # finds rows on which last col is extreme
ecdf_without <- apply(data_for_fitting, FUN = ecdf, MARGIN = 2)
ecdf_extreme <- apply(cond_data[extreme_index,], FUN = ecdf, MARGIN = 2)

epd <- apply(data_for_fitting, 
             FUN = function(x){
               (x - q_extremes) * (x > q_extremes)
             }, MARGIN = 1)
epd <- t(epd)
dim(epd)
sd_epd <- apply(epd, sd, MARGIN = 2)
epd <- apply(X = epd, MARGIN = 2, FUN = function(x){return(x/sd(x))})

db_unif <- vapply(data_for_pred[,conditional_on], 
                  function(x){
                    if(x <= q_extremes[conditional_on]){ecdf_without[[conditional_on]](x)}
                    else{
                  min(0.9999,ecdf_without[[conditional_on]](q_extremes[conditional_on]) + (1-q_extremes[conditional_on])*(1-(1+sign(alpha)*(x-q_extremes[conditional_on])/sd_epd[conditional_on]/(beta+kappa))^{-alpha}))
    }}, FUN.VALUE = 1)
hist(db_unif, breaks = 50)

set.seed(42)
q_extremes_cond <- vapply(1:ncol(data_train), function(x){mean(cond_data[,x] <= q_extremes[x])}, FUN.VALUE = 1)

ssmple <- RVineCondSim(target = conditional_on, samples = db_unif[extreme_index_test[1:200]], RVM = v_temp, nsub=100)

GetPredictionResults(simulations = ssmple,
                     conditional_data = cond_data,
                     conditional_data_pred = cond_data_test,
                     extreme_index_pred = extreme_index_test,
                     extreme_quantiles = q_extremes,
                     conditional_extreme_quantiles = q_extremes_cond)

GetMCPredictions <- function(data_train, data_test, marginal, p.zero, h, conditional_on, bootstrap, n_extreme=NA){
  alpha <- marginal$alpha
  beta <- marginal$beta
  kappa <- marginal$alpha
  
  q_extremes <- apply(data_for_fitting, FUN=function(x){quantile(x, p.zero)}, MARGIN = 2)
  q_extremes_plus_one <- c(q_extremes, q_extremes[conditional_on])
  
  cond_data <- GetConditionalMatrix(data_train, h=h_selected, conditional_on = conditional_on)
  cond_data_test <- GetConditionalMatrix(data_test, h=h_selected, conditional_on = conditional_on)
  
  
  extreme_index <- which(cond_data[,ncol(cond_data)] > q_extremes[conditional_on]) # finds rows on which last col is extreme
  extreme_index_test <- which(cond_data_test[,ncol(cond_data_test)] > q_extremes[conditional_on]) # finds rows on which last col is extreme
  ecdf_without <- apply(data_train, FUN = ecdf, MARGIN = 2)
  ecdf_extreme <- apply(cond_data[extreme_index,], FUN = ecdf, MARGIN = 2)
  
  epd <- apply(data_train, 
               FUN = function(x){
                 (x - q_extremes) * (x > q_extremes)
               }, MARGIN = 1)
  epd <- t(epd)
  dim(epd)
  sd_epd <- apply(epd, sd, MARGIN = 2)
  epd <- apply(X = epd, MARGIN = 2, FUN = function(x){return(x/sd(x))})
  
  db_unif <- vapply(data_for_pred[,conditional_on], 
                    function(x){
                      if(x <= q_extremes[conditional_on]){ecdf_without[[conditional_on]](x)}
                      else{
                        min(0.9999,ecdf_without[[conditional_on]](q_extremes[conditional_on]) + (1-q_extremes[conditional_on])*(1-(1+sign(alpha)*(x-q_extremes[conditional_on])/sd_epd[conditional_on]/(beta+kappa))^{-alpha}))
                      }}, FUN.VALUE = 1)
  q_extremes_cond <- vapply(1:ncol(data_train), function(x){mean(cond_data[,x] <= q_extremes[x])}, FUN.VALUE = 1)
  
  if(is.na(n_extreme)){
    ssmple <- RVineCondSim(target = conditional_on, samples = db_unif[extreme_index_test], RVM = v_temp, nsub=boostrap)
  }else{
    ssmple <- RVineCondSim(target = conditional_on, samples = db_unif[extreme_index_test[1:n_extreme]], RVM = v_temp, nsub=100)
  }
  
  res <- GetPredictionResults(simulations = ssmple,
                         conditional_data = cond_data,
                         conditional_data_pred = cond_data_test,
                         extreme_index_pred = extreme_index_test,
                         extreme_quantiles = q_extremes,
                         conditional_extreme_quantiles = q_extremes_cond)
  return(list(results=res, bootstrap=ssmaple)) 
}

set.seed(42)
GetMCPredictions(data_train = data_for_fitting, data_test = data_for_pred,
                 marginal = list(alpha=alpha, beta=beta, kappa=kappa),
                 p.zero = 0.95,
                 h = 12,
                 conditional_on = 2,
                 bootstrap = 50,
                 n_extreme = 50)



for(var in 1:ncol(data_for_pred)){
  counter <- 0
  counter_f <- 0
  counter_t <- 0
  n_examples <- dim(ssmple)[3]
  
  for(k in 1:n_examples){
    binary_response <- ssmple[,var,k] > q_extremes_cond[var]
    nsub <- length(ssmple[,var,k])
    proportion_extreme <- length(which(binary_response)) / nsub
    
    mean_extreme_response <- mean(binary_response) # same as proportion_extreme
    sd_extreme_response <- sd(binary_response)
    
    #cat(mean_extreme_response, ' - ', cond_data[k,3] > q_extremes_cond[3], '\n')
    lower_bound <- mean_extreme_response - 1.96*sd_extreme_response/sqrt(nsub)
    upper_bound <- mean_extreme_response + 1.96*sd_extreme_response/sqrt(nsub)
    
    prediction <- lower_bound > 1-q_extremes_cond[var]
    pred_upper <- upper_bound < 1-q_extremes_cond[var]
    
    actual <- cond_data[extreme_index[k], var] > q_extremes[var]
    #cat(prediction, ' - ', actual, '\n')
    #cat('Pred:', prediction, ' - actual:', actual, '\n')
    if(prediction == actual){
      if(actual == F){
        counter_f <- counter_f + 1
      }
      else{
        if(actual == T){
          counter_t <- counter_t + 1
        }
      }
      counter <- counter + 1
    }else{
      if(pred_upper == actual){
        if(actual == F){
          counter_f <- counter_f + 1
        }
        else{
          if(actual == T){
            counter_t <- counter_t + 1
          }
        }
        counter <- counter + 1
      }
    }
  }
  cat('Extreme in', colnames(data_for_pred)[var], '\t h =', h_selected, ' \n')
  cat('Accuracy \t \t \t', counter / n_examples, ' \n')
  n_no_extreme <- sum(cond_data[extreme_index[1:n_examples],var] <= q_extremes[var])
  n_extreme <- sum(cond_data[extreme_index[1:n_examples],var] > q_extremes[var])
  cat('Specificity \t\t\t', round(counter_f / n_no_extreme, 2), '\t found ', counter_f, '\t total', n_no_extreme,  '\n')
  cat('Sensitivity\t\t\t', round(counter_t / n_extreme, 2), '\t found ', counter_t, '\t total', n_extreme,  '\n')
}

# find all params
# find new quantile under extreme in other var
# use the vine condition on CO?
       