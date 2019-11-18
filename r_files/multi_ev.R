# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)
require(stats)
require(evir)
source("prep_univariate_latent_trawl_fit.R")
source("utils.R")
source("infer_latent_value.R")
source("auto_threshold_selection.R")



#' Wrapper from \code{ev.trawl} GenerateParameters to fit univariate latent
#' trawl extreme values model.
#'
#' @param data cleaned dataset.
#' @param clusters_size Extreme clusters size (integer) vector .
#' @param thresholds Thresholds values for extremes.
#' @param optim Logical (default=TRUE). Whether to execute BFGS optimisation
#'   with fixed alpha.
#' @param name a string for the file name once saved. Default is NA, if NA, it
#'   creates a timestamp format.
#' @param save Logical flag (default is TRUE). Whether to save the result as
#'   .RData.
#' @return Set of thresholds values as large as \code{data} after which values
#'   in data are considered extremes.
findUnivariateParams <- function(data, clusters_size, thresholds, optim=T, name=NA, save=T){
  n_vars <- length(data[1,])
  n_rows <- length(data[,1])
  
  # temp_params <-vapply(rlist::list.load("2019-3-1-15-43-27_params.RData"), as.vector, c(1,1,1,1)) %>% t
  # print(dim(temp_params))
  # return(vapply(rlist::list.load("2019-3-1-15-43-27_params.RData"), as.vector, c(1,1,1,1)) %>% t)
  # 
  if(any(clusters_size <= 0)){
    stop('clusters_size should have positive entries.')
  }
  
  if(length(clusters_size) == 1){
    clusters_size <- rep(clusters_size, n_vars)
  }else{
    if(n_vars != length(clusters_size)){
      stop('clusters_size should be either a scalar to apply to every column or as large as data.')
    }
  }
  
  val_params <- matrix(0, nrow = n_vars, ncol = 4)
  exc <- makeExceedances(data = data, thresholds = thresholds, normalize = TRUE)
  for(i_agent in 1:n_vars){
      cat("Generating univ MLEs for", colnames(data)[i_agent], "...")

    print(eva::gpdFit(exc[,i_agent][exc[,i_agent] > 0], 0.0))
    val_params[i_agent,1:2] <- eva::gpdFit(exc[,i_agent][exc[,i_agent] > 0], 0.0)$par.sum$Estimate[2:1]
    p <- length(which(exc[,i_agent] > 0))/length(exc[,i_agent])
    val_params[i_agent,4] <- (1 - p^{val_params[i_agent,1]}) * val_params[i_agent,2]/abs(val_params[i_agent,1])
    
    val_params[i_agent, 3] <- ev.trawl::GetEstimateRho(alpha = 1/val_params[i_agent, 1],
                                                       beta = val_params[i_agent, 2]/abs(val_params[i_agent, 1]) - val_params[i_agent, 4],
                                                       kappa = val_params[i_agent, 4],
                                                       cluster.size = 20,
                                                       data = exc[,i_agent]) 
    
    val_params[i_agent, 3] <- (val_params[i_agent, 3])
    
    #val_params[i_agent, 3] <- val_copy[i,3]
    
    #val_params[i_agent, 3] <- log(val_params[i_agent, 3])
      cat(" done.\n")
      if(optim){
        cat("|---> Preparing optimisation...")
        marginal_values <- exc[1:15000,i_agent] # TODO 1000!
        marginal_times <- (1:length(marginal_values))
        
        exp_params_names <- c("alpha", "beta", "rho", "kappa")
        fixed_params_names <- c('alpha', 'beta', 'kappa')
        fixed_params_index <- which(exp_params_names %in% fixed_params_names)
        initial_guess <- list()
        for(j in 1:length(val_params[1,])){
          if(!(j %in% fixed_params_index)){
            initial_guess[[exp_params_names[j]]] <- val_params[i_agent,j] 
          }
        }
        
        fn_to_optim <- function(x){
          x_list <- list(alpha=1/val_params[i_agent,1],
                         beta=val_params[i_agent,2]/abs(val_params[i_agent,1]) - val_params[i_agent,4],
                         rho=x,
                         kappa=val_params[i_agent,4])

          if(x_list[["beta"]] <= 0){
            return(1e8)
          }
          
          return(-log(abs(x[1])^{-3})  - 
                   ParamsListFullPL(
                    times = marginal_times, 
                    values = marginal_values, 
                    delta = clusters_size[i_agent],
                    params = x_list, logscale = T, 
                    transformation = T))
        }
        
        
        multiplier_bottom <- 0.8
        multiplier_top <- 1.2
        
        # if(val_params[i_agent,1] > 0){
        #   lower_b <-c(
        #     0.9 * abs(val_params[i_agent, 1]),
        #     multiplier_bottom * abs(val_params[i_agent, 2]),
        #             0.001,
        #             multiplier_bottom * abs(val_params[i_agent, 4]))
        #   upper_b <-c(
        #     1.1 * abs(val_params[i_agent, 1]),
        #     multiplier_top * abs(val_params[i_agent, 2]) ,
        #             5*val_params[i_agent, 3],
        #     multiplier_top * abs(val_params[i_agent, 4]))
        # }else{
        #   lower_b <-c(
        #     - multiplier_top * abs(val_params[i_agent, 1]),
        #     -multiplier_top*abs(val_params[i_agent, 1]),
        #     -4,
        #     0.1)
        #   upper_b <-c(
        #     -multiplier_bottom * abs(val_params[i_agent, 1]),
        #     multiplier_top * abs(val_params[i_agent, 2]),
        #     -0.5,
        #     multiplier_top * abs(val_params[i_agent, 2]))
        # }
        
        lower_b <- exp(-4)
        upper_b <- exp(1)
        
        cat(" initialise...")
        time_to_cv <- proc.time()[3]
        cat('\n')
        print(val_params[i_agent,])
        print(upper_b)
        print(lower_b)
        # res <- DEoptim::DEoptim(fn = fn_to_optim,
        #                         lower = lower_b,
        #                         upper = upper_b,
        #                         control = DEoptim::DEoptim.control(
        #                           NP=100,
        #                           itermax=3,
        #                           strategy = 1)
        #                         )$optim$bestmem
        # res_old <- res
        
        #upper_b[3] <- min(4 * abs(res[3]), 0.3)
        
        res <- val_params[i_agent,3]
        res <- stats::optim(fn_to_optim,
                            par = res,
                            method = "L-BFGS-B",
                            lower = lower_b, upper = upper_b,
                            control = list(maxit=100,
                                           factr=1e7,
                                           trace=3)
                      )$par
          
        cat("done in ", proc.time()[3]-time_to_cv, "s.\n", sep = "")
        res_concat <- rep(0, 4)
        res_concat[fixed_params_index] <- val_params[i_agent, fixed_params_index]
        if(length(fixed_params_index) > 0){
        res_concat[-fixed_params_index] <- res
        }else{
          res_concat <- res
        }
        cat("old", val_params[i_agent,] %>% as.vector,"\n")
        val_params[i_agent, ] <- res_concat
        cat("new",val_params[i_agent, ] %>% as.vector,"\n")
      }
      
  }
  
  cols_names <- colnames(data)
  if(save){
    val_params_list <- list()
    for(i in 1:nrow(val_params)){
      val_params_list[[cols_names[i]]] <- val_params[i,]
    }
    rlist::list.save(val_params_list, 
                     makeFileName(file_name = name, tag = "_params", extension = ".RData") )
  }
  
  return(val_params)
}

CustomMarginalMLE <- function(data){
  init_guess <- eva::gpdFit(data, 0.0)$par.sum$Estimate

  fn_mle <- function(par){
    data_for_mle <- data
    p <- length(which(data_for_mle>0)) / length(data_for_mle)
    kap <- (1-p^{par[1]})*par[2]/abs(par[1])
    p_non_zero <- 1-(1+kap/(par[2]/abs(par[1])-kap))^{-par[1]}
    like <- eva::dgpd(data_for_mle[data_for_mle>0.0], scale = par[2], shape=par[1],
                      loc=0, log.d = F) * p_non_zero

    log.like <- sum(log(like)) + length(data_for_mle == 0.0) * log(1-p_non_zero)
    return(-log.like)
  }
  
  return(stats::optim(par = init_guess[2:1], fn_mle, method='L-BFGS-B', lower=c(1e-03, 0.5), upper=c(2,20))$par)
}

#' Wrapper from \code{ev.trawl} GenerateParameters to fit univariate latent
#' trawl extreme values model.
#'
#' @param data cleaned dataset.
#' @param clusters_size Extreme clusters size (integer) vector .
#' @param thresholds Thresholds values for extremes.
#' @param optim Logical (default=TRUE). Whether to execute BFGS optimisation
#'   with fixed alpha.
#' @param name a string for the file name once saved. Default is NA, if NA, it
#'   creates a timestamp format.
#' @param save Logical flag (default is TRUE). Whether to save the result as
#'   .RData.
#' @return Set of thresholds values as large as \code{data} after which values
#'   in data are considered extremes.
findUnivariateParamsv2 <- function(data, clusters_size, thresholds, optim=T, name=NA, save=T){
  n_vars <- length(data[1,])
  n_rows <- length(data[,1])
  
  # temp_params <-vapply(rlist::list.load("2019-3-1-15-43-27_params.RData"), as.vector, c(1,1,1,1)) %>% t
  # print(dim(temp_params))
  # return(vapply(rlist::list.load("2019-3-1-15-43-27_params.RData"), as.vector, c(1,1,1,1)) %>% t)
  # 
  if(any(clusters_size <= 0)){
    stop('clusters_size should have positive entries.')
  }
  
  if(length(clusters_size) == 1){
    clusters_size <- rep(clusters_size, n_vars)
  }else{
    if(n_vars != length(clusters_size)){
      stop('clusters_size should be either a scalar to apply to every column or as large as data.')
    }
  }
  
  val_params <- matrix(0, nrow = n_vars, ncol = 4)
  exc <- makeExceedances(data = data, thresholds = thresholds, normalize = TRUE)
  for(i_agent in 1:n_vars){
    cat("Generating univ MLEs for", colnames(data)[i_agent], "...")
    
    #print(eva::gpdFit(exc[,i_agent][exc[,i_agent] > 0], 0.0))
    #val_params[i_agent,1:2] <- eva::gpdFit(exc[,i_agent][exc[,i_agent] > 0], 0.0)$par.sum$Estimate[2:1]
    
    val_params[i_agent,1:2] <- CustomMarginalMLE(exc[,i_agent])

    p <- length(which(exc[,i_agent] > 0))/length(exc[,i_agent])
    val_params[i_agent,4] <- (1 - p^{val_params[i_agent,1]}) * val_params[i_agent,2]/abs(val_params[i_agent,1])
    
    val_params[i_agent, 3] <- ev.trawl::GetEstimateRho(alpha = 1/val_params[i_agent, 1],
                                                       beta = val_params[i_agent, 2]/abs(val_params[i_agent, 1]) - val_params[i_agent, 4],
                                                       kappa = val_params[i_agent, 4],
                                                       cluster.size = 20,
                                                       data = exc[,i_agent]) 
    
    val_params[i_agent, 3] <- val_params[i_agent, 3]
  
    cat(" done.\n")
    if(optim){
      cat("|---> Preparing optimisation...")
      marginal_values <- exc[1:min(15000, n_rows), i_agent] # TODO 1000!
      marginal_times <- 1:length(marginal_values)
      
      exp_params_names <- c("alpha", "beta", "rho", "kappa")
      fixed_params_names <- c('alpha', 'beta', 'kappa')
      fixed_params_index <- which(exp_params_names %in% fixed_params_names)
      initial_guess <- list()
      for(j in 1:length(val_params[1,])){
        if(!(j %in% fixed_params_index)){
          initial_guess[[exp_params_names[j]]] <- val_params[i_agent,j] 
        }
      }
      
      fn_to_optim <- function(x){
        x_list <- list(alpha=1/val_params[i_agent,1],
                       beta=val_params[i_agent,2]/abs(val_params[i_agent,1]) - val_params[i_agent,4],
                       rho=x,
                       kappa=val_params[i_agent,4])
        
        if(x_list[["beta"]] <= 0){
          return(1e8)
        }
        
        return(-log(abs(x_list[['alpha']])^{-3})  - 
                 ev.trawl::ParamsListFullPL(
                   times = marginal_times, 
                   values = marginal_values, 
                   delta = clusters_size[i_agent],
                   params = x_list, logscale = T, 
                   transformation = T))
      }
      
      # plot(seq(0.001, 1, length.out = 10),
      #      vapply(seq(0.001, 1, length.out = 10), fn_to_optim, 1))
      
      
      multiplier_bottom <- 0.8
      multiplier_top <- 1.2
      
      # if(val_params[i_agent,1] > 0){
      #   lower_b <-c(
      #     0.9 * abs(val_params[i_agent, 1]),
      #     multiplier_bottom * abs(val_params[i_agent, 2]),
      #             0.001,
      #             multiplier_bottom * abs(val_params[i_agent, 4]))
      #   upper_b <-c(
      #     1.1 * abs(val_params[i_agent, 1]),
      #     multiplier_top * abs(val_params[i_agent, 2]) ,
      #             5*val_params[i_agent, 3],
      #     multiplier_top * abs(val_params[i_agent, 4]))
      # }else{
      #   lower_b <-c(
      #     - multiplier_top * abs(val_params[i_agent, 1]),
      #     -multiplier_top*abs(val_params[i_agent, 1]),
      #     -4,
      #     0.1)
      #   upper_b <-c(
      #     -multiplier_bottom * abs(val_params[i_agent, 1]),
      #     multiplier_top * abs(val_params[i_agent, 2]),
      #     -0.5,
      #     multiplier_top * abs(val_params[i_agent, 2]))
      # }
      
      lower_b <- exp(-4)
      upper_b <- exp(1)
      
      cat(" initialise...")
      time_to_cv <- proc.time()[3]
      cat('\n')
      print(val_params[i_agent,])
      print(upper_b)
      print(lower_b)
      # res <- DEoptim::DEoptim(fn = fn_to_optim,
      #                         lower = lower_b,
      #                         upper = upper_b,
      #                         control = DEoptim::DEoptim.control(
      #                           NP=100,
      #                           itermax=3,
      #                           strategy = 1)
      #                         )$optim$bestmem
      # res_old <- res
      
      #upper_b[3] <- min(4 * abs(res[3]), 0.3)
      
      res <- val_params[i_agent,3]
      res <- stats::optim(fn = fn_to_optim,
                          par = res,
                          method = "L-BFGS-B",
                          lower = lower_b, upper = upper_b,
                          control = list(maxit=100,
                                         factr=1e7,
                                         trace=0)
      )$par
      
      cat("done in ", proc.time()[3]-time_to_cv, "s.\n", sep = "")
      res_concat <- rep(0, 4)
      res_concat[fixed_params_index] <- val_params[i_agent, fixed_params_index]
      if(length(fixed_params_index) > 0){
        res_concat[-fixed_params_index] <- res
      }else{
        res_concat <- res
      }
      cat("old", val_params[i_agent,] %>% as.vector,"\n")
      val_params[i_agent, 3] <- res
      cat("new",val_params[i_agent, ] %>% as.vector,"\n")
    }
    
  }
  
  cols_names <- colnames(data)
  if(save){
    val_params_list <- list()
    for(i in 1:nrow(val_params)){
      val_params_list[[cols_names[i]]] <- val_params[i,]
    }
    rlist::list.save(val_params_list, 
                     makeFileName(file_name = name, tag = "_params", extension = ".RData") )
  }
  
  return(val_params)
}




#' This function returns the threshold(s) corresponding the quantile of probability p.exceed
#' @param data cleaned dataset
#' @param p.exceed threshold in probability to be considered an extreme (ex: 0.95)
#' @examples
#' d <- 3 
#' n <- 100
#' data <- matrix(runif(n*d), ncol=d)
#' getThresholds(data, rep(0.9, d))
#' getThresholds(data, 0.9)
#' getThresholds(data[,1], 0.9)
getThresholds <- function(data, p.exceed){
  if(any(p.exceed < 0) | any(p.exceed > 1)){
    stop('p.exceed should be between 0 and 1.')
  }
  
  # deals with data as a vector
  if(is.vector(data)){
    return(as.numeric(quantile(data, p.exceed[1])))
  }
  
  if(length(p.exceed) == 1){
    return(apply(data, MARGIN = 2, 
                 FUN = function(x){
                   as.numeric(quantile(x, p.exceed))
                   }))
  }else{
    if(length(data[1,]) == length(p.exceed)){
      return(sapply(1:length(data[1,]), 
                    function(i){as.numeric(quantile(data[,i], p.exceed[i]))
                      }))
    }else{
      stop('p.exceed should either be a scalar or as wide as data.')
    }
  }
}

#getThresholds(core_energy_data[,1:20], 0.95)

#' Wrapper to makeMatrix in the special case of extreme value thresholds.
#' We consider p.exceed as the threshold cut-off probability.
#' @param data cleaned dataset
#' @param p.exceed Threshold probability for extreme values (if higher, 
#'                 it is an extreme) 
#' @examples 
#' d <- 3
#' n <- 100
#' data <- matrix(runif(n*d), ncol=d)
#' p.exc <- 0.9
#' makeThresholdsMatrix(data = data, p.exceed = p.exc)
makeThresholdsMatrix <- function(data, p.exceed){
  return(makeMatrix(data = data, 
                    vector_to_rep = getThresholds(data, p.exceed = p.exceed)))
}



#' @param data clean dataset;
#' @param thresholds d dimensional vector of threshold values;
#' @param normalize Logical. Normalise extremes to have sd = 1;
#' @examples
#' threshold_data <- makeExceedances(test, thresholds = getThresholds(test, 0.8))
makeExceedances <- function(data, thresholds, normalize=TRUE){
  # TODO change thresholds to p.exceed?
  if(!is.vector(thresholds)){
    stop('Thresholds should be a vector of extreme threshold value.')
  }
  if(!is.vector(data)){
    if(length(data[1,]) != length(thresholds)){
      stop('thresholds and data have non-comforting width. 
         Tip: Use getThresholds to get the right size.')
    }
    rep_thres <- makeMatrix(data = data, vector_to_rep = thresholds) # TODO change to make makeThresholdsMatrix
    epd <- (data - rep_thres) * (data > rep_thres)
  }else{
    epd <- (data - thresholds) * (data > thresholds)
    
  }
  
  if(normalize){
    epd <- apply(X = as.matrix(epd), MARGIN = 2, 
                 FUN = function(x){return(x/sd(x))})
    
  }  
  epd <- data.frame(epd)
  colnames(epd) <- colnames(data)
  return(epd)
}


plgpd_unif_at_zero <- function(x, p.zero, alpha, beta, kappa){
  if(p.zero < 0 | p.zero > 1) stop("p.zero should be between 0 and 1.")
  if(x == 0)
    return(runif(n = 1, min = 0, max = p.zero))
  else{
    return(p.zero + (1-p.zero)*(1-max(0, (1+sign(alpha)*x/(beta+kappa))^{-(alpha)})))
  }
}

plgpd_unif_at_zero.row <- function(xs, p.zeroes, params.mat){
  # params.mat contains alpha beta rho kappa
  res <- rep(0, length(xs))
  for(i in 1:length(xs)){
    res[i] <- plgpd_unif_at_zero(x = xs[i],
                                 p.zero = p.zeroes[i],
                                 alpha = params.mat[i,1],
                                 beta = params.mat[i,2],
                                 kappa = params.mat[i,4])
  }

  return(res)
}

computePZero <- function(params){
  # Computes probability of having zero given univariate extreme value model
  return(1-(1+params[,4]/params[,2])^{-abs(params[,1])})
}

#'Constructs a collecttion of matrices with the primary set of keys being the
#'horizons and the second one being the column indices. Those matrices are the
#'data conditional on extreme value in the respectivel column.
#'
#'@param data clean dataset.
#'@param p.zeroes threshold in probability to be considered an extreme (ex:
#'  0.95).
#'@param conditional_on name or index of variable to condition the extremes on.
#'  Default if \code{NA} and in this case, it creates the full collection of
#'  matrices.
#'@param horizons (integer) sequence of selected extreme horizons.
#'@param clusters_size Extremes cluster size for each marginal.
#'@param save Logical (flag, default to FALSE). Whether we save the conditonal
#'  matrices.
#'@param n_samples Number of rows selected to create the matrices (from 1 to
#'  \code{n_samples}).
#'@param name Filename if the collection of matrices is saved. Default is NA
#'  (creates a timestamp on filename).
#'@param optim Logical (flag, default to TRUE). Whether to perform univ model
#'  optimisation fit on top of MLEs.
#'@return Double-leveled list with horizons as first set of keys and columns
#'  indices as second set of keys. E.g. \code{result[[2]][[4]]} would pick the
#'  second horizon on the 4th variable.
makeConditionalMatrices <- function(data, p.zeroes, conditional_on=NA,
                                    horizons, clusters_size, 
                                    n_samples = length(data[,1]), 
                                    save=F, name=NA, optim=T){
  n_vars <- length(data[1,])
  # adapt the size of p.zeroes to the number of cols
  if(length(p.zeroes) == 1){
    p.zeroes <- rep(p.zeroes, n_vars)
  }
  
  if(length(conditional_on) == 1){
    if(is.na(conditional_on)){
      conditional_on <- 1:n_vars
    }
  }else{
    if(any(!vapply(conditional_on, FUN = is.numeric, FUN.VALUE = T))){
      stop('conditional_on should be a list of column indices to select from.')
    }
    if(any(conditional_on < 1) | any(conditional_on > n_vars)){
      stop('conditional_on should be a list of index between 1 and the number of variables in data.')
    }
  }

  thres <- getThresholds(data = data, p.exceed = p.zeroes)
  print("thres")
  print(thres)
  
  # params <- findUnivariateParamsv2(data = data, clusters_size = clusters_size,
  #                                 thresholds = thres, name = name, save = T, optim = optim) # TODO WARNING save?!
  # write.table(params, "energy_weather_params")
  params <- read.table('energy_weather_params')
  params_copy <- params
  params_copy[,1] <- 1 / params[,1]
  params_copy[,2] <- params[,2] / params[,1] - params[,4]
  params_copy[,4] <- params[,4]
  
  params <- params_copy
  
  #params <- val_p3
  # params[,1] <- 1 / val_p3[,1]
  # params[,2] <- val_p3[,2] / abs(val_p3[,1])- val_p3[,4]
  # params[,4] <- val_p3[,4]
  # params <- matrix(0, ncol=4, nrow=6)
  # params[1,] <- c(8.07,14.7,0.26,6.62)
  # params[2,] <- c(4.07,4.12,0.04,4.50)
  # params[3,] <- c(10.1,20.4,1.90,7.02)
  # params[4,] <- c(7.73,13.5,0.15,6.38)
  # params[5,] <- c(3.11,2.20,0.03,3.56)
  # params[6,] <- c(8.17,13.6,0.27,6.01)
  # 
  exceedeances <- makeExceedances(data = data, thresholds = thres, normalize = T)

  # TODO refactor
  exceedeances_cdf_ecdf <- exceedeances
  exceedeances_cdf_ecdf <- t(apply(exceedeances_cdf_ecdf, MARGIN = 1,
        FUN = function(x){
          return(plgpd_unif_at_zero.row(xs=x, p.zeroes = p.zeroes, params.mat = params))
          }))
  exceedeances_cdf_ecdf <- as.matrix(exceedeances_cdf_ecdf)
  
  
  for(i in 1:n_vars){
    # This creates ordered uniform samples in the same order 
    # as in data.
    exceedeances_cdf_ecdf[which(exceedeances[,i]==0), i] <- 
      ecdf(data[which(exceedeances[,i]==0), i])(data[which(exceedeances[,i]==0), i]) * p.zeroes[i]
  }

  list_of_list_horizons <- list()
  cat("horizons", horizons, "\n")
  for(h in horizons){
    list_of_matrices_conditional <- list()
    # creates a square matrix with all the vars in
    quantile.update.values <- matrix(0, 
                                     nrow = n_vars, 
                                     ncol = n_vars)
    colnames(quantile.update.values) <- colnames(exceedeances)
    rownames(quantile.update.values) <- colnames(exceedeances)
    
    for(i in conditional_on){
      # creates a temporary matrix with cols equal to number of nvars + 1
      # and rows such that the i-th component is an extreme
      cat("n_samples", n_samples, " h ", h, "\n")
      mat_temp <- matrix(0,
                         nrow = length(which(exceedeances[1:(n_samples-h), i] > 0)),
                         ncol = n_vars+1)
      # filtering the i-th eCDF column with extremes 
      print(head(exceedeances_cdf_ecdf[,1]))
      temp <- exceedeances_cdf_ecdf[which(exceedeances[1:(n_samples-h), i] > 0), i]

      # addting those values in the nvars + 1 column
      mat_temp[,n_vars+1] <- ecdf(temp)(temp)
      
      # looping on the first nvars columns
      for(j in 1:n_vars){
        # filtering data of j-th vars when i-th is extreme h timesteps before
        # print(n_samples)
        # print(h)
        # print(max(which(exceedeances[1:(n_samples-h), i] > 0)+h))
        # print(min(which(exceedeances[1:(n_samples-h), i] > 0)+h))
        # print(dim(exceedeances_cdf_ecdf))
        
        data_j <- exceedeances_cdf_ecdf[which(exceedeances[1:(n_samples-h), i] > 0)+h, j]
        # computing the probability that j-th was an extreme as well 
        # h timesteps after i-th was an extreme
        quantile.update.values[i, j] <- mean(data_j <= p.zeroes[j]) # WARNING TODO <= or >=
        #hist(data_j, breaks=50)
        data_j <- ecdf(data_j)(data_j)
        # saving the unif values of j-th var used here.
        mat_temp[,j] <- data_j 
      }
      
      colnames(mat_temp) <- c(colnames(exceedeances), colnames(exceedeances)[i])
      list_of_matrices_conditional[[i]] <- mat_temp
      
      print(colSums(is.na(mat_temp)))
      
      if(sum(is.na(mat_temp))>0){
        print('NA in conditional matrices')
      }
    }
    
    list_of_list_horizons[[h]] <- list(
        unif.values = list_of_matrices_conditional,
        quantiles.values = quantile.update.values
      )
  }
  
  if(save){
    file_name <- paste(name, ".RData", sep="")
    cat(paste("Backing up the conditional matrices as", file_name, "..."))
    rlist::list.save(list_of_list_horizons,
                     file=file_name)
    cat("done\n")
  }
  return(list_of_list_horizons)
}


#' @param horizons (integer) sequence of selected extreme horizons.
#' @param list_of_matrix Collection of conditional matrices (as created by
#'   \code{makeConditionalMatrices}).
#' @param save Logical (flag, default to TRUE). Whether we save the conditonal
#'   matrices.
#' @param sparse Logical flag (default is FALSE). Whether to perform mBICV 
#' sparse vine computation. 
#' @seealso \code{makeConditionalMatrices}.
#' @examples fitExceedancesVines(threshold_data[,100:102], list_of_list_horizons)
fitExceedancesVines <- function(horizons, list_of_matrix, save=F, sparse=F){
  #list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")
  #cat("nvars?", length(list_of_matrix[[horizons[1]]]$unif.values[[1]][1,]))
  list_of_list_horizons <- list_of_matrix
  list_of_list_horizons_vines <- list()
  
  k <- 1
  while(list_of_list_horizons[[horizons[1]]]$unif.values[[k]] %>% is.null){
    k <- k + 1
  }
  
  n_vars <- length(list_of_list_horizons[[horizons[1]]]$unif.values[[k]][1,]) - 1 
  col_names <- colnames(list_of_list_horizons[[horizons[1]]]$unif.values[[k]])
 
  for(h in horizons){
    list_of_vines_mat <- list()
    cat("Horizon: ", h, "\n")
    i <- 1
    for(sub_matrix_list in list_of_matrix[[h]]$unif.values){
      if(!is.null(sub_matrix_list)){
          cat("---> fit for", col_names[i], "\n")
          time_proc <- proc.time()[3]
          list_of_vines_mat[[i]] <- fitExceedancesSingleVine(data_matrix = 
                                                               sub_matrix_list,
                                                             save = save, 
                                                             sparse = sparse)
          cat("       |-----> done in", round((proc.time()[3] - time_proc), 2), "s. \n")
      }
      i <- i + 1
    }
    
    list_of_list_horizons_vines[[h]] <- list_of_vines_mat
  }
  
  if(save){
    rlist::list.save(list_of_list_horizons_vines, file = "cond-mat-vines-12361224-v2.RData")
  }
  return(list_of_list_horizons_vines)
}


#' @param data_matrix Collection of conditional matrices (as created by
#'   \code{makeConditionalMatrices}), for one specific horizon and one
#'   conditional variable.
#' @param save Logical (flag, default to TRUE). Whether we save the conditonal
#'   matrices.
#' @param sparse Logical flag (default is FALSE). Whether to perform mBICV
#'   sparse vine computation.
#' @seealso \code{makeConditionalMatrices} and \code{fitExceedancesVines}.
#' @examples fitExceedancesVines(threshold_data[,100:102], list_of_list_horizons)
fitExceedancesSingleVine <- function(data_matrix, save=F, sparse=F){
  num_cores <- parallel::detectCores()
  if(sparse){
    result <- rvinecopulib::vinecop( # TODO warning include vinecop
      data = data_matrix,
      family_set = c("gumbel", "indep", "clayton"),  psi0 = 0.95,
      selcrit = "mbicv", trunc_lvl = Inf, tree_crit = "tau", threshold = 0,
      par_method = "mle", cores = num_cores-1)
  }else{
    result <- VineCopula::RVineStructureSelect( # TODO warning include vinecop
      data = data_matrix, familyset = c(0, 3, 4), type = 0,
      selectioncrit = "BIC", indeptest = TRUE, level = 0.05,
      trunclevel = NA, progress = FALSE, weights = NA, treecrit = "tau",
      se = FALSE, rotations = TRUE, method = "mle", cores = num_cores-1)
  }
  
  return(result)
}

#source("multi_ev.R")
#'computeTRON allows to compute TRON probabilities very easily!
#'@param data clean dataset
#'@param p.zeroes a scalar or vector (as large as the number of columns of
#'  data). of probabilities to be an exceedance of zero (proba of NOT being an
#'  extreme).
#'@param horizons a integer or vector (of integers) of look-ahead horizons for
#'  extremes.
#'@param clusters a integer or vector (of integers) of clusters size in the
#'  autocorrelation sense. See \code{\link[ev.trawl]{GenerateParameters}}.
#'@param n_samples Number of samples to compute the TRON probabilites via
#'  Monte-Carlo.
#'@param conditional_on name or index of variable to condition the extremes on.
#'  Default if \code{NA} and in this case, it creates the full collection of
#'  matrices and vines.
#'@param name_matrices_file Default is NA. If NA, we use \code{makeFileName(NA,
#'  "_matrix", ".RData")}, otherwise we use \code{name_matrices_file.RData}.
#'@param name_vine_file Default is NA. If NA, we use \code{makeFileName(NA,
#'  "_matrix", ".RData")}, otherwise we use \code{name_vine_file.RData}.
#'@param name_tron_file Default is NA. If NA, we use \code{makeFileName(NA,
#'  "_matrix", ".RData")}, otherwise we use \code{name_tron_file.RData}.
#'@param save Logical (default is TRUE) to save matrices, vines and tron
#'  probabilities as RData files.
#'@param sparse Logical flag (default is FALSE). Whether to perform mBICV 
#'  sparse vine computation. 
#'@return Returns a list of TRON probabilities with horizons as keys.
#'@seealso \code{\link[ev.trawl]{GenerateParameters}} for \code{clusters} and 
#'  paper sparse vine computation \link{https://arxiv.org/abs/1801.09739}.
computeTRON <- function(data, p.zeroes, horizons, clusters, n_samples, conditional_on=NA,
                        name_matrices_file=NA, name_vine_file=NA, 
                        name_tron_file=NA, save=TRUE, sparse=FALSE){
  name_matrices_file <- makeFileName(name_matrices_file, 
                                     tag = "_matrix",
                                     extension = ".RData")
  name_vine_file <- makeFileName(name_vine_file, 
                                 tag = "_vines",
                                 extension = ".RData")
  name_tron_file <- makeFileName(name_tron_file, 
                                 tag = "_tron",
                                 extension = ".RData")
  
  exceendances <- makeExceedances(data = data,
                                  thresholds = 
                                    getThresholds(data, p.exceed = p.zeroes),
                                  normalize = TRUE)
  
  # fit univ models + compute the matrices
  #print(horizons)
  list_of_mat <- makeConditionalMatrices(data = data,
                                         p.zeroes = p.zeroes,
                                         horizons = horizons,
                                         conditional_on = conditional_on,
                                         clusters_size = clusters,
                                         n_samples = n_samples,
                                         name = 'back_up_matrix',
                                         save = T,
                                         optim = T)
  if(save){
    rlist::list.save(list_of_mat, name_matrices_file) # save
  }
  #
  # # # compute the vines
  print("Fitting vines:")
  list_vines <- fitExceedancesVines(horizons = horizons,
                                    list_of_matrix = list_of_mat,
                                    save = F,
                                    sparse = sparse)
  if(save){
    rlist::list.save(list_vines, name_vine_file) #save
  }


  # list_of_mat <- rlist::list.load("matrix_east_light_1_to_72_matrix.RData")
  # list_vines <- rlist::list.load("vine_east_light_1_to_72_vines.RData")
  # # compute TRON
  print("Computing TRONs:")
  #list_of_mat <- rlist::list.load("matrix_east_light_1_to_72_2nd_matrix.RData")
  #list_vines <- rlist::list.load("vine_east_light_1_to_72_2nd_vines.RData")
  tron <- computeTRONwithLists(data = data,
                               horizons = horizons,
                               list_vines = list_vines,
                               list_of_matrix = list_of_mat,
                               sparse = sparse)
  if(save){
    rlist::list.save(tron, name_tron_file) # save
  }
  return(tron)
}


#' @param list_vines Fitted (single) vine as created by
#'   \code{fitExceedancesVines}.
#' @param list_of_matrix Collection of conditional matrices (as created by
#'   \code{makeConditionalMatrices}).
#' @param N number of samples to estimate TRON probabilities.
#' @param sparse Logical flag (default is FALSE). Whether to perform mBICV
#'   sparse vine computation.
#' @seealso \code{makeConditionalMatrices} and \code{fitExceedancesVines}.
#' @examples TODO WARNING
computeTRONwithListSingle <- function(vine, quantile_values, N, sparse){
  if(sparse){
    vine_sim <- rvinecopulib::rvinecop(vine = vine, 
                                    n = N,  
                                    cores = parallel::detectCores()-1)
  }else{
    vine_sim <- VineCopula::RVineSim(RVM = vine,
                                     N = N)
  }
  
  vine_sim <- vine_sim[,1:(length(vine_sim[1,])-1)]
  vine_sim <- t(apply(vine_sim, MARGIN = 1, FUN = function(x){x>quantile_values}))
  
  # exporting mean and sd
  results <- list()
  results[["mean"]] <- t(apply(vine_sim, MARGIN = 2, FUN = mean, na.rm=TRUE))
  results[["sd"]] <- t(apply(vine_sim, MARGIN = 2, FUN = sd, na.rm=TRUE))#/sqrt(length(which(!is.na(vine_sim[,1]))) - 1)
  
  return(results)
}


#' @param data dataset.
#' @param horizons (integer) sequence of selected extreme horizons.
#' @param list_vines Collection of vines as created by
#'   \code{fitExceedancesVines}.
#' @param list_of_matrix Collection of conditional matrices (as created by
#'   \code{makeConditionalMatrices}).
#' @param N number of samples to estimate TRON probabilities.
#' @param save Logical (flag, default to TRUE). Whether we save the conditonal
#'   matrices.
#' @param sparse Logical flag (default is FALSE). Whether to perform mBICV
#'   sparse vine computation.
#' @seealso \code{makeConditionalMatrices} and \code{fitExceedancesVines}.
#' @examples TODO WARNING
computeTRONwithLists <- function(data, horizons, list_vines, list_of_matrix, N=100000, save=F, sparse=F){
  tron_probabilities <- list()
  set.seed(42)
  n_vars <- length(data[1,])
  for(h in horizons){
    tron_proba_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
    colnames(tron_proba_matrix) <- colnames(data)
    rownames(tron_proba_matrix) <- colnames(data)
    tron_proba_matrix_sd <- tron_proba_matrix
    cat(paste("Horizon", h, "\n"))
    #for(i in 1:n_vars){
    i <- 1
    actual_conditional <- c()
    for(current_vine in list_vines[[h]]){
      if(!is.null(current_vine)){
        actual_conditional <- c(actual_conditional, i)
      #current_vine <- list_vines[[h]][[i]]
        current_quantiles <- list_of_matrix[[h]]$quantiles.values[i,]
        cat(paste("--> extreme in", colnames(tron_proba_matrix)[i]), "...")
        
        vine_sim_statistics <- 
          computeTRONwithListSingle(
                                  vine = current_vine,
                                  quantile_values = current_quantiles,
                                  N = N,
                                  sparse = sparse)
        tron_proba_matrix[i,] <- vine_sim_statistics$mean
        tron_proba_matrix_sd[i,] <- vine_sim_statistics$sd
        # print(tron_proba_matrix[i,])
        # print(tron_proba_matrix_sd[i,])
        cat("\t done\n")
      }
      i <- i + 1
    }
    tron_probabilities[[h]] <- list(mean=tron_proba_matrix[actual_conditional,], sd=tron_proba_matrix_sd[actual_conditional,])
  }
  # tron_probabilities[[1]]$mean
  # tron_probabilities[[1]]$sd
  if(save){
    list.save(tron_probabilities, file = "tron-cond-test.RData")
  }
  return(tron_probabilities)
}


#' @example XSRKtoABRK(params)
XSRKtoABRK <- function(params){
  # params should be (xi, sigma, rho, kappa)
  if(length(params) != 4) stop('In', match.call[1], 'Wrong dims: should be 4.')
  abrk <- params
  abrk[1] <- 1.0 / abrk[1]
  abrk[2] <- params[2]/abs(params[1]) - params[4]
  
  return(abrk)
}

#' @example XSRKtoABRK(params)
ABRKtoXSRK <- function(params){
  # params should be (alpha, beta, rho, kappa)
  if(length(params) != 4) stop('In', match.call[1], 'Wrong dims: should be 4.')
  abrk <- params
  abrk[1] <- 1.0 / abrk[1]
  abrk[2] <- (params[2] + params[4])/abs(params[1])
  
  return(abrk)
}
