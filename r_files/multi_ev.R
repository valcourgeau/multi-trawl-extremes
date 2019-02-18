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


#' Wrapper from \code{ev.trawl} GenerateParameters to fit univariate latent
#' trawl extreme values model.
#'
#' @param data cleaned dataset.
#' @param clusters_size Extreme clusters size (integer) vector .
#' @param thresholds Thresholds values for extremes.
#' @param optim Logical (default=TRUE). Whether to execute BFGS optimisation with fixed alpha.
#' @return Set of thresholds values as large as \code{data} after which values
#'   in data are considered extremes.
findUnivariateParams <- function(data, clusters_size, thresholds, optim=T){
  n_vars <- length(data[1,])
  n_rows <- length(data[,1])
  
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
      val_params[i_agent,] <- ev.trawl::GenerateParameters(data = exc[,i_agent],
                                                           cluster.size = clusters_size[i_agent])
      cat(" done.\n")
      if(optim){
        cat("|---> Preparing optimisation...")
        marginal_values <- exc[1:2000,i_agent] # TODO 1000!
        marginal_times <- (1:n_rows)/n_rows
        exp_params_names <- c("alpha", "beta", "rho", "kappa")
        fixed_params_names <- c("alpha")
        fixed_params_index <- which(exp_params_names %in% fixed_params_names)
        
        fn_to_optim <- function(x) {
          return(-UnivariateFullPL(values = marginal_values, 
                                   params = x, 
                                   times = marginal_times,
                                   delta = clusters_size[i_agent], 
                                   model_vars_names = exp_params_names,
                                   fixed_names = fixed_params_names,
                                   fixed_params = x[fixed_params_index],
                                   transformation = T))
        }
        lower <-c(1,
                  0.01,
                  0.1)
        upper <-c(100,
                  2,
                  100)
        
        initial_guess <- list()
        for(j in 1:length(val_params[1,])){
          if(!(j %in% fixed_params_index)){
            initial_guess[[exp_params_names[j]]] <- val_params[i_agent,j] 
          }
        }
        cat(" initialise...")
        time_to_cv <- proc.time()[3]
        res <- stats::optim(fn_to_optim, 
                            par = initial_guess,
                            method = "L-BFGS-B",
                            lower = lower, upper = upper,
                            control = list(maxit=100, 
                                           factr=1e13, 
                                           trace=0)
                      )$par
        cat("done in ", proc.time()[3]-time_to_cv, "s.\n", sep = "")
        res_concat <- rep(0, 4)
        res_concat[fixed_params_index] <- val_params[i_agent, fixed_params_index]
        res_concat[-fixed_params_index] <- res
        val_params[i_agent, ] <- res_concat
      }
      
  }
  return(val_params)
}

source("infer_latent_value.R")

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

getThresholds(core_energy_data[,1:3], c(0.1,0.5,0.8))

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

#' Constructs a collecttion of matrices with the primary set of keys being the
#' horizons and the second one being the column indices. Those matrices are the
#' data conditional on extreme value in the respectivel column.
#'
#' @param data clean dataset.
#' @param p.zeroes threshold in probability to be considered an extreme (ex:
#'   0.95).
#' @param horizons (integer) sequence of selected extreme horizons.
#' @param clusters_size Extremes cluster size for each marginal.
#' @param save Logical (flag, default to TRUE). Whether we save the conditonal
#'   matrices.
#' @param n_samples Number of rows selected to create the matrices (from 1 to
#'   \code{n_samples}).
#' @param name Filename if the collection of matrices is saved.
#' @return Double-leveled list with horizons as first set of keys and columns
#'   indices as second set of keys. E.g. \code{result[[2]][[4]]} would pick the
#'   second horizon on the 4th variable.
makeConditionalMatrices <- function(data, p.zeroes,
                                    horizons, clusters_size, n_samples = length(data[,1]), save=T, 
                                    name="conditional-matrices"){
  # adapt the size of p.zeroes to the number of cols
  if(length(p.zeroes) == 1){
    p.zeroes <- rep(p.zeroes, length(data[1,]))
  }
  
  #p.zeroes <- computePZero(params)
  thres <- getThresholds(data = data, p.exceed = p.zeroes)
  
  params <- findUnivariateParams(data = data, clusters_size = clusters_size, thresholds = thres)
  exceedeances <- makeExceedances(data = data, thresholds = thres, normalize = T)

  # TODO refactor
  exceedeances_cdf_ecdf <- exceedeances
  exceedeances_cdf_ecdf <- t(apply(exceedeances_cdf_ecdf, MARGIN = 1,
        FUN = function(x){
          return(plgpd_unif_at_zero.row(xs=x, p.zeroes = p.zeroes, params.mat = params))
          }))
  exceedeances_cdf_ecdf <- as.matrix(exceedeances_cdf_ecdf)
  
  
  for(i in 1:length(exceedeances[1,])){
    # This creates ordered uniform samples in the same order 
    # as in data.
    exceedeances_cdf_ecdf[which(exceedeances[,i]==0), i] <- 
      ecdf(data[which(exceedeances[,i]==0), i])(data[which(exceedeances[,i]==0), i]) * p.zeroes[i]
  }

  list_of_list_horizons <- list()
  n_vars <- length(data[1,])
  for(h in horizons){
    list_of_matrices_conditional <- list()
    # creates a square matrix with all the vars in
    quantile.update.values <- matrix(0, 
                                     nrow = length(exceedeances[1,]), 
                                     ncol = length(exceedeances[1,]))
    colnames(quantile.update.values) <- colnames(exceedeances)
    rownames(quantile.update.values) <- colnames(exceedeances)
    
    for(i in 1:n_vars){
      # creates a temporary matrix with cols equal to number of nvars + 1
      # and rows such that the i-th component is an extreme
      mat_temp <- matrix(0,
                         nrow = length(which(exceedeances[1:(n_samples-h), i] > 0)),
                         ncol = n_vars+1)
      # filtering the i-th eCDF column with extremes 
      temp <- exceedeances_cdf_ecdf[which(exceedeances[1:(n_samples-h), i] > 0), i]
      
      # addting those values in the nvars + 1 column
      mat_temp[,n_vars+1] <- ecdf(temp)(temp)
      
      # looping on the first nvars columns
      for(j in 1:n_vars){
        # filtering data of j-th vars when i-th is extreme h timesteps before
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
    }
    
    list_of_list_horizons[[h]] <- list(
        unif.values = list_of_matrices_conditional,
        quantiles.values = quantile.update.values
      )
  }
  
  if(save){
    file_name <- paste(name, ".RData", sep="")
    cat(paste("Backing up the conditional matrices as", file_name, "..."))
    list.save(list_of_list_horizons,
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
#' @seealso \code{makeConditionalMatrices}.
#' @examples fitExceedancesVines(threshold_data[,100:102], list_of_list_horizons)
fitExceedancesVines <- function(horizons, list_of_matrix, save=F){
  #list_of_list_horizons <- list.load(file = "conditional-mat-test.RData")
  paste("nvars?", length(list_of_matrix[[horizons[1]]]$unif.values[[1]][1,])) %>% print
  list_of_list_horizons <- list_of_matrix
  list_of_list_horizons_vines <- list()
  
  n_vars <- length(list_of_list_horizons[[horizons[1]]]$unif.values[[1]][1,]) - 1 
  col_names <- colnames(list_of_list_horizons[[horizons[1]]]$unif.values[[1]])
  
  for(h in horizons){
    list_of_vines_mat <- list()
    cat("Horizon: ", h, "\n")
    for(i in 1:n_vars){
      cat("--->", col_names[i], "\n")
      time_proc <- proc.time()[3]
      list_of_vines_mat[[i]] <- RVineStructureSelect(
        data = list_of_list_horizons[[h]]$unif.values[[i]], familyset = c(3, 4), type = 0,
        selectioncrit = "AIC", indeptest = TRUE, level = 0.05,
        trunclevel = NA, progress = FALSE, weights = NA, treecrit = "tau",
        se = FALSE, rotations = TRUE, method = "mle", cores = 7)
      cat("       |-----> done in", round((proc.time()[3] - time_proc), 2), "s. \n")
    }
    list_of_list_horizons_vines[[h]] <- list_of_vines_mat
  }
  
  if(save){
    list.save(list_of_list_horizons_vines, file = "cond-mat-vines-12361224-v2.RData")
  }
  return(list_of_list_horizons_vines)
}

