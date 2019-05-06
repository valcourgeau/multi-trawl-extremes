library(parallel)

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)

a <- parLapply(cl, 2:100000,
          function(exponent)
            exponent)

stopCluster(cl)

pairwise_likelihood_single_full_parallel <- function(times, values, alpha, beta, kappa, rho, delta, logscale=T, transformation=F){
  ok_ind <- which(!is.na(values))
  values <- values[ok_ind]
  times <- times[ok_ind]
  k <- length(values)
  
  temp <- 0.0
  upper <- pmin(1:(k-1)+delta, k)
  lower <- 2:k
  
  accepted <- 0
  total <- 0
  
  # Constructing time_val
  job_list <- pmin(rep(2:(delta), length(times)-1)+rep(0:(length(times)-2), each=delta-1), length(times))
  job_list <- rbind(rep(1:(length(times)-1), each=delta-1), job_list)
  job_list <- rbind(matrix(times[job_list], nrow=delta-1), matrix(values[job_list], nrow=delta-1))
  job_list <- job_list[,1:(length(job_list[1,])-delta+2)]
  job_list <- as.list(as.data.frame(job_list))
  
  no_cores <- 5 - 1
  cl <- makeCluster(no_cores, type="PSOCK")
  clusterExport(cl,c("pairwise_likelihood_single_pair"))
  job_results <- parLapply(cl, job_list,
                 function(job_params){
                   return(pairwise_likelihood_single_pair(t1 = as.numeric(job_params[[1]][1]),
                                                          x1 = as.numeric(job_params[[1]][3]), 
                                                          t2 = as.numeric(job_params[[1]][2]),
                                                          x2 = as.numeric(job_params[[1]][4]),
                                                          alpha = as.numeric(alpha),
                                                          beta = as.numeric(beta),
                                                          kappa = as.numeric(kappa),
                                                          rho = as.numeric(rho),
                                                          transformation = as.logical(transformation))) # TODO add transformation
                 }
                 )
  
  base <- 4
  test <- function (exponent) {
    foreach(exponent = 2:4, 
            .combine = c,
            .export = "base")  %dopar%  
        pairwise_likelihood_single_pair(t1 = job_params[1],
                                         x1 = job_params[3], 
                                         t2 = job_params[2],
                                         x2 = job_params[4],
                                         alpha = alpha,
                                         beta = beta,
                                         kappa = kappa,
                                         rho = rho,
                                         transformation = transformation) # TODO add transformation
  }
  test()
  foreach(exponent = 2:4, 
          .combine = c)  %dopar%  
    base^exponent
  stopImplicitCluster()
  stopCluster(cl)
  #accepted <- length(which(job_results > 0))
  
  job_results[which(job_results > 0 & job_results < 1e-12)] <- -1000
  temp <- sum(log(job_results))
  
  #cat("Accepted: ", accepted/total, "\n")
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

set.seed(42)
times <- 1:10
values <- runif(length(times))
alpha <- 6
beta <- 10
kappa <- 0.3
rho <- 0.4
delta <- 3
transformation <- F

pairwise_likelihood_single_full_parallel(times = times,
                                         values = values,
                                         alpha = alpha,
                                         beta = beta,
                                         kappa = kappa,
                                         rho = rho,
                                         delta = delta)

