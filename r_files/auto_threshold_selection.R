## load the eva package
if (!require ("eva") ) {
  install.packages ("eva")
  library(eva)
}
## Find source function for the method of Northrop and Coleman (2014) here :
## www.homepages .ucl .ac.uk/~ ucakpjn / SOFTWARE / NorthropColeman2014.zip
## After unzipping , the functions can be sourced in
source("NorthropColeman2014.fns")
## Moran ’s test ( outputs p- value )
gpd.morans <- function (data) {
  m <- length (data) + 1
  fit <- gpdFit(data, threshold = NA, nextremes = length (data), method = "mps")
  mu_m <- m*( log (m) - digamma (1) ) - (1 / 2) - (1 / (12 *m) )
  sigma_m <- sqrt( m* (( pi ^2 / 6) - 1) - (1 / 2) - (1 / (6 *m )))
  C1 <- mu_m - ( sqrt ( m/ 2) ) * sigma_m
  C2 <- sigma_m / sqrt (2 * m )
  stat <- ( fit $ moran + 1 - C1 ) / C2
  1 - as.numeric(pchisq(stat, length (data)))
}


# WARNING: dat + gpdAd !- score.fitrange

threshold_test <- function(dat, p.zero){
  result <- list()
  #result <- rep(0, 4)
  ## If shape outside range of critical values table , use bootstrap
  result[["Ad"]] <- tryCatch(gpdAd(dat[dat>quantile(dat, p.zero)], allowParallel = T, numCores = 7)$p.value , error = function(e) {
    print("error in Ad")
    tryCatch ( gpdAd(dat , bootstrap = TRUE , bootnum = 1000, allowParallel = T, numCores = 7)$p.value , error = function (e) NA)
  })
  result[["Cvm"]] <- tryCatch(gpdCvm(dat[dat>quantile(dat, p.zero)], allowParallel = T, numCores = 7)$p.value , error = function(e) {
    print("error in Cvm")
    tryCatch ( gpdCvm(dat , bootstrap = TRUE , bootnum = 1000, allowParallel = T, numCores = 7)$p.value , error = function (e) NA
    ) })
  ## Need to give Northrop - Coleman ’s score test a set of thresholds
  proba <- seq(from=0.8, to=0.99, by=0.02)
  threshes <- quantile(dat, probs = proba)
  
  ## Sometimes score test fails
  temp <- tryCatch(score.fitrange(dat , threshes)$e.p.values , error = function ( e) NA)
  result[["score"]] <- list(proba=proba[2:length(proba)], p.values=temp)
  
  temp <- tryCatch(gpd.morans(dat[dat > dat_threholds]) , error = function (e ) NA )
  result[["morans"]] <- list(proba=proba[2:length(proba)], p.values=temp)
  return(result)
}
# 
# threshold_test(pmax(0,ts(core_energy_data[,100])-quantile(ts(core_energy_data[,100]), 0.98)), 0.0)
# 
# 
# threshold_test(solcl[,3], thr_stl.sol[1])
# abline(v=thr_stl.sol[1], col="red")
# 
# threshold_test(solcl[,5], thr_stl.sol[3])
# abline(v=thr_stl.sol[3], col="red")
# 
# threshold_test(esol[,4][esol[,4] > 0])
# threshold_test(solcl[,3], thr_stl.sol[1])
# 
# ## Function to generate GPD mixture
# rgpd_mix <- function (n , loc = 0, scale = 1 , shape1 , shape2 ) {
#   c( rgpd ( trunc (n/ 2) , loc = loc , scale = scale , shape = shape1 ) ,
#      rgpd ( trunc (n / 2) , loc = loc , scale = scale , shape = shape2 ))
# }
# 
# ## Simulation settings
# n <- c(50 , 100 , 200 , 400)
# nsim <- 1000
# set.seed (7)
# ## sim _fun outputs the power and percentage of failed simulations
# ## gen _fun is the data generating function , followed by required args
# weibull_75 <- sim_fun ( nsim = nsim , n=n , gen_fun = rweibull , shape =0.75)
# weibull_125 <- sim_fun ( nsim = nsim , n =n , gen_fun = rweibull , shape =1.25)
# lognorm_0_1 <- sim_fun ( nsim = nsim , n =n , gen_fun = rlnorm , meanlog =0 , sdlog =1)
# gamma_2_1 <- sim_fun ( nsim = nsim , n=n , gen_fun = rgamma , shape =2 , scale =1)
# gpd_mix1 <- sim_fun ( nsim = nsim , n =n , gen_fun = rgpd_mix , shape1 =0.4 , shape2 = -0.4)
# gpd_mix2 <- sim_fun ( nsim = nsim , n =n , gen_fun = rgpd_mix , shape1 =0.4 , shape2 =0)
# gpd_mix3 <- sim_fun ( nsim = nsim , n =n , gen_fun = rgpd_mix , shape1 =0.25 , shape2 = -0.25)
# gpd_size <- sim_fun ( nsim = nsim , n =n , gen_fun = rgpd , shape =0.25)
