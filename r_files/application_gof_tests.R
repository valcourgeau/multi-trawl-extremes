library(rlist)
library(VineCopula)
library(magrittr)

setwd('~/GitHub/multi-trawl-extremes/r_files/')
tron_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/r_files/2019-5-4-10-48-17_tron.RData')
vines_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/r_files/2019-5-4-10-48-17_vines.RData')
data_ap <- read.table('air_pollution_data')
data_ap_matrix <- rlist::list.load('air_pollution_rerun_3_matrix.RData')
  

tron_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/data/merged-datasets/tron_east_ligh_1_to_72_2nd_tron.RData')
vines_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/data/merged-datasets/vine_east_light_1_to_72_2nd_vines.RData')
data_ap <- read.table('air_pollution_data')
data_ap_matrix <- rlist::list.load('~/GitHub/multi-trawl-extremes/data/merged-datasets/matrix_east_light_1_to_72_2nd_matrix.RData')

horizon_set <- c(1,2,3,4,5,6,12,24,48,72)
gof_test_res <- matrix(nrow=6, ncol=2)
l_test_data <- 1
n.boostrap <- 200
set.seed(43)

gof_test_set <- c('CvM', 'KS')

GoFTesting <- function(data, matrices, horizons, n_random, n_bootstrap){stop('Not yet implemented')}

res_gof_energy <- list()
for(h in horizon_set){
  vine_set <- vines_ap[[h]]
  matrix_set <- data_ap_matrix[[h]]$unif.values
  n_vars <- ncol(matrix_set[[horizon_set[1]]]) - 1
  print('===== GoF tests  =====')
  cat('Horizon', h, '\n')
  res_gof_energy[[h]] <- list()
  for(vine_index in 1:n_vars){
    curr_vine <- vine_set[[vine_index]]
    curr_data <- matrix_set[[vine_index]]
    curr_rvm <- VineCopula::RVineMatrix(Matrix = curr_vine$Matrix,
                                        family = curr_vine$family,
                                        par = curr_vine$par,
                                        par2 = curr_vine$par2,
                                        names = curr_vine$names
    )
    index_rndm <- sample(x = 1:nrow(curr_data), replace = F, size = l_test_data)
    assertthat::are_equal(VineCopula::RVineMatrixCheck(M = curr_vine$Matrix), 1)
    cat('\t *', colnames(curr_data)[n_vars+1], 'with statistics', gof_test_set, '...')
    tmp_time <- proc.time()
    gof_white <- vapply(gof_test_set, function(test_stat){
      res_test <- VineCopula::RVineGofTest(RVM = curr_rvm, 
                                           data = curr_data[index_rndm,],
                                           B = n.boostrap,
                                           method='White',
                                           statistic = test_stat)
                    return(c(res_test[[1]]))#, res_test[[2]]))
                  },
                  rep(0,1))
    
    res_gof_energy[[h]][[vine_index]] <- gof_white
    cat('done in', (proc.time()-tmp_time)[3], 's.\n')
  }
}

res_got
res_gof_energy

