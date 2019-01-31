source('pairwise_latent_trawl.R')
source('gaver_lewis_simulation.R')
source('warren_simulation.R')
source('latent_novel_simulation.R')

n_sims <- 1000
n_timesteps <- 4400

# NO TRF
## Example 1 alpha = 4, beta = 4, kappa: 10% exceedances, rho = 0.2
alpha <- 4
beta <- 4
kappa <- beta * (0.1^{-1/alpha} - 1)
rho <- 0.2
### Warren process
set.seed(40)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rwprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "warren_4_4_3.113118_0.2_1000.csv", row.names = F, col.names = F)

### Gaver-Lewis process
set.seed(40)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rlgprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "gaverlewis_4_4_3.113118_0.2_1000.csv", row.names = F, col.names = F)

### Latent Noven
## Trawl process simulation
library(gPdtest)
times <- 1:4430
f_times <- 4400

alpha <- 4
beta <- 4
rho <- 0.2
kappa <- 3.113118

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
exceeds <- matrix(0, nrow = f_times, ncol = n_sims) #  TODO 30 to n_sims
for(i in 1:n_sims){ # TODO 30 again...
  exceeds[,i] <- rlexceed(alpha = alpha,
                      beta = beta,
                      kappa = kappa,
                      trawl_fs = trawl_1,
                      trawl_fs_prim = trawl_1_prim,
                      times = times,
                      n = 1,
                      transformation = F)
  print(i)
}

write.table(exceeds, file =  "latent_4_4_3.113118_0.2_1000_v3.csv", row.names = F, col.names = F) 

data_test <- exceeds[,2][exceeds[,2]>0]
data_test
gpdFit(exceeds[,5], u=0)
length(data_test)/length(times)


# length(which(exceeds[,2]>0.0))/length(exceeds[,1])
plot(density(exceeds[,1][exceeds[,1] < 2]))
hist(exceeds[,1][exceeds[,1] < 2], probability = T, breaks= 30)
lines(seq(0,3,length.out = 500), dgamma(seq(0,3,length.out = 500), shape = alpha, rate = beta))
acf(exceeds[,1])
lines(0:20, 1/exp(0:20 * rho))

hist(exceeds[,1][exceeds[,1] > 1e-6], probability = T, breaks = 100)
lines(seq(0.1,20,length.out = 500), dgpd(seq(0.1,20,length.out = 500), xi = 1/alpha, beta = (beta+kappa)/alpha))

## Example 2 alpha = 9, beta = 1, kappa: 10% exceedances, rho = 0.05
alpha <- 9
beta <- 1
rho <- 0.05
kappa <- 0.75

### Warren process
set.seed(40)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rwprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "warren_9_1_0.2915497_0.05_1000.csv", row.names = F, col.names = F)

### Gaver-Lewis process
set.seed(40)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rlgprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "gaverlewis_9_1_0.2915497_0.05_1000.csv", row.names = F, col.names = F)

### Latent Noven
## Trawl process simulation
library(gPdtest)
times <- 1:(n_timesteps+30)
alpha <- 4
beta <- 1
rho <- 0.05
kappa <- 0.75

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
n_sims <- 500
exceeds <- matrix(0, nrow = length(times)-30, ncol = n_sims)
for(i in 1:n_sims){
  exceeds[,i] <- rlexceed(alpha = alpha,
                          beta = beta,
                          kappa = kappa,
                          trawl_fs = trawl_1,
                          trawl_fs_prim = trawl_1_prim,
                          times = times,
                          n = 1,
                          transformation = F)
  print(i)
  print(length(which(exceeds[,i] > 0))/n_timesteps)
}


write.table(exceeds, file =  "latent_9_1_0.75_0.05_500.csv", row.names = F, col.names = F) 


### TRF 
library(gPdtest)
times <- 1:(n_timesteps+30)
alpha <- -4
beta <- 4
rho <- 0.2
kappa <- 9

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
n_sims <- 500
exceeds <- matrix(0, nrow = length(times)-30, ncol = n_sims)
for(i in 1:n_sims){ # TODO put nsims back in place
  exceeds[,i] <- rlexceed(alpha = alpha,
                          beta = beta,
                          kappa = kappa,
                          trawl_fs = trawl_1,
                          trawl_fs_prim = trawl_1_prim,
                          times = times,
                          n = 1,
                          n_moments = 0,
                          transformation = T)
  print(i)
  #print(length(which(exceeds[,i] > 0))/n_timesteps)
}

#trawls
hist(exceeds[,1][exceeds[,1] < 10], probability = T, breaks = 40)
lines(1:40/10, dgamma(1:40/10, shape = 1, rate = 1))

hist(exceeds[,1][exceeds[,1]>0 & exceeds[,1]<50], probability = T)
lines(1:100/10, dgpd(x = 1:100/10, xi = 1/alpha, beta = (beta+kappa)/abs(alpha)))
evir::gpd(exceeds[,1][exceeds[,1]>0], threshold = 0.0)$par.ests

write.table(exceeds, file =  "latent_minus_4_4_9_0.2_500.csv", row.names = F, col.names = F, sep = ",") 


## TRF 2
times <- 1:(n_timesteps+30)
alpha <- -4
beta <- 1
rho <- 0.10
kappa <- 15

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
n_sims <- 500
exceeds <- matrix(0, nrow = length(times)-30, ncol = n_sims)
for(i in 1:n_sims){ # TODO put nsims back in place
  exceeds[,i] <- rlexceed(alpha = alpha,
                          beta = beta,
                          kappa = kappa,
                          trawl_fs = trawl_1,
                          trawl_fs_prim = trawl_1_prim,
                          times = times,
                          n = 1,
                          n_moments = 0,
                          transformation = T)
  print(i)
  print(length(which(exceeds[,i] > 0))/n_timesteps)
}

#trawls
hist(exceeds[,1][exceeds[,1] < 10], probability = T, breaks = 40)
lines(1:40/10, dgamma(1:40/10, shape = 1, rate = 1))

hist(exceeds[,1][exceeds[,1]>0 & exceeds[,1]<200], probability = T)
lines(1:1000/10, dgpd(x = 1:1000/10, xi = 1/alpha, beta = (beta+kappa)/abs(alpha)))
evir::gpd(exceeds[,1][exceeds[,1]>0], threshold = 0.0)$par.ests

write.table(exceeds, file =  "latent_minus_4_1_15_0.1_500.csv", row.names = F, col.names = F, sep = ",") 


## TRF Example UNIV Solar Oklahoma
n_timesteps <- 5000
times <- 1:(n_timesteps+30)
alpha <- -5.5549
beta <- 0.603
rho <- 0.165
kappa <- 0.5047
n_sims <- 2

### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(40)
n_sims <- 1000
exceeds <- matrix(0, nrow = length(times)-30, ncol = n_sims)
for(i in 1:n_sims){ # TODO put nsims back in place
  exceeds[,i] <- rlexceed(alpha = alpha,
                          beta = beta,
                          kappa = kappa,
                          trawl_fs = trawl_1,
                          trawl_fs_prim = trawl_1_prim,
                          times = times,
                          n = 1,
                          n_moments = 0,
                          transformation = T)
  print(i)
  print(length(which(exceeds[,i] > 0))/n_timesteps)
}

#trawls

hist(exceeds[,1][exceeds[,1]>0 & exceeds[,1]<200], probability = T)
lines(1:1000/10, dgpd(x = 1:1000/10, xi = 1/alpha, beta = (beta+kappa)/abs(alpha)))
evir::gpd(exceeds[,1][exceeds[,1]>0], threshold = 0.0)$par.ests

write.table(exceeds, file =  "latent_univ_solar.csv", row.names = F, col.names = F, sep = ",") 
