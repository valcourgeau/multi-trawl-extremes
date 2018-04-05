source('gaver_lewis_simulation.R')
source('warren_simulation.R')

n_sims <- 1000
n_timesteps <- 4400

# NO TRF
## Example 1 alpha = 4, beta = 4, kappa: 10% exceedances, rho = 0.2
alpha <- 4
beta <- 4
kappa <- beta * (0.1^{-1/alpha} - 1)
rho <- 0.2
### Warren process
set.seed(42)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rwprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "warren_4_4_3.113118_0.2.csv", row.names = F, col.names = F)

### Gaver-Lewis process
set.seed(42)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rlgprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "gaverlewis_4_4_3.113118_0.2.csv", row.names = F, col.names = F)

### Latent Noven
## Trawl process simulation
library(gPdtest)
times <- 1:n_timesteps
### Generating the functions
trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)

set.seed(42)
exceeds <- matrix(0, nrow = length(times), ncol = n_sims)
for(i in 1:n_sims){
  exceeds[,i] <- rlexceed(alpha = alpha,
                      beta = beta,
                      kappa = kappa,
                      trawl_fs = trawl_1,
                      trawl_fs_prim = trawl_1_prim,
                      times = times,
                      n = 1,
                      transformation = F)
  if(i %% 10 == 0) print(i)
}
# length(which(exceeds[,2]>0.0))/length(exceeds[,1])
# plot(density(exceeds[,2]))
# hist(exceeds[,1], probability = T)
# lines(seq(0,3,length.out = 500), dgamma(seq(0,3,length.out = 500), shape = alpha, rate = beta))
# 
# hist(exceeds[,1][exceeds[,1] > 0.1], probability = T)
# lines(seq(0,30,length.out = 500), dgpd(seq(0,30,length.out = 500), xi = alpha, beta = (beta+kappa)/alpha))

## Example 2 alpha = 9, beta = 1, kappa: 10% exceedances, rho = 0.2
alpha <- 9
beta <- 1
kappa <- beta * (0.1^{-1/alpha} - 1)
rho <- 0.2
### Warren process
set.seed(42)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rwprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "warren_9_1_0.2915497_0.2.csv", row.names = F, col.names = F)

### Gaver-Lewis process
set.seed(42)
exceeds <- matrix(0, nrow = n_timesteps, ncol = n_sims)
for(sims in 1:n_sims){
  exceeds[, sims] <- rlgprocess(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = n_timesteps, n = 1)
}

write.table(exceeds, file =  "gaverlewis_9_1_0.2915497_0.2.csv", row.names = F, col.names = F)
