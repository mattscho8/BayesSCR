library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)

### get the blackbear data ###
tmp = readRDS('../data.rds') # gets data summaries that we use in all models
list2env(tmp, .GlobalEnv)

#### Plot
plot(grid, xlim = xmask, ylim = ymask, pch = 4, lwd = 2)

## other variables
M = 150

yfull = rbind(xobs,matrix(0,M-n,K))

##### constants
# x.max, y.max are same
# det.xy is grid
# K is K
# M is M
# J is J

##### data
# y is yfull

###### inits
# psi
# p0
# sigma
# z
# sxy
inits = function(){
  return(list(psi = runif(1, 0.25, 0.75),
             p0 = runif(1, 0.25, 0.75),
             sigma = runif(1, 0.25, 1),
             sx = runif(M, 0,x.max),
             sy = runif(M,0,y.max)))
}
dat = list(x_max = x.max, y_max = y.max, det_x = grid[,1], det_y = grid[,2], K = K, 
           M = M, J = J, xobs = xobs, n = n)


aaa = "

data {
  int<lower=0> M; // upper bound
  int<lower=0> n; // the number of observed individuals
  int<lower=0> K; // the number of traps
  int<lower=0> J; // the number of sampling occasions
  real<lower=0> x_max;
  real<lower=0> y_max;
  row_vector[K] det_x; // x-coordinates for traps
  row_vector[K] det_y;
  array[n] row_vector[K] xobs; // observed number of captures at each detector
}

parameters {
  real<lower=0, upper=1> p0; 
  real<lower=0, upper = 50> sigma;
  real<lower=0, upper=1> psi;
  row_vector<lower=0,upper=x_max>[M] sx;
  row_vector<lower=0,upper=y_max>[M] sy;
}


model {
  array[M] row_vector[K] logp;
  array[M] row_vector[K] d2;
  
  for(i in 1:M){
    d2[i] = (sx[i] - det_x).*(sx[i] - det_x) + (sy[i] - det_y).*(sy[i] - det_y);
    logp[i] = log(p0) - d2[i]/(2*sigma^2);
  }
  
  for(i in 1:n){
    target += log(psi);
    target += xobs[i].*logp[i] + (J-xobs[i]).*log1m_exp(logp[i]); 
  }
  for(i in (n+1):M){
    target += log_mix(psi, J*sum(log1m_exp(logp[i])), 0);
  }
      
  // uniform priors as in NIMBLE
}

"

cat(aaa, file = "mod_da.stan")

stan_fit = cmdstan_model("mod_da.stan")

fit <- stan_fit$sample(
  data = dat,
  init = inits,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 50
)

parout = as_draws_df(fit$draws(c("p0", "sigma","psi")))
mcmc_trace(parout)
summary(parout)
