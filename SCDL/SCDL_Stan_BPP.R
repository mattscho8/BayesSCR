library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)

### get the blackbear data ###
tmp = readRDS('../data.rds') # gets data summaries that we use in all models
tmp2 = readRDS('../int_data.rds')
list2env(tmp, .GlobalEnv)
list2env(tmp2, .GlobalEnv)

#### Plot
plot(grid, xlim = xmask, ylim = ymask, pch = 4, lwd = 2)
points(int_grid, pch = 16, col = "blue", cex = 0.5)

##### constants
# x.max, y.max are same
# det.xy is grid
# K is K
# M is M
# J is J

##### data
# y is yfull

M = 150

###### inits
# psi
# p0
# sigma
# z
# sxy
inits = function(){
  return(list(p0 = runif(1),
             sigma = rlnorm(1),
             sx = runif(n,0,x.max),
             sy = runif(n,0,y.max),
             N = runif(1,n,150)))
}
dat = list(x_max = x.max, y_max = y.max, det_x = grid[,1], det_y = grid[,2], K = K, 
           M = M, J = J, n=n, d2int = d2int, n_i = n_i, tot_area = tot_area, 
           int_area = int_area, xobs = xobs)

aaa = "

data {
  int<lower=0> M;
  int<lower=0> n; // the number of observed individuals
  int<lower=0> K; // the number of traps
  int<lower=0> J; // the number of sampling occasions
  real<lower=0> x_max;
  real<lower=0> y_max;
  row_vector[K] det_x; // x-coordinates for traps
  row_vector[K] det_y;
  array[n] row_vector[K] xobs; // observed number of captures at each detector
  int<lower=0> n_i;
  matrix[n_i,K] d2int;
  real tot_area;
  real int_area;
}

parameters {
  real<lower=0, upper=1> p0; 
  real<lower=0, upper = 50> sigma;
  real<lower=n, upper=M> N;
  row_vector<lower=0,upper=x_max>[n] sx;
  row_vector<lower=0,upper=y_max>[n] sy;
}

model {
  array[n] row_vector[K] logp;
  array[n] row_vector[K] d2;
  matrix[n_i,K] logpint;
  vector[n_i] logpdot;
  real logint;
  vector[K] onevec_d = rep_vector(1.0, K);
  
  for(i in 1:n){
    d2[i] = (sx[i] - det_x).*(sx[i] - det_x) + (sy[i] - det_y).*(sy[i] - det_y);
    logp[i] = log(p0) - d2[i]/(2*sigma^2);
  }
  
  for(i in 1:n){
    target += xobs[i].*logp[i] + (J-xobs[i]).*log1m_exp(logp[i]); 
  }
  
  logpint = log(p0) - d2int / (2*sigma^2);    
  logpdot = J*log1m_exp(logpint) * onevec_d;
  logint = log_sum_exp(logpdot) - log(tot_area) + log(int_area);
  target += lgamma(N+1) - lgamma(N-n+1) + (N-n)*logint;
  
  // uniform priors as in NIMBLE
  
}



"

cat(aaa, file = "mod_scdl_bpp.stan")

stan_fit = cmdstan_model("mod_scdl_bpp.stan")

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

parout2 = as_draws_df(fit$draws(c("p0", "sigma", "N")))
mcmc_trace(parout2)
summary(parout2)
