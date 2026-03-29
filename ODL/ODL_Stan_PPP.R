library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)

### get the blackbear data ###
tmp = readRDS(here("data","data.rds")) # gets data summaries that we use in all models
tmp2 = readRDS(here("data","int_data.rds"))
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

###### inits
# psi
# p0
# sigma
# z
# sxy
inits = function(){
  return(list(p0 = runif(1, 0.25, 0.75),
              sigma = rlnorm(1),
              phi = runif(1,0.25,1)))
}
dat = list(x_max = x.max, y_max = y.max, det_x = grid[,1], det_y = grid[,2], K = K, 
           J = J, n=n, d2int = d2int, n_i = n_i, tot_area = tot_area, 
           int_area = int_area, xobs = t(xobs), phimax = 1)



aaa = "

data {
  int<lower=0> n; // the number of observed individuals
  int<lower=0> K; // the number of traps
  int<lower=0> J; // the number of sampling occasions
  real<lower=0> x_max;
  real<lower=0> y_max;
  row_vector[K] det_x; // x-coordinates for traps
  row_vector[K] det_y;
  matrix[K,n] xobs; // observed number of captures at each detector (it is transposed to previous representations)
  int<lower=0> n_i; // number of integration points
  matrix[n_i,K] d2int; // distances between integration points and each trap
  real tot_area;
  real int_area;
  real phimax;
}


parameters {
  real<lower=0, upper=1> p0; 
  real<lower=0, upper = 50> sigma;
  real<lower=0, upper = phimax> phi;
}

transformed parameters {
  real EN= phi*tot_area;
}

model {
  matrix[n_i, K] logp;
  matrix[n_i, K] log1mp;
  matrix[n_i,n] capt_obs;
  real some;
  vector[n_i] capt_mis;
  real Lambda;

  //  likelihood component for n observed individuals
  logp = log(p0) - d2int / (2*sigma*sigma); // n_i x K
  log1mp = log1m_exp(logp);
  capt_obs = int_area*phi*exp(logp * xobs + log1mp * (J - xobs));  // xobs is K x n matrix (transpose of typical)
  some = log(rep_row_vector(1.0, n_i) * capt_obs)*rep_vector(1.0, n); // takes the sum over integration points: capt is n x nint, ones is nint x 1
  target += some; // some gives the log likelihood contribution for each individual

  // likelihood component for e^-Lambda
  capt_mis = int_area*phi*(1-exp(log1mp * rep_vector(J,K))); //log1mp is n_i x K
  Lambda = sum(capt_mis);
  target += -Lambda;
 
  // priors: for now continue with uniform priors

}



"

cat(aaa, file = here("ODL","mod_odl_ppp.stan"))

stan_fit = cmdstan_model(here("ODL","mod_odl_ppp.stan"))

t1 = proc.time()
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
t2 = proc.time()
partime4 = t2-t1

parout4 = as_draws_df(fit$draws(c("p0", "sigma", "EN", "phi")))
mcmc_trace(parout4)
summary(parout4)
