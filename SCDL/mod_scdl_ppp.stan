

data {
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
  real phimax;
}

parameters {
  real<lower=0, upper=1> p0; 
  real<lower=0, upper = 50> sigma;
  real<lower=0, upper = phimax> phi;
  row_vector<lower=0,upper=x_max>[n] sx;
  row_vector<lower=0,upper=y_max>[n] sy;
}

transformed parameters {
  real EN= phi*tot_area;
}

model {
  array[n] row_vector[K] logp;
  array[n] row_vector[K] d2;
  matrix[n_i,K] logpint;
  vector[n_i] logpdot;
  vector[K] onevec_d  = rep_vector(1.0, K);
  
  for(i in 1:n){
    d2[i] = (sx[i] - det_x).*(sx[i] - det_x) + (sy[i] - det_y).*(sy[i] - det_y);
    logp[i] = log(p0) - d2[i]/(2*sigma^2);
  }
  
  for(i in 1:n){
    target += xobs[i].*logp[i] + (J-xobs[i]).*log1m_exp(logp[i]); 
  }
  
  logpint = log(p0) - d2int / (2*sigma^2);    
  logpdot = J*log1m_exp(logpint) * onevec_d;
  target += -(n_i - sum(exp(logpdot)))*phi*int_area;
  
  target += n*log(phi);
}



