

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

