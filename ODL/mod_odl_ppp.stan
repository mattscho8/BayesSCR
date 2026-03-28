

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



