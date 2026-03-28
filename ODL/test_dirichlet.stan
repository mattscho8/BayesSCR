

data {
  integer<lower = 0> k; // number of categories
  integer<lower = 0> n_s; // number of studies
  array[n_s, k] integer y1; // could make this a 3D array
  array[n_s, k] integer y2; 
}

parameters{
  simplex[k] pie_overall;
  array[2] simplex[k] pie_region;
  array[n_s] simplex[k] pie_study_1;
  array[n_s] simplex[k] pie_study_2;
  vector<lower=0>[2] gamma;
}

model{
  for(i in 1:n_s){
    y1[i] ~ multinomial(pie_study_1[i]);
    y2[i] ~ multinomial(pie_study_2[i]);
    
    pie_study_1[i] ~ dirichlet(gamma[1]*pie_region[1]);
    pie_study_2[i] ~ dirichlet(gamma[1]*pie_region[2]);
  }
  
  pie_region[1] ~ dirichlet(gamma[2]*pie_overall);
  pie_region[2] ~ dirichlet(gamma[2]*pie_overall);
  
  pie_overall ~ dirichlet(rep_vector(0.5, k));
  
}

