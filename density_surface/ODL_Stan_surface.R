library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(here)
library(fmesher)

### get the blackbear data ###
tmp = readRDS(here("data","data.rds")) # gets data summaries that we use in all models
tmp2 = readRDS(here("data","int_data.rds"))
list2env(tmp, .GlobalEnv)
list2env(tmp2, .GlobalEnv)


#### Plot
plot(grid, xlim = xmask, ylim = ymask, pch = 4, lwd = 2)
points(int_grid, pch = 16, col = "blue", cex = 0.5)


#### Set up the tensor spline
m = 4 ## cubic spline -- standard notation (note: Wood uses m = 2)
bk = 4 ## Seber's K (# of interior knots; different from Wood's k).  The number of parameters in the B-spline is K + m

bx = bk + 2*m
tmpdiffx = diff(xmask)/(bk+1)
tmpdiffy = diff(ymask)/(bk+1)
bxgrid = seq(from = xmask[1] - (m-1)*tmpdiffx, to = xmask[2] + (m-1)*tmpdiffx, by = tmpdiffx)
bygrid = seq(from = ymask[1] - (m-1)*tmpdiffy, to = ymask[2] + (m-1)*tmpdiffy, by = tmpdiffy)

bbasis = function(x,knot,m=4){
  # x is a vector of values to evaluate basis
  # knot is a vector of knot points
  nn = length(x)
  nk = length(knot) # total number of knots
  tmpbasis = array(0,c(nn,nk+2*m-1,m))
  for(j in 1:(nk-1)){
    tmpbasis[,j,1] = (x >= knot[j]) & (x < knot[j+1])
  }
  for(h in 2:m){
    for(j in 1:(nk-h)){
      tmpbasis[,j,h] = (x - knot[j])/(knot[j+h-1] - knot[j])*tmpbasis[,j,h-1] + (knot[j+h]-x)/(knot[j+h]-knot[j+1])*tmpbasis[,j+1,h-1]
    }
  }
  tmp = tmpbasis[,1:(nk-m),m]
  return(tmp)
}
xbasis = bbasis(int_grid[,1], bxgrid)
ybasis = bbasis(int_grid[,2], bygrid)
npar_spl = ncol(xbasis)
P = diff(diag(npar_spl))
Sx = Sy = t(P)%*%P
Sx_marg = kronecker(Sx, diag(npar_spl))
Sy_marg = kronecker(diag(npar_spl), Sy)
# 
# ### Get the basis
tot_basis_full = as.matrix(fm_row_kron(xbasis,ybasis))
# 
### Get sum to zero constraint
# Impose a sum-to-zero constraint on \sum(X\beta) (i.e. 1 %*% X %*% beta). Need QR on 1 *%* X
tmp = qr(t(matrix(1,nrow = 1, ncol = n_i) %*% tot_basis_full))
tmp2 = qr.Q(tmp, complete = TRUE)

# # ## sum to zero across the beta parameters
# tmp = qr(rep(1,npar_spl^2))
# tmp2 = qr.Q(tmp, complete = TRUE)

BZmat = tmp2[,2:(npar_spl^2)]
tot_basis = tot_basis_full %*% BZmat

Sx_use = t(BZmat) %*% Sx_marg %*% BZmat
Sy_use = t(BZmat) %*% Sy_marg %*% BZmat


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
  return(list(p0 = runif(1),
              sigma = rlnorm(1),
              sx = runif(n,0,x.max),
              sy = runif(n,0,y.max),
              beta0 = rnorm(1),
              lambda = rlnorm(2),
              delta = rnorm(npar_spl^2 - 1)))
}
dat = list(x_max = x.max, y_max = y.max, det_x = grid[,1], det_y = grid[,2], K = K, 
           J = J, n=n, d2int = d2int, n_i = n_i, tot_area = tot_area, 
           int_area = int_area, xobs = t(xobs), sk = npar_spl^2 - 1, 
           basis_mat = tot_basis, prior_x = Sx_use, prior_y = Sy_use)

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
  int<lower=0> sk; // the number of basis functions
  matrix[n_i,sk] basis_mat;
  matrix[sk,sk] prior_x;
  matrix[sk,sk] prior_y;
}


parameters {
  real<lower=0, upper=1> p0; 
  real<lower=0, upper = 50> sigma;
  real beta0;
  vector[sk] delta;
  vector<lower=0>[2] lambda;
}

model {
  matrix[n_i, K] logp;
  matrix[n_i, K] log1mp;
  matrix[n_i,n] capt_obs;
  vector[n_i] dens; 
  real some;
  vector[n_i] capt_mis;
  real Lambda;

  //  likelihood component for n observed individuals
  dens = exp(beta0 + basis_mat * delta);
  logp = log(p0) - d2int / (2*sigma*sigma); // n_i x K
  log1mp = log1m_exp(logp);
  capt_obs = int_area*exp(logp * xobs + log1mp * (J - xobs));  // xobs is K x n matrix (transpose of typical)
  some = log(dens' * capt_obs)*rep_vector(1.0, n); // takes the sum over integration points: capt is n x nint, ones is nint x 1
  target += some; // some gives the log likelihood contribution for each individual

  // likelihood component for e^-Lambda
  capt_mis = int_area*dens.*(1-exp(log1mp * rep_vector(J,K))); //log1mp is n_i x K
  Lambda = sum(capt_mis);
  target += -Lambda;
 
  // priors: for now continue with uniform priors
  delta ~ multi_normal_prec(rep_vector(0,sk), 1/lambda[1]/lambda[1]*prior_x + 1/lambda[2]/lambda[2]*prior_y);
  lambda ~ student_t(3,0,5);
}



"

cat(aaa, file = here("density_surface","mod_dens.stan"))

stan_fit = cmdstan_model(here("density_surface","mod_dens.stan"))

t1 = proc.time()
fit <- stan_fit$sample(
  data = dat,
  init = inits,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 5000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 50
)
t2 = proc.time()
partimed1 = t2-t1

parout = as_draws_df(fit$draws(c("p0", "sigma", "beta0", "lambda")))
mcmc_trace(parout)
summary(parout)


parout_del = as_draws_df(fit$draws(c("delta")))
deluse = select(parout_del,contains('delta'))
delmat = as.matrix(deluse)
storemat = exp(rep(1,n_i)%*%t(parout$beta0) + tot_basis%*%t(delmat))
meddens = apply(storemat,1,median)

plotdf = data.frame(x = (int_grid[,1]-w)*1000 + min(rawgrid[,1]), y = (int_grid[,2]-w)*1000 + min(rawgrid[,2]), z = meddens)
#plotdf = data.frame(x = int_grid[,1]*1000, y = int_grid[,2]*1000, z = meddens)
plotdf2 = data.frame(x = rawgrid[,1], y = rawgrid[,2], nuniq = colSums(xobs>0), ntot = colSums(xobs))
plotdf2$nuniq[plotdf2$nuniq==0] = NA
plotdf2$ntot[plotdf2$ntot==0] = NA

# ### # unique
# ggplot(plotdf, aes(x = x, y = y, z = z)) +
#   geom_raster(data = plotdf, mapping = aes(fill = z)) +
#   # geom_contour(color= 'black', bins = 8, alpha = 0.25) +
#   scale_fill_viridis_c() +
#   xlab('UTM Zone 18 Easting') +
#   ylab('UTM Zone 18 Northing') +
#   coord_cartesian(expand = FALSE) +
#   geom_point(data = plotdf2, mapping = aes(x = x, y = y, z = NULL, size=nuniq), shape = 21, fill = "coral", col = "coral") +
#   labs(size='# unique ind') +
#   geom_point(data = plotdf2, mapping = aes(x = x, y = y, z = NULL), shape = 21, size = 5, col = "white", stroke = 1) +
#   scale_size(range = c(1,5)) +
#   labs(fill='Density') +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank())
# #ggsave('../output/densityplot_uniq_bspl.pdf',width = 10, height = 8)

###### General surface plot ###########
ggplot(plotdf, aes(x = x, y = y, z = z)) +
  geom_raster(data = plotdf, mapping = aes(fill = z)) +
  # geom_contour(color= 'black', bins = 8, alpha = 0.25) +
  scale_fill_viridis_c() +
  xlab('UTM Zone 18 Easting') +
  ylab('UTM Zone 18 Northing') +
  coord_cartesian(expand = FALSE) +
  geom_point(data = plotdf2, mapping = aes(x = x, y = y, z = NULL), shape = 21, size = 3, col = "white", stroke = 1) +
  scale_size(range = c(1,5)) +
  labs(fill='Density') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
ggsave(here("output","density_surface.pdf"),width = 10, height = 8)

####### Plotting each individual separately
# Convert xobs to long format
# Long format: one row per individual x detector
xlong <- data.frame(
  ind  = as.vector(row(xobs)),
  det  = as.vector(col(xobs)),
  ndet = as.vector(xobs)
) %>%
  filter(ndet > 0) %>%                      # keep only detectors where the individual was seen
  left_join(plotdf2 %>% mutate(det = row_number()), by = "det")

ggplot(plotdf, aes(x = x, y = y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_viridis_c() +
  # detections for this individual only
  geom_point(data = xlong, aes(x = x, y = y, size = ndet),
             shape = 21, fill = "coral", col = "coral") +
  scale_size(range = c(1, 5)) +
  facet_wrap(~ ind, ncol = 8) +
  coord_cartesian(expand = FALSE) +
  labs(x = "UTM Zone 18 Easting", y = "UTM Zone 18 Northing",
       fill = "Density", size = "# detections") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(here("output","density_ind.pdf"),width = 10, height = 8)
