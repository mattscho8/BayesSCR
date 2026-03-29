library(nimble)
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

M = 150

##### constants
# x.max, y.max are same
# det.xy is grid
# n.d is ntrap
# M is M
# trials is k

##### data
# y is yfull

###### inits
# psi
# p0
# sigma
# z
# sxy
inits = list(p0 = runif(1),
             sigma = rlnorm(1),
             sxy = cbind(runif(M,0,x.max), runif(M,0,y.max),
             N = runif(1,n,M))
)
con = list(x.max = x.max, y.max = y.max, det.xy = grid, K = K, M = M, J = J, 
           n=n, d2int = d2int, n.i = n_i, tot_area = tot_area, int_area = int_area)
dat = list(y = xobs, one = 1, zero = 0)

code <- nimbleCode({
  for(i in 1:n) {
    sxy[i,1] ~ dunif(0, x.max)
    sxy[i,2] ~ dunif(0, y.max)
  }

  sigma ~ dunif(0, 50)
  alpha <- -1 / (2 * sigma^2)
  p0 ~ dunif(0, 1)
  for(i in 1:n) {
    d2[i, 1:K] <- (sxy[i,1] - det.xy[1:K,1])^2 + (sxy[i,2] - det.xy[1:K,2])^2
    p[i, 1:K] <- p0 * exp(alpha * d2[i,1:K])
    for(j in 1:K){
      y[i,j] ~ dbinom(prob = p[i,j], size = J)
    }
  }
  
  for(i in 1:n.i){
    pint[i,1:K] <- p0 * exp(alpha * d2int[i,1:K])
    ptmp[i,1:K] <- (1-pint[i,1:K])^J
    pdot[i] <- prod(ptmp[i,1:K])
  }
  intuse <- sum(pdot[1:n.i])*int_area/tot_area 
  one ~ dbinom(prob = intuse^(N-n), size = 1) # include the integral
  zero ~  dpois(loggam(M+1)-loggam(N+1)+loggam((N-n+1))) # include N!/(N-n)!
  
  N ~ dunif(n,M)
})

## Create and compile the NIMBLE model
model <- nimbleModel( code = code,
                      constants = con,
                      data = dat,
                      inits = inits,
                      check = F,
                      calculate = F)

# ## Check the initial log-likelihood 
# model$calculate()

cmodel <- compileNimble(model)

MCMCconf <- configureMCMC(model = model,
                          monitors  = c("N","sigma","p0"),
                          control = list(reflective = TRUE),
                          thin = 10)

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC)

t1 = proc.time()
samples <- runMCMC( mcmc = cMCMC,
                    nburnin = 1000,
                    niter = 5000,
                    nchains = 4,
                    samplesAsCodaMCMC = TRUE))
t2 = proc.time()
outtime2 = t2-t1

out_df2 = as_draws_df(samples)
mcmc_trace(out_df2)
summary(out_df2)
