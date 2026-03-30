library(nimble)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(here)

### get the blackbear data ###
tmp = readRDS(here("data","data.rds"))
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
inits = list(psi = runif(1),
             p0 = runif(1),
             sigma = rlnorm(1),
             z = c(rep(1, n), rbinom(M-n,1,0.5)),
             sxy = cbind(runif(M,0,x.max), runif(M,0,y.max))
)
con = list(x.max = x.max, y.max = y.max, det.xy = grid, K = K, M = M, J = J)
dat = list(y = yfull)

code <- nimbleCode({
  for(i in 1:M) {
    sxy[i,1] ~ dunif(0, x.max)
    sxy[i,2] ~ dunif(0, y.max)
  }
  psi ~ dunif(0,1)
  for (i in 1:M) {
    z[i] ~ dbern(psi)
  }
  
  sigma ~ dunif(0, 50)
  alpha <- -1 / (2 * sigma^2)
  p0 ~ dunif(0, 1)
  for(i in 1:M) {
    d2[i, 1:K] <- (sxy[i,1] - det.xy[1:K,1])^2 + (sxy[i,2] - det.xy[1:K,2])^2
    p[i, 1:K] <- p0 * exp(alpha * d2[i,1:K])
    for(j in 1:K){
      y[i,j] ~ dbinom(prob = z[i]*p[i,j], size = J)
    }
  }
  N <- sum(z[1:M])
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
                          monitors  = c("N","sigma","psi","p0"),
                          control = list(reflective = TRUE))

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC)

t1 = proc.time()
samples <- runMCMC( mcmc = cMCMC,
                    nburnin = 1000,
                    niter = 6000,
                    nchains = 4,
                    samplesAsCodaMCMC = TRUE)
t2 = proc.time()
outtime = t2-t1

out_df = as_draws_df(samples)
mcmc_trace(out_df)
summary(out_df)
