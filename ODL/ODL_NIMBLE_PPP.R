library(nimble)
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

###### inits
# psi
# p0
# sigma
# z
# sxy
inits = list(p0 = runif(1),
             sigma = rlnorm(1))

con = list(K = K, J = J, 
           n=n, d2int = d2int, n_i = n_i, tot_area = tot_area, int_area = int_area, 
           phimax = 1)
dat = list(y = xobs, one = 1, zero = 0, vone = rep(1,n))


code <- nimbleCode({
  
  sigma ~ dunif(0, 50)
  alpha <- -1 / (2 * sigma^2)
  p0 ~ dunif(0, 1)
  p[1:n_i,1:K] <- p0 * exp(alpha * d2int[1:n_i,1:K])
  for(i in 1:n) {
    for(j in 1:K){
      tmp1[i,1:n_i,j] <- p[1:n_i,j]^y[i,j] * (1-p[1:n_i,j])^(J-y[i,j])
    }
    for(h in 1:n_i){
      tmp2[i,h] <- prod(tmp1[i,h,1:K])
    }
    like[i] <- int_area*phi*sum(tmp2[i,1:n_i])
    vone[i] ~ dbern(like[i])
  }
  
  for(i in 1:n_i){
    pint[i,1:K] <- p0 * exp(alpha * d2int[i,1:K])
    ptmp[i,1:K] <- (1-pint[i,1:K])^J
    ptot[i] <- prod(ptmp[i,1:K])
    pdot[i] <- 1-ptot[i]
  }
  lam <- sum(pdot[1:n_i])*phi*int_area # area 1 gives total area, area 2 is area of each cell
  zero ~ dpois(lam)
  
  one ~ dbinom(p = (phi/phimax)^n, size = 1) # gets lambda() into model
  
  phi ~ dunif(0,phimax)
  
  EN <- phi*tot_area
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
                          monitors  = c("EN","sigma","phi","p0"),
                          control = list(reflective = TRUE))

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC)

MCMCRuntime <- system.time(samples <- runMCMC( mcmc = cMCMC,
                                               nburnin = 1000,
                                               niter = 5000,
                                               nchains = 4,
                                               samplesAsCodaMCMC = TRUE))

out_df3 = as_draws_df(samples)
mcmc_trace(out_df3)
summary(out_df3)