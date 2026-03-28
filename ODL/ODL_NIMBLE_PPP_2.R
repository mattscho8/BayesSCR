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
           n=n, n_i = n_i, tot_area = tot_area, int_area = int_area, 
           phimax = 1)
dat = list(y = xobs, one = 1, zero = 0, zero2 = 0, d2int = d2int, one_row = matrix(1,1,n_i), one_col = matrix(1,n,1))


code <- nimbleCode({
  
  sigma ~ dunif(0, 50)
  alpha <- -1 / (2 * sigma^2)
  p0 ~ dunif(0, 1)
  
  logp <- log(p0) + alpha * d2int
  log1mp <- log(1-exp(logp))
  capt_obs <- int_area*phi*exp(logp %*% xobs + log1mp %*% (J - xobs))
  #some <- -log(one_row %*% capt_obs) %*% one_col
  #zero2 ~ dpois(some)
  
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