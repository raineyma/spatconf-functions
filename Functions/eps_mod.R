# Author: Maddie Rainey, 2024 #

library(mgcv)

## Exposure-Penalized Spline model ##

# Input:
#   x = exposure from sim.data()
#   y = response from sim.data()
#   samp.grid = sampled data points
#   k = df used in TPRS basis
#   beta = true value of exposure coefficient
#   k = number of knots for splines
#   LOGISITIC_DAT = indicates whether or not outcome is binary
#   alpha = error rate for confidence interval
# Output:
#   beta = the estimated exposure coefficient of the spatial+ model
#   se = the standard error of the exposure coefficient of the spatial+ model
#   coverage = 1 if 95% confidence interval covers true value of beta, 0 otherwise 
#   k.used = number of knots used in the tprs
#   coverage = indicates if (1-alpha)% CI covers the the true beta 
#   k.used.e = number of df used in exposure model 
#   k.used.r = number of df used in outcome model
#   reject = indicates if p-value is less than alpha (1 = p-value < alpha, 0 p-value >= alpha)

eps.mod <- function(x, y, samp.grid, k = 301, beta, LOGISTIC_DAT = FALSE, alpha = 0.05){
  
  ## Obtain smoothing penalty for second step
  x.gam <- mgcv::gam(x~s(samp.grid[,1], samp.grid[,2], k = k, fx = F))
  
  ## Fit E-PS model
  if(LOGISTIC_DAT){
    eps.mod <- mgcv::gam(y ~ x + s(samp.grid[,1], samp.grid[,2], k = k, fx = F, sp = x.gam$sp), family = binomial)
  }else{
    eps.mod <- mgcv::gam(y ~ x + s(samp.grid[,1], samp.grid[,2], k = k, fx = F, sp = x.gam$sp))
  }
  #Get summary values to get standard errors
  eps.summary <- summary(eps.mod)
  
  #Get Coverage
  lb <- as.numeric(eps.mod$coefficients[2]) - qnorm(1 - alpha/2)*as.numeric(eps.summary$se[2])
  ub <- as.numeric(eps.mod$coefficients[2]) + qnorm(1 - alpha/2)*as.numeric(eps.summary$se[2])
  
  covers <- ifelse(lb <= beta & beta <= ub, 1, 0)
  
  #gets approximate knots used

  k.used <- sum(eps.mod$edf) - 1

  #get indicator if p-value rejects or fails to reject
  
  reject <- ifelse(as.numeric(eps.summary$p.pv[2]) < alpha, 1, 0)
  
  return(list("beta" = as.numeric(eps.mod$coefficients[2]),
              "se" = as.numeric(eps.summary$se[2]),
              "coverage" = covers,
              "k.used.e" = k.used,
              "k.used.r" = k.used,
              "reject" = reject))
}
