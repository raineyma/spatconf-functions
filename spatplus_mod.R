# Author: Maddie Rainey, 2025 #

library(mgcv)

## Spatial+ model ##

# Input:
#   x = exposure from sim.data()
#   y = response from sim.data()
#   samp.grid = sampled data points
#   ind = indicator for additional measured covariates
#   cov = list of additional measured covariates
#   k = df used for TPRS basis
#   fx = Indicator for spatial smoothing (0 = smoothing penalty, 1 = no penalty)
#   beta = true value of exposure coefficient
#   LOGISTIC_DAT = indicates whether or not outcome is binary
#   alpha = error rate for confidence interval
# Output:
#   beta = the estimated exposure coefficient of the spatial+ model
#   se = the standard error of the exposure coefficient of the spatial+ model
#   coverage = indicates if (1-alpha)% CI covers the the true beta 
#   k.used.e = number of df used in exposure model 
#   k.used.r = number of df used in outcome model
#   reject = indicates if p-value is less than alpha (1 = p-value < alpha, 0 p-value >= alpha)

spat.plus.mod <- function(x, y, samp.grid, ind, cov, k = 301, fx, beta, LOGISTIC_DAT = FALSE, alpha = 0.05){
  
  #obtain spatial residuals for x
  ifelse(fx == 0,
         x.gam <- mgcv::gam(x~s(samp.grid[,1], samp.grid[,2], k = k, fx = F)),
         x.gam <- mgcv::gam(x~s(samp.grid[,1], samp.grid[,2], k = k, fx = T)))
  r.x <- x.gam$residuals
  
  k.used.exp <- sum(x.gam$edf) - 1
  
  #fit the Spatial+ model and output summary statistics
  if(LOGISTIC_DAT){
    if(ind == 1){
      ifelse(fx == 0,
             spat.plus <- mgcv::gam(y ~ r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + s(samp.grid[,1], samp.grid[,2], k = k, fx = F), family = binomial),
             spat.plus <- mgcv::gam(y ~ r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + s(samp.grid[,1], samp.grid[,2], k = k, fx = T), family = binomial))
    }else{
      ifelse(fx == 0,
             spat.plus <- mgcv::gam(y ~ r.x + s(samp.grid[,1], samp.grid[,2], k = k, fx = F), family = binomial),
             spat.plus <- mgcv::gam(y ~ r.x + s(samp.grid[,1], samp.grid[,2], k = k, fx = T), family = binomial))
    }
  }else{
    if(ind == 1){
      ifelse(fx == 0,
             spat.plus <- mgcv::gam(y ~ r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + s(samp.grid[,1], samp.grid[,2], k = k, fx = F)),
             spat.plus <- mgcv::gam(y ~ r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + s(samp.grid[,1], samp.grid[,2], k = k, fx = T)))
    }else{
      ifelse(fx == 0,
             spat.plus <- mgcv::gam(y ~ r.x + s(samp.grid[,1], samp.grid[,2], k = k, fx = F)),
             spat.plus <- mgcv::gam(y ~ r.x + s(samp.grid[,1], samp.grid[,2], k = k, fx = T)))
    }
  }
  #Get summary values to get standard errors
  spat.plus.summary <- summary(spat.plus)
  
  #Get Coverage
  lb <- as.numeric(spat.plus$coefficients[2]) - qnorm(1 - alpha/2)*as.numeric(spat.plus.summary$se[2])
  ub <- as.numeric(spat.plus$coefficients[2]) + qnorm(1 - alpha/2)*as.numeric(spat.plus.summary$se[2])
  
  covers <- ifelse(lb <= beta & beta <= ub, 1, 0)
  
  #gets approximate knots used
  
  if(fx == 0){
    k.used <- sum(spat.plus$edf) - 1
  }else{
    k.used <- k-1
  }
  #get indicator if p-value rejects or fails to reject
  
  reject <- ifelse(as.numeric(spat.plus.summary$p.pv[2]) < alpha, 1, 0)
  
  #return beta
  return(list("beta" = as.numeric(spat.plus$coefficients[2]),
              "se" = as.numeric(spat.plus.summary$se[2]),
              "coverage" = covers,
              "k.used.e" = k.used.exp,
              "k.used.r" = k.used,
              "reject" = reject))
  
}
