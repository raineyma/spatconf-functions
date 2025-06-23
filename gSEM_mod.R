# Author: Maddie Rainey, 2025 #

library(mgcv)

## gSEM model ##

# Input:
#   x = exposure from sim.data()
#   y = outcome from sim.data()
#   samp.grid = sampled data points
#   ind = indicator for additional measured covariates
#   cov = list of additional measured covariates
#   k = number of knots for TPRS basis
#   fx = Indicator for spatial smoothing (0 = smoothing penalty, 1 = no penalty)
#   beta = true value of exposure coefficient
#   LOGISITC_DAT = indicates whether or not outcome is binary
#   alpha = error rate for confidence interval
# Output:
#   beta = the estimated exposure coefficient of the gSEM model
#   se = the standard error of the exposure coefficient of the gSEM model
#   coverage = indicates if (1-alpha)% CI covers the the true beta
#   k.used.e = number of df used in exposure model 
#   k.used.r = number of df used in outcome model
#   reject = indicates if p-value is less than alpha (1 = p-value < alpha, 0 p-value >= alpha)


gSEM.mod <- function(x, y, samp.grid, ind, cov, k = 301, fx, beta, LOGISTIC_DAT = FALSE, alpha = 0.05){
  
  #obtain spatial residuals for x
  ifelse(fx == 0,
         x.gam <- mgcv::gam(x~s(samp.grid[,1], samp.grid[,2], k = k, fx = F)),
         x.gam <- mgcv::gam(x~s(samp.grid[,1], samp.grid[,2], k = 301, fx = T)))
  r.x <- x.gam$residuals
  
  k.used.exp <- sum(x.gam$edf) - 1
  
  #obtain spatial residuals for y
  if(LOGISTIC_DAT){
    ifelse(fx == 0,
           y.gam <- mgcv::gam(y~s(samp.grid[,1], samp.grid[,2], k = k, fx = F), family = binomial),
           y.gam <- mgcv::gam(y~s(samp.grid[,1], samp.grid[,2], k = k, fx = T), family = binomial))
    r.y <- y.gam$residuals
   
  }else{
    ifelse(fx == 0,
           y.gam <- mgcv::gam(y~s(samp.grid[,1], samp.grid[,2], k = k, fx = F)),
           y.gam <- mgcv::gam(y~s(samp.grid[,1], samp.grid[,2], k = k, fx = T)))
    r.y <- y.gam$residuals
  }
  #fit the gSEM model and output summary statistics
  if(ind == 1){
    gSEM.lm <- lm(r.y~r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + 0)
  }else{
    gSEM.lm <- lm(r.y~r.x + 0)
  }
  
  #Get summary information
  summary.gSEM <- summary(gSEM.lm)
  
  #get coverage
  lb <- as.numeric(gSEM.lm$coefficients[1]) - qnorm(1 - alpha/2)*as.numeric(summary.gSEM$coefficients[1,2])
  ub <- as.numeric(gSEM.lm$coefficients[1]) + qnorm(1 - alpha/2)*as.numeric(summary.gSEM$coefficients[1,2])
  
  
  covers <- ifelse(lb <= beta & beta <= ub, 1, 0)
  
  #get indicator if p-value rejects or fails to reject
  
  reject <- ifelse(as.numeric(summary.gSEM$coefficients[1,4]) < alpha, 1, 0)
  
  #gets approximate knots used
  
  if(fx == 0){
    k.used <- sum(y.gam$edf) - 1
  }else{
    k.used <- k-1
  }
  
  #Return beta
  return(list("beta" = as.numeric(gSEM.lm$coefficients[1]),
              "se" = as.numeric(summary.gSEM$coefficients[1,2]),
              "coverage" = covers,
              "k.used.e" = k.used.exp,
              "k.used.r" = k.used,
              "reject" = reject))
}
