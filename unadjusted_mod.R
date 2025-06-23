# Author: Maddie Rainey, 2025 #


## Spatally-Unadjusted model ##

# Input:
#   x = exposure from sim.data()
#   y = outcome from sim.data()
#   beta = true value of exposure coefficient
#   ind = indicator for additional measured covariates
#   cov = additional measured covariates
#   LOGISITC_DAT = indicates whether or not outcome is binary
#   alpha = error rate for confidence interval
# Output:
#   beta = the estimated exposure coefficient of the lm model
#   se = the standard error of the exposure coefficient of the lm model
#   coverage = indicates if (1-alpha)% CI covers the the true beta
#   k.used = NA, no knot selection is used for this regression
#   reject = indicates if p-value is less than alpha (1 = p-value < alpha, 0 p-value >= alpha)

null.mod <- function(x, y, beta, ind, cov, LOGISTIC_DAT = FALSE, alpha = 0.05){
  
  #Fit base model
  if(LOGISTIC_DAT){
    if(ind == 1){
      mod <- glm(y ~ x + cov$z4 + as.factor(cov$z5) + cov$z6, family = binomial)
    }else{
      mod <- glm(y ~ x, family = binomial)
    }
  }else{
    if(ind == 1){
      mod <- lm(y ~ x + cov$z4 + as.factor(cov$z5) + cov$z6)
    }else{
      mod <- lm(y ~ x)
    }
  }
  #Get summary information
  summary.mod <- summary(mod)
  
  #get coverage
  lb <- as.numeric(coef(mod)[2]) - qnorm(1 - alpha/2)*as.numeric(summary.mod$coefficients[2,2])
  ub <- as.numeric(coef(mod)[2]) + qnorm(1 - alpha/2)*as.numeric(summary.mod$coefficients[2,2])
  
  covers <- ifelse(lb <= beta && beta <= ub, 1, 0)
  
  #get indicator if p-value rejects or fails to reject
  
  reject <- ifelse(as.numeric(summary.mod$coefficients[2,4]) < alpha, 1, 0)
  
  #Return beta
  return(list("beta" = as.numeric(coef(mod)[2]),
              "se" = as.numeric(summary.mod$coefficients[2,2]),
              "coverage" = covers,
              "k.used.e" = NA,
              "k.used.r" = NA,
              "reject" = reject))
  
}
