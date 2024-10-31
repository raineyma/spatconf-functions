# Author: Maddie Rainey, 2024 #

library(mgcv)

## KS model ##

# Input:
#   x = exposure from sim.data()
#   y = outcome from sim.data()
#   samp.grid = sampled data points
#   AIC = Indicator for AIC or BIC. (TRUE = AIC, FALSE = BIC)
#   beta = true value of exposure coefficient
#   LOGSITIC_DAT = indicates whether or not outcome is binary
#   alpha = error rate for confidence interval
#   max.df = maximum degrees of freedom to check
#   fx = indicates whether to select number of df using information criterion
#   df = number of df to use in model if fx = TRUE
# Output:
#   beta = The estimated exposure coefficient from tprs adjusted model
#   se = The estimated standard error of the exposure from tprs adjusted model
#   coverage = indicates if (1-alpha)% CI covers the the true beta
#   k.used.e = number of df used in exposure model (WILL BE NA)
#   k.used.r = number of df used in outcome model
#   reject = indicates if p-value is less than alpha (1 = p-value < alpha, 0 p-value >= alpha)


tprs.adjust.mod <- function(x, y, samp.grid, AIC = FALSE, beta, LOGISTIC_DAT = FALSE, alpha = 0.05, max.df = 400, fx = F, df = 0){
  
  if(!fx){
  #dimension of tprs
  k.tprs <- .75*dim(samp.grid)[1]
  
  #transform data to list
  x.1 <- samp.grid[,1]
  x.2 <- samp.grid[,2]
  df.list <- list(x.1 = x.1, x.2 = x.2)
  
  #Get TPRS Basis
  tprs <- mgcv::smooth.construct(mgcv::s(x.1, x.2, k = k.tprs+1, fx = T), data = df.list, knots = NULL)$X
  tprs <- tprs[, c(k.tprs, k.tprs+1, 1:(k.tprs-2))]
  
  #Fit models different df of Basis
  if(max.df <= 50){
    df.check <- seq(3,max.df,1)
  } else if(max.df <= 100){
    df.check <- c(seq(3, 50, 1), seq(55, max.df, 5))
  }else{
    df.check <- c(seq(3, 50, 1), seq(55, 100, 5), seq(100, max.df, 25))
  }
  
  val <- numeric(length(df.check))
  
  for(i in 1: length(df.check)){
    df <- df.check[i]
   # print(df)
    if(LOGISTIC_DAT){
      mod <- glm(y ~ tprs[, 1:df], family = binomial)
    }else{
      mod <- lm(y ~ tprs[, 1:df])
    }
    if(AIC){
      val[i] <- AIC(mod)
    }else{
      val[i] <- BIC(mod)
    }
  }
  
  df <- df.check[which(val == min(val))[[1]]]
  }
  
  #Fit final model
  if(LOGISTIC_DAT){
    tprs.adjust <- mgcv::gam(y ~ x + s(samp.grid[,1], samp.grid[,2], k = df+1, fx = T), family = binomial)
  }else{
    tprs.adjust <- mgcv::gam(y ~ x + s(samp.grid[,1], samp.grid[,2], k = df+1, fx = T))
  }
  
  #Get summary values to get standard errors
  tprs.summary <- summary(tprs.adjust)
  
  #Compute coverage
  lb <- as.numeric(tprs.adjust$coefficients[2]) - qnorm(1 - alpha/2)*as.numeric(tprs.summary$se[2])
  ub <- as.numeric(tprs.adjust$coefficients[2]) + qnorm(1 - alpha/2)*as.numeric(tprs.summary$se[2])
  
  covers <- ifelse(lb <= beta && beta <= ub, 1, 0)
  
  #get indicator if p-value rejects or fails to reject
  
  reject <- ifelse(as.numeric(tprs.summary$p.pv[2]) < alpha, 1, 0)
  

  #return betas, standard errors and coverage
  return(list("beta" = as.numeric(tprs.adjust$coefficients[2]),
              "se" = as.numeric(tprs.summary$se[2]),
              "coverage" = covers,
              "k.used.e" = NA,
              "k.used.r" = df,
              "reject" = reject))
  
}


