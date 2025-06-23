# Author: Maddie Rainey, 2025 #

library(mgcv)

## df-Spatial+ model ##

# Input:
#   x = exposure from sim.data()
#   y = response from sim.data()
#   samp.grid = sampled data points
#   ind = indicator for additional measured covariates
#   cov = list of additional measured covariates
#   AIC = Indicator for AIC (TRUE) or BIC (FALSE).
#   beta = true value of exposure covariate
#   LOGISTIC_DAT = indicates whether or not outcome is binary
#   alpha = error rate for confidence interval
#   max.df = maximum df used in models where information criterion are computed
# Output:
#   beta = the estimated exposure coefficient of the df-spatial+ model
#   se = the standard error of the exposure coefficient of the df-spatial+ model
#   coverage = indicates if (1-alpha)% CI covers the the true beta 
#   k.used.e = number of df used in exposure model 
#   k.used.r = number of df used in outcome model
#   reject = indicates if p-value is less than alpha (1 = p-value < alpha, 0 p-value >= alpha)

knot.spat.plus.mod <- function(x, y, samp.grid, ind, cov, AIC = TRUE, beta, LOGISTIC_DAT = FALSE, alpha = 0.05, max.df = 400){
  
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
  #  print(df)
    mod <- lm(x ~ tprs[, 1:df])
    if(AIC){
      val[i] <- AIC(mod)
    }else{
      val[i] <- BIC(mod)
    }
  }
  
  df.min <- df.check[which(val == min(val))[[1]]]
  
  #obtain spatial residuals for x
  x.gam <- mgcv::gam(x~s(samp.grid[,1], samp.grid[,2], k = df.min+1, fx = T))
  r.x <- x.gam$residuals
  
  #fit the model and output summary statistics
  if(LOGISTIC_DAT){
    if(ind == 1){
      knot.spat.plus <- mgcv::gam(y ~ r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + s(samp.grid[,1], samp.grid[,2], k = df.min+1, fx = T), family = binomial)
    }else{
      knot.spat.plus <- mgcv::gam(y ~ r.x + s(samp.grid[,1], samp.grid[,2], k = df.min+1, fx = T), family = binomial)
    }
  }else{
    if(ind == 1){
      knot.spat.plus <- mgcv::gam(y ~ r.x + cov$z4 + as.factor(cov$z5) + cov$z6 + s(samp.grid[,1], samp.grid[,2], k = df.min+1, fx = T))
    }else{
      knot.spat.plus <- mgcv::gam(y ~ r.x + s(samp.grid[,1], samp.grid[,2], k = df.min+1, fx = T))
    }
  }
  #Get summary values to get standard errors
  knot.spat.plus.summary <- summary(knot.spat.plus)
  
  #Get Coverage
  lb <- as.numeric(knot.spat.plus$coefficients[2]) - qnorm(1 - alpha/2)*as.numeric(knot.spat.plus.summary$se[2])
  ub <- as.numeric(knot.spat.plus$coefficients[2]) + qnorm(1 - alpha/2)*as.numeric(knot.spat.plus.summary$se[2])
  
  covers <- ifelse(lb <= beta & beta <= ub, 1, 0)
  
  #get indicator if p-value rejects or fails to reject
  
  reject <- ifelse(as.numeric(knot.spat.plus.summary$p.pv[2]) < alpha, 1, 0)
  
  #return beta
  return(list("beta" = as.numeric(knot.spat.plus$coefficients[2]),
              "se" = as.numeric(knot.spat.plus.summary$se[2]),
              "coverage" = covers,
              "k.used.e" = df.min,
              "k.used.r" = df.min,
              "reject" = reject))
  
}
