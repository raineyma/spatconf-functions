# Author: Maddie Rainey, 2025 #

library(mvtnorm)
library(mgcv)

## Functions to Simulate data for methods ##

## Function 1 - Create and sample from grid, outputting sample grid ##
## Function 2 - Create distance matrix ##
## Function 3 - Calculate covariance matrices and simulate spatial field and build data ##

## Function 1 ##
# Input:
#   n = number of sampled grid points
#   x.size = length of one side of grid
#   y.size = length of one side of grid
#   step = step size increment for grid spacing
# Output:
#   samp.grid = nxn matrix of sampled points

sim.grid <- function(n, x.size, y.size, step, seed = seed){
  
  set.seed(seed)
  
  #create sampling grid
  grid <- expand.grid(seq(0, x.size, step), seq(0, y.size, step))
  
  #sample grid
  obs <- sample(c(1:nrow(grid)), size = n, replace = F)
  samp.grid <- grid[obs,]
  
  return(samp.grid)
  
}


## Function 2 ##
# Input:
#  samp.grid = sampled point matrix simulated from the sim.grid() function
# Output: 
#  h.mat = distance matrix, euclidean distance between two point

dist.mat <- function(samp.grid){
  as.matrix(dist(samp.grid))
}

## Function 3 ##
# Input:
#   p = probability input for logit function
# Output:
#   x = real line value output from logit calculation

logit <- function(p){
  x <- log(p/(1-p))
  return(x)
}

## Function 4 ##
# Input:
#   x = real line value input for expit (inverse logit) function
# Output:
#   p = probability output from expit calculation

expit <- function(x){
  p <- 1/(exp(-x) + 1)
  return(p)
}


## Function 5 ##
# Input:
#   n = number of sampled grid points
#   samp.grid = sampling grid simulated to make the distance matrix
#   h.mat = distance matrix computed from dist.mat()
#   R.1 = range parameter for first exponential covaraince matrix
#   R.2 = range parameter for second exponential covaraince matrix
#   R.3 = range parameter for third exponential covaraince matrix
#   ind.1 = Indicate whether to included non-spatial variation in exposure (0 = No, 1 = Yes)
#   ind.2 = Indicate whether second spatial structure is used (0 = No, 1 = Yes)
#   ind.3 = Indicate whether third spatial structure is used (0 = No, 1 = Yes)
#   ind.4 = Indicate whether spatial confounding is present (0 = No, 1 = Yes)
#   ind.5 = Indicate whether additional covariates are to be simulated and incorporated (0 = No, 1 = Yes)
#   d.1 = dampening parameter for x
#   d.2 = dampening parameter for f
#   sigma.x = standard deviation of the non-spatial errors for X
#   sigma.y = standard deviation of the non-spatial errors for Y
#   beta = true coefficient for exposure
#   beta.f = coefficient for confounder
#   PROJECT = Indicator if spatial fields should be fitted to tprs
#   k = number of knots for tprs used when PROJECT = TRUE
#   LOGISTIC_DAT = indicates whether or not simulated data is Logistic
#   p = baseline probability of experience an event in logistic regression
# Output:
#   y = outcome variable
#   x = exposure variable
#   f = confounder variable
#   z1 = confounded spatial field
#   z2 = spatial field unique to exposure
#   z3 = spatial field unique to confounder
#   z4 = spatial continuous measured confounder
#   z5 = non-spatial binary measured confounder
#   z6 = non-spatial continuous measured confounder
#   prob = probability used to sample a binary outcome


sim.data <- function(n = 1000, samp.grid = NULL, h.mat, 
                     R.1, R.2, R.3, ind.1, ind.2, ind.3, ind.4, ind.5, d.1, d.2,
                     sigma.x, sigma.y = NULL, sigma.p = NA, beta, beta.f, 
                     PROJECT=FALSE, k=301, LOGISTIC_DAT = FALSE, p = 0.05, seed = seed){
  
  set.seed(seed)
  
  #Create covariance matrices
  if(ind.4 == 1){
    Cov.mat.1 <- exp(-(h.mat/R.1))
  }
  if(ind.2 == 1){
    Cov.mat.2 <- exp(-(h.mat/R.2))
  }
  if(ind.3 == 1){
    Cov.mat.3 <- exp(-(h.mat/R.3))
  }
  if(ind.5 == 1){
    Cov.mat.4 <- exp(-(h.mat/25))
  }
  
  #Create spatial fields for data simulation
  
  if(ind.4 == 1){
    z1 <- t(mvtnorm::rmvnorm(1, mean = rep(0, n), sigma = Cov.mat.1))
   # z1 <- scale(z1, scale = F)
  }else{
    z1 <- NA
  }
  if(ind.2 == 1){
    z2 <- t(mvtnorm::rmvnorm(1, mean = rep(0, n), sigma = Cov.mat.2)) 
   # z2 <- scale(z2, scale = F)
  }else{
    z2 <- NA
  }
  if(ind.3 == 1){
    z3 <- t(mvtnorm::rmvnorm(1, mean = rep(0, n), sigma = Cov.mat.3)) 
   # z3 <- scale(z3, scale = F)
  }else{
    z3 <- NA
  }
  if(ind.5 == 1){
    z4 <- t(mvtnorm::rmvnorm(1, mean = rep(0, n), sigma = Cov.mat.4)) 
    z5 <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(.45,.45,.1))
    z6 <- rnorm(n = n, mean = 0, sd = 2)
    if(ind.4 == 1){
      cov <- 0.7*(z4+z1) + 0.5*as.numeric(I(z5 == 1)) - 0.45*as.numeric(I(z5 == 2)) + 1.25*z6
    }else{
      cov <- 0.7*z4 + 0.5*as.numeric(I(z5 == 1)) - 0.45*as.numeric(I(z5 == 2)) + 1.25*z6
    }
  }else{
    z4 <- NA
    z5 <- NA
    z6 <- NA
    cov <- rep(0,n)
  }
  
  
  
  #If PROJECT = TRUE, the project spatial fields onto tprs and use fitted values
  #when simulating data
  
  if(PROJECT){
    if(ind.4 == 1){
      z1 <- fitted(mgcv::gam(z1~s(samp.grid[,1], samp.grid[,2], k = k, fx = F)))
    }else{
      z1 <- NA
    }
    if(ind.2 == 1){
      z2 <- fitted(mgcv::gam(z2~s(samp.grid[,1], samp.grid[,2], k = k, fx = F))) 
    }else{
      z2 <- NA
    }
    if(ind.3 == 1){
      z3 <- fitted(mgcv::gam(z3~s(samp.grid[,1], samp.grid[,2], k = k, fx = F)))
    }else{
      z3 <- NA
    }
  }
  
  #Simulate exposure and outcome
  
  x.dat <- rep(0, n)
  if(ind.1 == 1){
    err.x <- rnorm(n, mean = 0, sd = sigma.x)
    x.dat <- x.dat + err.x
  }
  if(ind.4 == 1){
    x.dat <- x.dat + d.2*z1
  }
  if(ind.2 == 1){
    x.dat <- x.dat + z2
  }
  x.dat <- d.1*x.dat

  f <- rep(0, n)
  if(ind.4 == 1){
    f <- f + z1
  }
  if(ind.3 == 1){
    f <- f + z3 
  }
  f <- d.2*f
  
  
  if(LOGISTIC_DAT){
    beta.0 <- logit(p)
    prob <- expit(beta.0 + beta*x.dat + beta.f*f)
    y.dat <- rbinom(n, 1, expit(beta.0 + beta*x.dat + beta.f*f))
  }else{
    prob <- NA
    err.y <- rnorm(n, mean = 0, sd = sigma.y)
    y.dat <- beta*x.dat + beta.f*f + cov +  err.y
  }
  
  #Return data matrix
  return(list(y = y.dat, x = x.dat, f = f, z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, prob = prob))
}


