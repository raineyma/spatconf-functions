# Author: Maddie Rainey, 2024 #

source("data_generating_functions.R")
source("unadjusted_mod.R")
source("gSEM_mod.R")
source("spatplus_mod.R")
source("ks_mod.R")
source("dfspatplus_mod.R")
source("eps_mod.R")


## Variables needed to run simulation ##
# n = number of sampled grid points
# num.sim = number of iteration of the simulation
# max.grid.x = size of one side of grid
# max.grid.y = size of one side of grid
# step.size = step size increment for grid spacing
# R.1 = range parameter for first exponential covaraince matrix
# R.2 = range parameter for second exponential covaraince matrix
# R.3 = range parameter for third exponential covaraince matrix
# ind.1 = Indicate whether sigma.x is used (0 = No, 1 = Yes)
# ind.2 = Indicate whether second spatial structure is used (0 = No, 1 = Yes)
# ind.3 = Indicate whether third spatial structure is used (0 = No, 1 = Yes)
# ind.4 = Indicate whether spatial confounding is present (0 = No, 1 = Yes)
# sigma.x = standard deviation used to simulate non-spatial error in x when ind.1 = 1
# sigma.y = standard deviation used to simulate non-spatial error in y
# beta = True coefficient of exposure
# beta.f = True coefficient of confounder
# PROJECT = Indicator if spatial fields should be fit to TPRS basis
# LOGISTIC_DAT = indicates whether or not simulated data is Logistic
# p = expected value of response when covariate is at the mean (Logistic model)
# k.sim = number of knots for tprs used when PROJECT = TRUE for simulating data
# k = number of knots to use in the methods
# max.df = maximum number of splines to check for AIC/BIC selection methods
# gSEM.flag = Indicates whether gSEM models should be run
# spatplus.flag = Indicates whether Spatial+ models should be run
# tprs.flag = Indicates whether the KS models should be run
# tprs.fx.flag = Indicates whether the KS model with fixed 10 df should be run
# knotspat.flag = Indicates whether the df-Spatial+ models should be run
# eps.flag = Indicates whether the E-PS model should be run
# seed = set to be set for replicability purposes

n = 1000
num.sim = 10
max.grid.x = max.grid.y = 10
step.size = 0.2
ind.1 = 1
ind.2 = 1
ind.3 = 1
ind.4 = 1
R.1 = 150
R.2 = 10
R.3 = 100
sigma.x = 0.1
sigma.y = 1
k = 301
k.sim = 301
max.df = 300
beta = 1
beta.f = 1
gSEM.flag = TRUE
spatplus.flag = TRUE
tprs.flag = TRUE
tprs.fx.flag = TRUE
knotspat.flag = TRUE
eps.flag = TRUE
PROJECT = FALSE
LOGISTIC_DAT = FALSE
p = 0.15
seed = 11041997

## Set matrices to save results from simulation ##
lm.results <- matrix(NA, nrow = num.sim, ncol = 6)
gSEM.results <- matrix(NA, nrow = num.sim, ncol = 6)
gSEM.fx.results <- matrix(NA, nrow = num.sim, ncol = 6)
gSEM.fx10.results <- matrix(NA, nrow = num.sim, ncol = 6)
spatplus.results <- matrix(NA, nrow = num.sim, ncol = 6)
spatplus.fx.results <- matrix(NA, nrow = num.sim, ncol = 6)
spatplus.fx10.results <- matrix(NA, nrow = num.sim, ncol = 6)
tprs.AIC.results <- matrix(NA, nrow = num.sim, ncol = 6)
tprs.BIC.results <- matrix(NA, nrow = num.sim, ncol = 6)
tprs.df3.results <- matrix(NA, nrow = num.sim, ncol = 6)
tprs.df10.results <- matrix(NA, nrow = num.sim, ncol = 6)
knotspat.AIC.results <- matrix(NA, nrow = num.sim, ncol = 6)
knotspat.BIC.results <- matrix(NA, nrow = num.sim, ncol = 6)
eps.results <- matrix(NA, nrow = num.sim, ncol = 6)

col.names <- c("Beta", "SE", "Coverage", "knots exposure", "knots response", "reject")

colnames(lm.results) <- col.names
colnames(gSEM.results) <- col.names
colnames(gSEM.fx.results) <- col.names
colnames(gSEM.fx10.results) <- col.names
colnames(spatplus.results) <- col.names
colnames(spatplus.fx.results) <- col.names
colnames(spatplus.fx10.results) <- col.names
colnames(tprs.AIC.results) <- col.names
colnames(tprs.BIC.results) <- col.names
colnames(tprs.df3.results) <- col.names
colnames(tprs.df10.results) <- col.names
colnames(knotspat.AIC.results) <- col.names
colnames(knotspat.BIC.results) <- col.names
colnames(eps.results) <- col.names


## Run Simulation ##

for(i in 1:num.sim){
  
  # Simulate Data #
  test.grid <- sim.grid(n = n, x.size = max.grid.x, y.size = max.grid.y, step = step.size, seed = seed)
  test.dist <- dist.mat(samp.grid = test.grid)
  test.data <- sim.data(n = n, h.mat = test.dist, samp.grid = test.grid, 
                        R.1 = R.1, R.2 = R.2, R.3 = R.3, 
                        ind.1 = ind.1, ind.2 = ind.2, ind.3 = ind.3, ind.4 = ind.4, 
                        d.1 = 1, d.2 = 1,
                        sigma.x = sigma.x, sigma.y = sigma.y, beta = beta, beta.f = beta.f,
                        PROJECT = PROJECT, k=k.sim, LOGISTIC_DAT = LOGISTIC_DAT, p = p, seed = seed)
  
  #Run Each Model
  
  #Run spatially-unadjusted model
  lm <- null.mod(x = test.data$x, y = test.data$y, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)
  
  #Run gSEM models
  if(gSEM.flag){
    gSEM <- gSEM.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k = k, fx = 0, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)
    gSEM.fx <- gSEM.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k= k, fx = 1, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)
    gSEM.fx10 <- gSEM.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k= 11, fx = 1, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)
  }else{
    gSEM <- list("beta" = NA,
                 "se" = NA,
                 "coverage" = NA,
                 "k.used.e" = NA,
                 "k.used.r" = NA,
                 "reject" = NA)
    gSEM.fx <- list("beta" = NA,
                    "se" = NA,
                    "coverage" = NA,
                    "k.used.e" = NA,
                    "k.used.r" = NA,
                    "reject" = NA)
    gSEM.fx10 <- list("beta" = NA,
                      "se" = NA,
                      "coverage" = NA,
                      "k.used.e" = NA,
                      "k.used.r" = NA,
                      "reject" = NA)
  }
  
  #Run Spatial+ models
  if(spatplus.flag){
    print('Spatial+ GCV')
    spatplus <- spat.plus.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k = k, 
                              fx = 0, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)
    print('Spatial+ Fixed')
    spatplus.fx <- spat.plus.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k = k, 
                                 fx = 1, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)
    print('Spatial+ Fixed df 10')
    spatplus.fx10 <- spat.plus.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k = 11, 
                                   fx = 1, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT)  
  }else{
    spatplus <- list("beta" = NA,
                     "se" = NA,
                     "coverage" = NA,
                     "k.used.e" = NA,
                     "k.used.r" = NA,
                     "reject" = NA)
    spatplus.fx <- list("beta" = NA,
                        "se" = NA,
                        "coverage" = NA,
                        "k.used.e" = NA,
                        "k.used.r" = NA,
                        "reject" = NA)
    spatplus.fx10 <- list("beta" = NA,
                          "se" = NA,
                          "coverage" = NA,
                          "k.used.e" = NA,
                          "k.used.r" = NA,
                          "reject" = NA)
  }
  
  #Run KS models
  if(tprs.flag){
    print('KS AIC')
    tprs.adj.AIC <- tprs.adjust.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, 
                                    AIC = TRUE, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT, max.df = max.df, fx = F, df = 0)
    print('KS BIC')
    tprs.adj.BIC <- tprs.adjust.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, 
                                    AIC = FALSE, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT, max.df = max.df, fx = F, df = 0)
  }else{
    tprs.adj.AIC <- list("beta" = NA,
                         "se" = NA,
                         "coverage" = NA,
                         "k.used.e" = NA,
                         "k.used.r" = NA,
                         "reject" = NA)
    tprs.adj.BIC <- list("beta" = NA,
                         "se" = NA,
                         "coverage" = NA,
                         "k.used.e" = NA,
                         "k.used.r" = NA,
                         "reject" = NA)
  }
  
  if(tprs.fx.flag){
    print('KS 3')
    tprs.adj.df3 <- tprs.adjust.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, 
                                    AIC = TRUE, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT, max.df = max.df, fx = T, df = 3)
    print('KS 10')
    tprs.adj.df10 <- tprs.adjust.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, 
                                     AIC = TRUE, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT, max.df = max.df, fx = T, df = 10)
  }else{
    tprs.adj.df3 <- list("beta" = NA,
                         "se" = NA,
                         "coverage" = NA,
                         "k.used.e" = NA,
                         "k.used.r" = NA,
                         "reject" = NA)
    tprs.adj.df10 <- list("beta" = NA,
                          "se" = NA,
                          "coverage" = NA,
                          "k.used.e" = NA,
                          "k.used.r" = NA,
                          "reject" = NA)
  }
  
  #Run df-Spatial+ models
  if(knotspat.flag){
    print('Spatial+ AIC')
    knotspat.AIC <- knot.spat.plus.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, 
                                       AIC = TRUE, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT, max.df = max.df)
    print('Spatial+ BIC')
    knotspat.BIC <- knot.spat.plus.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, 
                                       AIC = FALSE, beta = beta, LOGISTIC_DAT = LOGISTIC_DAT, max.df = max.df)
  }else{
    knotspat.AIC <- list("beta" = NA,
                         "se" = NA,
                         "coverage" = NA,
                         "k.used.e" = NA,
                         "k.used.r" = NA,
                         "reject" = NA)
    knotspat.BIC <- list("beta" = NA,
                         "se" = NA,
                         "coverage" = NA,
                         "k.used.e" = NA,
                         "k.used.r" = NA,
                         "reject" = NA)
  }
  
  # Run E-PS models
  if(eps.flag){
    print('E-PS')
    eps <- eps.mod(x = test.data$x, y = test.data$y, samp.grid = test.grid, k = k, LOGISTIC_DAT = LOGISTIC_DAT, beta = beta)
  }else{
    eps <- list("beta" = NA,
                "se" = NA,
                "coverage" = NA,
                "k.used.e" = NA,
                "k.used.r" = NA,
                "reject" = NA)
  }
  
  lm.results[i,] <- c(lm$beta, lm$se, lm$coverage, lm$k.used.e, lm$k.used.r, lm$reject)
  gSEM.results[i,] <- c(gSEM$beta, gSEM$se, gSEM$coverage, gSEM$k.used.e, gSEM$k.used.r, gSEM$reject)
  gSEM.fx.results[i,] <- c(gSEM.fx$beta, gSEM.fx$se, gSEM.fx$coverage, gSEM.fx$k.used.e, gSEM.fx$k.used.r, gSEM.fx$reject)
  gSEM.fx10.results[i,] <- c(gSEM.fx10$beta, gSEM.fx10$se, gSEM.fx10$coverage, gSEM.fx10$k.used.e, gSEM.fx10$k.used.r, gSEM.fx10$reject)
  spatplus.results[i,] <- c(spatplus$beta, spatplus$se, spatplus$coverage, spatplus$k.used.e, spatplus$k.used.r, spatplus$reject)
  spatplus.fx.results[i,] <- c(spatplus.fx$beta, spatplus.fx$se, spatplus.fx$coverage, spatplus.fx$k.used.e,spatplus.fx$k.used.r, spatplus.fx$reject)
  spatplus.fx10.results[i,] <- c(spatplus.fx10$beta, spatplus.fx10$se, spatplus.fx10$coverage, spatplus.fx10$k.used.e,spatplus.fx10$k.used.r, spatplus.fx10$reject)
  tprs.AIC.results[i,] <- c(tprs.adj.AIC$beta, tprs.adj.AIC$se, tprs.adj.AIC$coverage, tprs.adj.AIC$k.used.e, tprs.adj.AIC$k.used.r, tprs.adj.AIC$reject)
  tprs.BIC.results[i,] <- c(tprs.adj.BIC$beta, tprs.adj.BIC$se, tprs.adj.BIC$coverage, tprs.adj.BIC$k.used.e, tprs.adj.BIC$k.used.r, tprs.adj.BIC$reject)
  tprs.df3.results[i,] <- c(tprs.adj.df3$beta, tprs.adj.df3$se, tprs.adj.df3$coverage, tprs.adj.df3$k.used.e, tprs.adj.df3$k.used.r, tprs.adj.df3$reject)
  tprs.df10.results[i,] <- c(tprs.adj.df10$beta, tprs.adj.df10$se, tprs.adj.df10$coverage, tprs.adj.df10$k.used.e, tprs.adj.df10$k.used.r, tprs.adj.df10$reject)
  knotspat.AIC.results[i,] <- c(knotspat.AIC$beta, knotspat.AIC$se, knotspat.AIC$coverage, knotspat.AIC$k.used.e, knotspat.AIC$k.used.r, knotspat.AIC$reject)
  knotspat.BIC.results[i,] <- c(knotspat.BIC$beta, knotspat.BIC$se, knotspat.BIC$coverage, knotspat.BIC$k.used.e, knotspat.BIC$k.used.r, knotspat.BIC$reject)
  eps.results[i,] <- c(eps$beta, eps$se, eps$coverage, eps$k.used.e, eps$k.used.r, eps$reject)
  
  seed = seed + 1
}

## Format Final Results Table ##
final.results <- matrix(NA, ncol = 8, nrow = 14)

final.results[1,] <- c(mean(lm.results[,1]), mean(lm.results[,2]), sd(lm.results[,1]), sqrt(mean((lm.results[,1]-beta)^2)), mean(lm.results[,3]), median(lm.results[,4]), median(lm.results[,5]), mean(lm.results[,6]))
final.results[2,] <- c(mean(gSEM.results[,1]), mean(gSEM.results[,2]), sd(gSEM.results[,1]), sqrt(mean((gSEM.results[,1]-beta)^2)), mean(gSEM.results[,3]), median(gSEM.results[,4]), median(gSEM.results[,5]), mean(gSEM.results[,6]))
final.results[3,] <- c(mean(gSEM.fx.results[,1]), mean(gSEM.fx.results[,2]), sd(gSEM.fx.results[,1]), sqrt(mean((gSEM.fx.results[,1]-beta)^2)), mean(gSEM.fx.results[,3]), median(gSEM.fx.results[,4]), median(gSEM.fx.results[,5]), mean(gSEM.fx.results[,6]))
final.results[4,] <- c(mean(gSEM.fx10.results[,1]), mean(gSEM.fx10.results[,2]), sd(gSEM.fx10.results[,1]), sqrt(mean((gSEM.fx10.results[,1]-beta)^2)), mean(gSEM.fx10.results[,3]), median(gSEM.fx10.results[,4]), median(gSEM.fx10.results[,5]), mean(gSEM.fx10.results[,6]))
final.results[5,] <- c(mean(spatplus.results[,1]), mean(spatplus.results[,2]), sd(spatplus.results[,1]), sqrt(mean((spatplus.results[,1]-beta)^2)), mean(spatplus.results[,3]), median(spatplus.results[,4]), median(spatplus.results[,5]), mean(spatplus.results[,6]))
final.results[6,] <- c(mean(spatplus.fx.results[,1]), mean(spatplus.fx.results[,2]), sd(spatplus.fx.results[,1]), sqrt(mean((spatplus.fx.results[,1]-beta)^2)), mean(spatplus.fx.results[,3]), median(spatplus.fx.results[,4]), median(spatplus.fx.results[,5]), mean(spatplus.fx.results[,6]))
final.results[7,] <- c(mean(spatplus.fx10.results[,1]), mean(spatplus.fx10.results[,2]), sd(spatplus.fx10.results[,1]), sqrt(mean((spatplus.fx10.results[,1]-beta)^2)), mean(spatplus.fx10.results[,3]), median(spatplus.fx10.results[,4]), median(spatplus.fx10.results[,5]), mean(spatplus.fx10.results[,6]))
final.results[8,] <- c(mean(tprs.AIC.results[,1]), mean(tprs.AIC.results[,2]), sd(tprs.AIC.results[,1]), sqrt(mean((tprs.AIC.results[,1]-beta)^2)), mean(tprs.AIC.results[,3]), median(tprs.AIC.results[,4]), median(tprs.AIC.results[,5]), mean(tprs.AIC.results[,6]))
final.results[9,] <- c(mean(tprs.BIC.results[,1]), mean(tprs.BIC.results[,2]), sd(tprs.BIC.results[,1]), sqrt(mean((tprs.BIC.results[,1]-beta)^2)), mean(tprs.BIC.results[,3]), median(tprs.BIC.results[,4]), median(tprs.BIC.results[,5]), mean(tprs.BIC.results[,6]))
final.results[10,] <- c(mean(tprs.df3.results[,1]), mean(tprs.df3.results[,2]), sd(tprs.df3.results[,1]), sqrt(mean((tprs.df3.results[,1]-beta)^2)), mean(tprs.df3.results[,3]), median(tprs.df3.results[,4]), median(tprs.df3.results[,5]), mean(tprs.df3.results[,6]))
final.results[11,] <- c(mean(tprs.df10.results[,1]), mean(tprs.df10.results[,2]), sd(tprs.df10.results[,1]), sqrt(mean((tprs.df10.results[,1]-beta)^2)), mean(tprs.df10.results[,3]), median(tprs.df10.results[,4]), median(tprs.df10.results[,5]), mean(tprs.df10.results[,6]))
final.results[12,] <- c(mean(knotspat.AIC.results[,1]), mean(knotspat.AIC.results[,2]), sd(knotspat.AIC.results[,1]), sqrt(mean((knotspat.AIC.results[,1]-beta)^2)), mean(knotspat.AIC.results[,3]), median(knotspat.AIC.results[,4]), median(knotspat.AIC.results[,5]), mean(knotspat.AIC.results[,6]))
final.results[13,] <- c(mean(knotspat.BIC.results[,1]), mean(knotspat.BIC.results[,2]), sd(knotspat.BIC.results[,1]), sqrt(mean((knotspat.BIC.results[,1]-beta)^2)), mean(knotspat.BIC.results[,3]), median(knotspat.BIC.results[,4]), median(knotspat.BIC.results[,5]), mean(knotspat.BIC.results[,6]))
final.results[14,] <- c(mean(eps.results[,1]), mean(eps.results[,2]), sd(eps.results[,1]), sqrt(mean((eps.results[,1]-beta)^2)), mean(eps.results[,3]), median(eps.results[,4]), median(eps.results[,5]), mean(eps.results[,6]))

rownames(final.results) <- c("linear", "gSEM GCV", "gSEM df = 300", "gSEM df = 10", "Spatial+ GCV", "Spatial+ df = 300", "Spatial+ df = 10", "tprs AIC", "tprs BIC", "tprs df3", "tprs df10", "AIC Spatial+", "BIC Spatial+", "E-PS")
colnames(final.results) <- c("Estimated Beta", "Averaged Standard Errors", "SD of Estimates", "RMSE", "Coverage Rate", "Med. Knot Used (Exp)", "Med. Knot Used (Resp)", "Power")

final.results
