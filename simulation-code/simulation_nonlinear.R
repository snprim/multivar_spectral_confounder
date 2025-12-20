#############################################
##### Author: Shih-Ni Prim ##################
#############################################

# This is the additional simulations for nonlinear confounding
# 500 datasets 

# MCMC_model4_8 full model with HS prior
# MCMC_model4_8_2 full model with IG prior
# MCMC_model4_9 univariate with HS prior
# MCMC_model4_9_2 univariate with IG prior
# MCMC_model3_2 does not allow coefficients to vary by frequency

args = commandArgs(TRUE)
ind = as.numeric(args[1])
start_i = as.numeric(args[2])

source("functions_sim.R")

iters <- 5000
burn <- 1000
sig2X <- c(3,3,3,3,3,3,3,3,3)
sig2Z <- rep(16,10)
sig2theta <- c(2,2,2,2,2)
tau2r <- c(4,4,4,4,4)
lambdaZ <- .9999
lambda1 <- .9
lambda2 <- .9
rho1 <- .9
rho2 <- .9
n <- 400
n1 <- 20
E <- 10
R <- 5
nz <- 10

L <- 10
K <- 5
b <- 8

Qmat <- makeQ(20)
W <- eigen(Qmat)$values

beta <- matrix(c(3,-3,2,3,5,4,-2,4,5,7,2,-3,3,4,6,1,-4,2,3,4,5,1,3,4,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow = E, ncol = R, byrow = T)/10

reps <- 100

load("simulation/inputs/B1.RData")
load("simulation/inputs/B2.RData")

###########################################################
################ main simulation ##########################
###########################################################


# extra smoothing, bw = 1, bXZ = 1
# k=0.5
if (ind == 1){
  sim1_par <- foreach (i=start_i:(start_i+reps-1), .combine = rbind) %dopar% {
    list(sim_func_nonlinear(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters=20000, burn=5000, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, k=0.5, b=b, i))
  }
  nonlinear_sim1 <- flatten(sim1_par, reps, 2, n, 7, E, R)
  filename <- paste0("nonlinear_sim1_",start_i,".RData")
  save(nonlinear_sim1, file = filename)
}


# extra smoothing, bw = 1, bXZ = 1
# k=1
if (ind == 2){
  sim2_par <- foreach (i=start_i:(start_i+reps-1), .combine = rbind) %dopar% {
    list(sim_func_nonlinear(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters=20000, burn=5000, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, k=1, b=b, i))
  }
  nonlinear_sim2 <- flatten(sim2_par, reps, 2, n, 7, E, R)
  filename <- paste0("nonlinear_sim2_",start_i,".RData")
  save(nonlinear_sim2, file = filename)
}


# extra smoothing, bw = 1, bXZ = 1
# k=1.5
if (ind == 3){
  sim3_par <- foreach (i=start_i:(start_i+reps-1), .combine = rbind) %dopar% {
    list(sim_func_nonlinear(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters=20000, burn=5000, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, k=1.5, b=b, i))
  }
  nonlinear_sim3 <- flatten(sim3_par, reps, 2, n, 7, E, R)
  filename <- paste0("nonlinear_sim3_",start_i,".RData")
  save(nonlinear_sim3, file = filename)
}


# extra smoothing, bw = 1, bXZ = 1
# k=2
if (ind == 4){
  sim4_par <- foreach (i=start_i:(start_i+reps-1), .combine = rbind) %dopar% {
    list(sim_func_nonlinear(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters=20000, burn=5000, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, k=2, b=b, i))
  }
  nonlinear_sim4 <- flatten(sim4_par, reps, 2, n, 7, E, R)
  filename <- paste0("nonlinear_sim4_",start_i,".RData")
  save(nonlinear_sim4, file = filename)
}


# extra smoothing, bw = 1, bXZ = 1
# k=2.5
if (ind == 5){
  sim5_par <- foreach (i=start_i:(start_i+reps-1), .combine = rbind) %dopar% {
    list(sim_func_nonlinear(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters=20000, burn=5000, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, k=2.5, b=b, i))
  }
  nonlinear_sim5 <- flatten(sim5_par, reps, 2, n, 7, E, R)
  filename <- paste0("nonlinear_sim5_",start_i,".RData")
  save(nonlinear_sim5, file = filename)
}


# extra smoothing, bw = 1, bXZ = 1
# k=3
if (ind == 6){
  sim6_par <- foreach (i=start_i:(start_i+reps-1), .combine = rbind) %dopar% {
    list(sim_func_nonlinear(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters=20000, burn=5000, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, k=3, b=b, i))
  }
  nonlinear_sim6 <- flatten(sim6_par, reps, 2, n, 7, E, R)
  filename <- paste0("nonlinear_sim6_",start_i,".RData")
  save(nonlinear_sim6, file = filename)
}
