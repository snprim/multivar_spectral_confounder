
####################################################################################
### gather data first in data_preparation.R (if you have access to Medicare data) ## 
### data_preparation_synthetic.R (if you don't have access to Medicare data) #######
####################################################################################

library(tidyverse)

# these scripts contain necessary functions 
source("data-analysis-code/functions_analysis.R")
source("data-analysis-code/cross_validation.R")

# read in data 

# Y_star <- readRDS("Y_star.rds")
# X_star <- readRDS("X_star.rds")
# Z_star <- readRDS("Z_star.rds")
# 
# Qmat_south <- readRDS("Qmat_south.rds")
# Gamma_south <- readRDS("Gamma_south.rds")
# W_south <- readRDS("W_south.rds")

exposure_names <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
outcome_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

set.seed(919)

B5_2 <-  bSpline(rank(W_south)/length(W_south),df=5, intercept = T)
B10_2 <-  bSpline(rank(W_south)/length(W_south),df=10, intercept = T)

# use standardized basis function
# exposures are 4 svi
# covariates are log(population) and urbanicity
# final model

iters <- 100000
burn <- 10000
thin <- 10
n_sample= (iters-burn)/thin


# this takes about 2 hours for 5000 iterations on a personal laptop
fit_svi_K10_L10_stB <- MCMC_model4_8_3_XZ(Y_star, X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)

saveRDS(fit_svi_K10_L10_stB, file = "fit_svi_K10_L10_stB.rds")
# fit_svi_K10_L10_stB <- readRDS("fit_svi_K10_L10_stB.rds")

# make trace plot

E <- 4
R <- 5
n <- nrow(Y_star)
final_beta <- array(0, c(n,E,R,n_sample))
for (i in 1:n_sample){
  final_beta[,,,i] <- bbygamma(B10_2, find_tol_gamma(fit_svi_K10_L10_stB$Tl[,,i], fit_svi_K10_L10_stB$Te[,,i], fit_svi_K10_L10_stB$Tr[,,i]), n) 
}

png("plots/trace_final.png", width=6, height=4, units="in", res=300)
par(mfrow = c(4, 5), mar = c(2,2,2,2))
for (i in 1:4){
  for (j in 1:5){
    plot(final_beta[1,i,j,], type = "l")
  }
}
dev.off()


# find ESS
ESS <- matrix(NA, nrow = 4, ncol = 5)
for (i in 1:4){
  for (j in 1:5){
    ESS[i,j] <- effectiveSize(final_beta[1,i,j,])
  }
}
ESS


### naive model 

# final naive model
# takes about 40 minutes for 5000 iterations
fit_svi_naive <- MCMC_model3_2_XZ(Y_star, X_star, Z_star, Qmat_south, p=4, n=nrow(Y_star), R=5, Q=5, iters=iters, burn=burn, 1, nrow(Y_star), W_south, thin = thin)

saveRDS(fit_svi_naive, file = "fit_svi_naive.rds")

#################################################################
########## make plots at final_results_plots.R ##################
#################################################################
