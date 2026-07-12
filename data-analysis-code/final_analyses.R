
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

Qmat_south <- readRDS("data/outputs/Qmat_south.rds")
# eig_south <- eigen(Qmat_south)
Gamma_south <- readRDS("data/outputs/Gamma_south.rds")
W_south <- readRDS("data/outputs/W_south.rds")

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
# fit_svi_K10_L10_stB <- MCMC_model4_8_3_XZ(Y_star, X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)

# saveRDS(fit_svi_K10_L10_stB, file = "fit_svi_K10_L10_stB.rds")
# fit_svi_K10_L10_stB <- readRDS("fit_svi_K10_L10_stB.rds")

# make trace plot
# 
# E <- 4
# R <- 5
# n <- nrow(Y_star)
# final_beta <- array(0, c(n,E,R,n_sample))
# for (i in 1:n_sample){
#   final_beta[,,,i] <- bbygamma(B10_2, find_tol_gamma(fit_svi_K10_L10_stB$Tl[,,i], fit_svi_K10_L10_stB$Te[,,i], fit_svi_K10_L10_stB$Tr[,,i]), n) 
# }
# 
# png("plots/trace_final.png", width=6, height=4, units="in", res=300)
# par(mfrow = c(4, 5), mar = c(2,2,2,2))
# for (i in 1:4){
#   for (j in 1:5){
#     plot(final_beta[1,i,j,], type = "l")
#   }
# }
# dev.off()

################################################
##### Rerun the analysis to save theta-hat #####
################################################

fit_svi_K10_L10 <- MCMC_model4_8_4_XZ(Y_star, X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)

saveRDS(fit_svi_K10_L10, file = "fit_svi_K10_L10.rds")

# fit_svi_K10_L10 <- readRDS("results/fit_svi_K10_L10.rds")

n <- 10149
E <- 4
R <- 5
final_beta <- array(0, c(n,E,R,n_sample))
for (i in 1:n_sample){
  final_beta[,,,i] <- bbygamma(B10_2, find_tol_gamma(fit_svi_K10_L10$Tl[,,i], fit_svi_K10_L10$Te[,,i], fit_svi_K10_L10$Tr[,,i]), n) 
}

for (i in 1:4){
  for (j in 1:5){
    plot(final_beta[1,i,j,], type = "l")
  }
}

for (i in 1:5){
  plot(fit_svi_K10_L10$theta[1,i,], type = "l", ylab = expression(theta), main = paste0("Theta-hat of ", i, "th Outcome"))  
}

round(rowMeans(fit_svi_K10_L10$theta[1,,]),3)
# [1] -0.098  0.045 -0.153 -0.041  0.124
round(rowMeans(final_beta[1,,,], dims = 2),3)
#        [,1]   [,2]   [,3]   [,4]   [,5]
# [1,]  0.042  0.296 -0.472  0.242  0.535
# [2,]  0.076  0.078  0.017  0.173  0.036
# [3,] -0.047 -0.104 -0.147 -0.148  0.140
# [4,] -0.011  0.146 -0.059  0.253 -0.059

####################################################
##### Rerun the analysis for Univariate models #####
####################################################

# fit_svi_K10_L10_univ1 <- MCMC_model4_9_XZ(Y_star[,1], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=1, Q=1, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)
# 
# saveRDS(fit_svi_K10_L10_univ1, file = "fit_svi_K10_L10_univ1.rds")
# 
# fit_svi_K10_L10_univ2 <- MCMC_model4_9_XZ(Y_star[,2], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=1, Q=1, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)
# 
# saveRDS(fit_svi_K10_L10_univ2, file = "fit_svi_K10_L10_univ2.rds")
# 
# fit_svi_K10_L10_univ3 <- MCMC_model4_9_XZ(Y_star[,3], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=1, Q=1, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)
# 
# saveRDS(fit_svi_K10_L10_univ3, file = "fit_svi_K10_L10_univ3.rds")
# 
# fit_svi_K10_L10_univ4 <- MCMC_model4_9_XZ(Y_star[,4], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=1, Q=1, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)
# 
# saveRDS(fit_svi_K10_L10_univ4, file = "fit_svi_K10_L10_univ4.rds")
# 
# fit_svi_K10_L10_univ5 <- MCMC_model4_9_XZ(Y_star[,5], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=1, Q=1, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_south, thin = thin)
# 
# saveRDS(fit_svi_K10_L10_univ5, file = "fit_svi_K10_L10_univ5.rds")


fit_svi_K10_L10_univ1 <- readRDS("results/fit_svi_K10_L10_univ1.rds")
fit_svi_K10_L10_univ2 <- readRDS("results/fit_svi_K10_L10_univ2.rds")
fit_svi_K10_L10_univ3 <- readRDS("results/fit_svi_K10_L10_univ3.rds")
fit_svi_K10_L10_univ4 <- readRDS("results/fit_svi_K10_L10_univ4.rds")
fit_svi_K10_L10_univ5 <- readRDS("results/fit_svi_K10_L10_univ5.rds")

# calculate beta 

# univ1
final_beta_univ1 <- array(0, c(n,E,n_sample))
for (i in 1:n_sample){
  tol_gamma <- find_tol_gamma3(fit_svi_K10_L10_univ1$Tl[,,i], fit_svi_K10_L10_univ1$Te[,,i])
  final_beta_univ1[,,i] <- B10_2 %*% tol_gamma
}

round(rowMeans(final_beta_univ1[1,,]),3)

# univ2
final_beta_univ2 <- array(0, c(n,E,n_sample))
for (i in 1:n_sample){
  tol_gamma <- find_tol_gamma3(fit_svi_K10_L10_univ2$Tl[,,i], fit_svi_K10_L10_univ2$Te[,,i])
  final_beta_univ2[,,i] <- B10_2 %*% tol_gamma
}

# univ3
final_beta_univ3 <- array(0, c(n,E,n_sample))
for (i in 1:n_sample){
  tol_gamma <- find_tol_gamma3(fit_svi_K10_L10_univ3$Tl[,,i], fit_svi_K10_L10_univ3$Te[,,i])
  final_beta_univ3[,,i] <- B10_2 %*% tol_gamma
}

# univ4
final_beta_univ4 <- array(0, c(n,E,n_sample))
for (i in 1:n_sample){
  tol_gamma <- find_tol_gamma3(fit_svi_K10_L10_univ4$Tl[,,i], fit_svi_K10_L10_univ4$Te[,,i])
  final_beta_univ4[,,i] <- B10_2 %*% tol_gamma
}

# univ5
final_beta_univ5 <- array(0, c(n,E,n_sample))
for (i in 1:n_sample){
  tol_gamma <- find_tol_gamma3(fit_svi_K10_L10_univ5$Tl[,,i], fit_svi_K10_L10_univ5$Te[,,i])
  final_beta_univ5[,,i] <- B10_2 %*% tol_gamma
}


cbind(round(rowMeans(final_beta_univ1[1,,]),3),
      round(rowMeans(final_beta_univ2[1,,]),3),
      round(rowMeans(final_beta_univ3[1,,]),3),
      round(rowMeans(final_beta_univ4[1,,]),3),
      round(rowMeans(final_beta_univ5[1,,]),3)
      )
#        [,1]   [,2]   [,3]   [,4]   [,5]
# [1,] -0.092  0.224 -0.650  0.023  0.290
# [2,]  0.071 -0.124 -0.119  0.232 -0.030
# [3,] -0.224 -0.229 -0.251 -0.211  0.028
# [4,]  0.059  0.231 -0.016  0.230 -0.034

cbind(round(mean(fit_svi_K10_L10_univ1$U[1,1,]),3),
      round(mean(fit_svi_K10_L10_univ2$U[1,1,]),3),
      round(mean(fit_svi_K10_L10_univ3$U[1,1,]),3),
      round(mean(fit_svi_K10_L10_univ4$U[1,1,]),3),
      round(mean(fit_svi_K10_L10_univ5$U[1,1,]),3))
#        [,1]  [,2]   [,3]   [,4]  [,5]
# [1,] -0.129 0.074 -0.237 -0.066 0.213

for (e in 1:4){
  plot(final_beta_univ1[1,e,], type = "l")
}

plot(fit_svi_K10_L10_univ1$U[1,1,], type="l")


#############################################
############### find ESS ####################
#############################################

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
