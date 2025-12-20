
library(tidyverse)

# these scripts contain necessary functions 
source("data-analysis-code/functions_analysis.R")

# "true" beta
sim_beta <- matrix(c(.2, .3, .4, .5, .3, .4, .5, .6, .4, .5, .2, .3, .3, .5, .2, .4, .5, .6, .1, .5), nrow = 4, ncol = 5)

# read in data 

Y_star <- readRDS("data/simulated/Y_star.rds")
X_star <- readRDS("data/simulated/X_star.rds")
Z_star <- readRDS("data/simulated/Z_star.rds")

Qmat_nc <- readRDS("data/simulated/Qmat_nc.rds")
Gamma_nc <- readRDS("data/simulated/Gamma_nc.rds")
W_nc <- readRDS("data/simulated/W_nc.rds")

exposure_names <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
outcome_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

set.seed(919)

B10_2 <-  bSpline(rank(W_nc)/length(W_nc),df=10, intercept = T)

iters <- 10000
burn <- 2000
thin <- 10
n_sample= (iters-burn)/thin

# MSM (full model)

fit_svi_K10_L10_stB <- MCMC_model4_8_3_XZ(Y_star, X_star, Z_star, Qmat_nc, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=10, B=B10_2, iters = iters, burn = burn, W_nc, thin = thin)
saveRDS(fit_svi_K10_L10_stB, file = "data-analysis-code/simulated_result/fit_svi_K10_L10_stB.rds")
# fit_svi_K10_L10_stB <- readRDS("data-analysis-code/simulated_result/fit_svi_K10_L10_stB.rds")


### naive model 

fit_svi_naive <- MCMC_model3_2_XZ(Y_star, X_star, Z_star, Qmat_nc, p=4, n=nrow(Y_star), R=5, Q=5, iters=iters, burn=burn, 1, nrow(Y_star), W_nc, thin = thin)
saveRDS(fit_svi_naive, file = "data-analysis-code/simulated_result/fit_svi_naive.rds")
# fit_svi_naive <- readRDS("data-analysis-code/simulated_result/fit_svi_naive.rds")

### Calculate coefficient from tensor margins

E <- 4
R <- 5
Model <- fit_svi_K10_L10_stB
n <- nrow(Y_star)
Beta <- array(0, c(n,E,R,n_sample))
for (i in 1:n_sample){
  Beta[,,,i] <- bbygamma(B10_2, find_tol_gamma(Model$Tl[,,i], Model$Te[,,i], Model$Tr[,,i]), n) 
}

mean <-  rowMeans(Beta[,,,], dims = 3)
lower <- upper <- array(NA, c(nrow(Y_star),E,R))
for (freq in 1:nrow(Y_star)){
  for (e in 1:E){
    for (r in 1:R){
      lower[freq,e,r] <- quantile(Beta[freq,e,r,], probs = .025)
      upper[freq,e,r] <- quantile(Beta[freq,e,r,], probs = .975)
    }
  }
}

prob <- confounding <- matrix(NA, nrow = E, ncol = R)
for (e in 1:E){
  for (r in 1:R){
    prob[e,r] <- mean(Beta[1,e,r,]>0)
    confounding[e,r] <- mean(Beta[n,e,r,]>Beta[1,e,r,])
  }
}
rownames(prob) <- exposure_names
colnames(prob) <- outcome_names

# The object prob make a table for Percentages of posterior estimates being greater than zero (like Table 1 in the paper)

prob


# The object confounding make a table for Probabilities of coefficient estimates at the lowest frequency greater than those at the highest frequency (like Table 2 in the paper)

confounding

# check result for all frequencies
par(mfrow = c(E,R), mar = c(2,2,2,2))
for (e in 1:E){
  for (r in 1:R){
    plot(W_nc, mean[,e,r], type = "l", main = paste0(exposure_names[e],", ", outcome_names[r]), ylim = c(min(lower[,e,r]),max(upper[,e,r])))
    lines(W_nc, upper[,e,r], col = "red")
    lines(W_nc, lower[,e,r], col = "red")
    abline(h=sim_beta[e,r], col = "green")
  }
}

# trace plots
par(mfrow = c(E,R), mar = c(2,2,2,2))
for (e in 1:E){
  for (r in 1:R){
    plot(Beta[1,e,r,], type = "l", main = paste0(exposure_names[e],", ", outcome_names[r]))
    abline(h = mean[1,e,r], col = "yellow")
    abline(h = upper[1,e,r], col = "red")
    abline(h = lower[1,e,r], col = "red")
    abline(h = sim_beta[e,r], col = "green")
  }
}

### OLS

fit_svi <- list()
for (r in 1:R){
  fit_svi[[r]] <- lm(Y_star[,r] ~ X_star + Z_star - 1)
  print(summary(fit_svi[[r]]))
}

### naive

# fit_svi_naive <- readRDS("fit_svi_naive.rds")

model_naive <- fit_svi_naive
mean_naive <-  upper_naive <- lower_naive <- matrix(NA, nrow = E, ncol = R)
for (e in 1:E){
  for (r in 1:R){
    mean_naive[e,r] <- mean(model_naive$beta[e,r,])
    lower_naive[e,r] <- quantile(model_naive$beta[e,r,], probs = .025)
    upper_naive[e,r] <- quantile(model_naive$beta[e,r,], probs = .975)
  }
}

rownames(mean_naive) <- rownames(lower_naive) <- rownames(upper_naive) <- exposure_names
colnames(mean_naive) <- colnames(lower_naive) <- colnames(upper_naive) <- outcome_names

######################################
########### make plots ###############
######################################

n <- nrow(X_star)
model_type = c("MSM", "Naive", "OLS", "Truth")
model_level <- c(paste0(exposure_names[1], ", ", model_type),
                 paste0(exposure_names[2], ", ", model_type),
                 paste0(exposure_names[3], ", ", model_type),
                 paste0(exposure_names[4], ", ", model_type))

sum <- list()
for (r in 1:R){
  sum_msm_high <- sum_naive <- sum_OLS <- sum_OLS_state <- truth <- list()
  for (e in 1:E){
    sum_msm_high[[e]] <- data.frame(
      Exposure = exposure_names[e], 
      Outcome = outcome_names[r],
      Model = paste0("MSM"),
      EffectSize = mean[1,e,r], 
      CI_lower = lower[1,e,r], 
      CI_upper = upper[1,e,r])
    
    sum_naive[[e]] <- data.frame(
      Exposure = exposure_names[e],
      Outcome = outcome_names[r],
      Model = paste0("Naive"),
      EffectSize = mean_naive[e,r], 
      CI_lower = lower_naive[e,r], 
      CI_upper = upper_naive[e,r])
    
    sum_OLS[[e]] <- data.frame(
      Exposure = exposure_names[e], 
      Outcome = outcome_names[r],
      Model = paste0("OLS"),
      EffectSize = summary(fit_svi[[r]])$coefficients[e,1],
      CI_lower = summary(fit_svi[[r]])$coefficients[e,1]-2*summary(fit_svi[[r]])$coefficients[e,2],
      CI_upper = summary(fit_svi[[r]])$coefficients[e,1]+2*summary(fit_svi[[r]])$coefficients[e,2]) 
    
    truth[[e]] <- data.frame(
      Exposure = exposure_names[e], 
      Outcome = outcome_names[r],
      Model = paste0("Truth"),
      EffectSize = sim_beta[e,r],
      CI_lower = sim_beta[e,r],
      CI_upper = sim_beta[e,r]) 
  }
  
  sum[[r]] <- rbind(
    sum_msm_high[[1]],
    sum_naive[[1]],
    sum_OLS[[1]],
    sum_msm_high[[2]],
    sum_naive[[2]],
    sum_OLS[[2]],
    sum_msm_high[[3]],
    sum_naive[[3]],
    sum_OLS[[3]],
    sum_msm_high[[4]],
    sum_naive[[4]],
    sum_OLS[[4]],
    truth[[1]],
    truth[[2]],
    truth[[3]],
    truth[[4]])
  
  sum[[r]] <- sum[[r]] %>% arrange(factor(Model, levels = model_level))       
}

full_sum <- rbind(sum[[1]],sum[[2]],sum[[3]],sum[[4]],sum[[5]])

ggplot(data = full_sum, aes(x = factor(Outcome, level=outcome_names), y = EffectSize, ymin = CI_lower, ymax = CI_upper)) +
  geom_pointrange(position=position_dodge(width=.5), size = .2, aes(color = factor(Model, levels = model_type))) + 
  scale_color_manual(name='Model', labels=c('MSM', 'Naive', 'OLS', 'Truth'), values = c("#440154FF", "#2A788EFF", "#FDE725FF", "red"))+
  labs(color = "Model", x = "", y = "Effect Size") + 
  # geom_hline(yintercept = 0, linetype = 3) +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(factor(Exposure,level=exposure_names)), nrow = 2) 
  # scale_colour_viridis_d()
ggsave("data-analysis-code/simulated_plots/forest.png")

# coefficient by frequency plots 

curr_data <- data.frame(W = W_nc, Mean = mean[,1,1], CI_upper = upper[,1,1], CI_lower = lower[,1,1], Exposure = exposure_names[1], Outcome = outcome_names[1], Model = "MSM")

for (r in 2:5){
  curr_data2 <- data.frame(W = W_nc, Mean = mean[,1,r], CI_upper = upper[,1,r], CI_lower = lower[,1,r], Exposure = exposure_names[1], Outcome = outcome_names[r], Model = "MSM")
  curr_data <- rbind(curr_data2, curr_data)
}

for (e in 2:E){
  for (r in 1:5){
    curr_data3 <- data.frame(W = W_nc, Mean = mean[,e,r], CI_upper = upper[,e,r], CI_lower = lower[,e,r], Exposure = exposure_names[e], Outcome = outcome_names[r], Model = "MSM")
    curr_data <- rbind(curr_data3, curr_data)
  }
}

for (e in 1:E){
  for (r in 1:R){
    curr_data_naive <- data.frame(
      W = W_nc, 
      Mean = mean_naive[e,r], 
      CI_upper = upper_naive[e,r],
      CI_lower = lower_naive[e,r], 
      Exposure = exposure_names[e],
      Outcome = outcome_names[r],
      Model = paste0("Naive")
    )
    curr_data <- rbind(curr_data_naive, curr_data)
  }
}

for (e in 1:E){
  for (r in 1:R){
    curr_data_ols <- data.frame(
      W = W_nc,
      Mean = summary(fit_svi[[r]])$coefficients[e,1],
      CI_upper = summary(fit_svi[[r]])$coefficients[e,1]+2*summary(fit_svi[[r]])$coefficients[e,2],
      CI_lower = summary(fit_svi[[r]])$coefficients[e,1]-2*summary(fit_svi[[r]])$coefficients[e,2],
      Exposure = exposure_names[e], 
      Outcome = outcome_names[r],
      Model = paste0("OLS")
    ) 
    curr_data <- rbind(curr_data_ols, curr_data)
  }
}

for (e in 1:E){
  for (r in 1:R){
    curr_data_true <- data.frame(
      W = W_nc,
      Mean = sim_beta[e,r],
      CI_upper = sim_beta[e,r],
      CI_lower = sim_beta[e,r],
      Exposure = exposure_names[e], 
      Outcome = outcome_names[r],
      Model = paste0("True")
    ) 
    curr_data <- rbind(curr_data_true, curr_data)
  }
}



curr_data_long <- curr_data %>% pivot_longer(cols = c("Mean","CI_upper","CI_lower"), names_to = "line_type", values_to = "estimate")

ggplot(data = curr_data_long, aes(x=W, y=estimate, linetype = line_type, color = Model)) + 
  geom_line() + 
  scale_color_manual(name='Model', labels=c('MSM', 'Naive', 'OLS', 'Truth'), values = c("#440154FF", "#2A788EFF", "#FDE725FF", "red")) +
  labs(x="Eigenvalue",y="Effect Size") + 
  scale_linetype_manual(values = c("Mean" = "solid", "CI_upper" = "dotted", "CI_lower" = "dotted")) +
  geom_hline(yintercept = 0, linetype = 3) + 
  facet_grid(cols=vars(factor(Outcome,level=outcome_names)),rows=vars(factor(Exposure,level=exposure_names)), scales = "free") +
  theme(strip.text = element_text(size=8),
        axis.text = element_text(size=7),
        legend.position = "right") + 
  # geom_hline(yintercept = , linetype = 3) +
  guides(linetype = "none")
ggsave("data-analysis-code/simulated_plots/freq.png")
