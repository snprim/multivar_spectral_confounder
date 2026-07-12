
#################################################################
############### This file is for make plots. ####################   
############### Run final_analyses.R first ######################
#################################################################


library(tidyverse)

# these scripts contain necessary functions 
source("functions_analysis.R")
source("cross_validation.R")

# read in data 

# Y_star <- readRDS("Y_star.rds")
# X_star <- readRDS("X_star.rds")
# Z_star <- readRDS("Z_star.rds")
# 
# Qmat_south <- readRDS("Qmat_south.rds")
# Gamma_south <- readRDS("Gamma_south.rds")
# W_south <- readRDS("W_south.rds")

iters <- 100000
burn <- 10000
thin <- 10
n_sample= (iters-burn)/thin

exposure_names <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
outcome_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

# read in results for MSM

# calculate beta from Tl, Te, Tr

# fit_svi_K10_L10_stB <- readRDS("fit_svi_K10_L10_stB.rds")
E <- 4
R <- 5
n <- 10149
Model <- fit_svi_K10_L10
n <- nrow(Y_star)
final_beta <- array(0, c(n,E,R,n_sample))
for (i in 1:n_sample){
  final_beta[,,,i] <- bbygamma(B10_2, find_tol_gamma(Model$Tl[,,i], Model$Te[,,i], Model$Tr[,,i]), n) 
}

mean <-  rowMeans(final_beta[,,,], dims = 3)
lower <- upper <- array(NA, c(n,E,R))
for (freq in 1:n){
  for (e in 1:E){
    for (r in 1:R){
      lower[freq,e,r] <- quantile(final_beta[freq,e,r,], probs = .025)
      upper[freq,e,r] <- quantile(final_beta[freq,e,r,], probs = .975)
    }
  }
}

saveRDS(mean, "data/outputs/mean.rds")
saveRDS(lower, "data/outputs/lower.rds")
saveRDS(upper, "data/outputs/upper.rds")

# prob <- confounding <- matrix(NA, nrow = E, ncol = R)
# for (e in 1:E){
#   for (r in 1:R){
#     prob[e,r] <- mean(Beta[1,e,r,]>0)
#     confounding[e,r] <- mean(Beta[n,e,r,]>Beta[1,e,r,])
#   }
# }
# rownames(prob) <- exposure_names
# colnames(prob) <- outcome_names

par(mfrow = c(E,R), mar = c(2,2,2,2))
for (e in 1:E){
  for (r in 1:R){
    plot(W_south, mean[,e,r], type = "l", main = paste0(exposure_names[e],", ", outcome_names[r]), ylim = c(min(lower[,e,r]),max(upper[,e,r])))
    lines(W_south, upper[,e,r], col = "red")
    lines(W_south, lower[,e,r], col = "red")
    abline(h=0, col = "green")
  }
}

par(mfrow = c(E,R), mar = c(2,2,2,2))
for (e in 1:E){
  for (r in 1:R){
    plot(Beta[1,e,r,], type = "l", main = paste0(exposure_names[e],", ", outcome_names[r]))
    abline(h = mean[1,e,r], col = "yellow")
    abline(h = upper[1,e,r], col = "red")
    abline(h = lower[1,e,r], col = "red")
    abline(h = 0, col = "green")
  }
}

### Univariate 

mean_univ <- upper_univ <- lower_univ <- array(0, c(n,E,R))

mean_univ[,,1] <- rowMeans(final_beta_univ1, dims = 2)
mean_univ[,,2] <- rowMeans(final_beta_univ2, dims = 2)
mean_univ[,,3] <- rowMeans(final_beta_univ3, dims = 2)
mean_univ[,,4] <- rowMeans(final_beta_univ4, dims = 2)
mean_univ[,,5] <- rowMeans(final_beta_univ5, dims = 2)

for (freq in 1:n){
  for (e in 1:E){
    lower_univ[freq,e,1] <- quantile(final_beta_univ1[freq,e,], probs = .025)
    upper_univ[freq,e,1] <- quantile(final_beta_univ1[freq,e,], probs = .975)
    lower_univ[freq,e,2] <- quantile(final_beta_univ2[freq,e,], probs = .025)
    upper_univ[freq,e,2] <- quantile(final_beta_univ2[freq,e,], probs = .975)
    lower_univ[freq,e,3] <- quantile(final_beta_univ3[freq,e,], probs = .025)
    upper_univ[freq,e,3] <- quantile(final_beta_univ3[freq,e,], probs = .975)
    lower_univ[freq,e,4] <- quantile(final_beta_univ4[freq,e,], probs = .025)
    upper_univ[freq,e,4] <- quantile(final_beta_univ4[freq,e,], probs = .975)
    lower_univ[freq,e,5] <- quantile(final_beta_univ5[freq,e,], probs = .025)
    upper_univ[freq,e,5] <- quantile(final_beta_univ5[freq,e,], probs = .975)
  }
}


saveRDS(mean_univ, "data/outputs/mean_univ.rds")
saveRDS(lower_univ, "data/outputs/lower_univ.rds")
saveRDS(upper_univ, "data/outputs/upper_univ.rds")


### OLS
# 
# fit_svi <- list()
# for (r in 1:R){
#   fit_svi[[r]] <- lm(Y_star[,r] ~ X_star + Z_star - 1)
#   print(summary(fit_svi[[r]]))
# }

fit_svi <- readRDS("results/fit_svi_ols.rds")

### naive

fit_svi_naive <- readRDS("results/fit_svi_naive.rds")

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

mean <- readRDS("data/outputs/mean.rds")
lower <- readRDS("data/outputs/lower.rds")
upper <- readRDS("data/outputs/upper.rds")

mean_univ <- readRDS("data/outputs/mean_univ.rds")
lower_univ <- readRDS("data/outputs/lower_univ.rds")
upper_univ <- readRDS("data/outputs/upper_univ.rds")



model_type = c("MSM", "USM", "Naive", "OLS")
model_level <- c(paste0(exposure_names[1], ", ", model_type),
                 paste0(exposure_names[2], ", ", model_type),
                 paste0(exposure_names[3], ", ", model_type),
                 paste0(exposure_names[4], ", ", model_type))

sum <- list()
for (r in 1:R){
  sum_msm_high <- sum_usm_high <- sum_naive <- sum_OLS <- list()
  for (e in 1:E){
    sum_msm_high[[e]] <- data.frame(
      Exposure = exposure_names[e], 
      Outcome = outcome_names[r],
      Model = paste0("MSM"),
      EffectSize = mean[1,e,r], 
      CI_lower = lower[1,e,r], 
      CI_upper = upper[1,e,r])
    
    sum_usm_high[[e]] <- data.frame(
      Exposure = exposure_names[e],
      Outcome = outcome_names[r],
      Model = paste0("USM"),
      EffectSize = mean_univ[1,e,r], 
      CI_lower = lower_univ[1,e,r], 
      CI_upper = upper_univ[1,e,r])
  
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
  }

  
  sum[[r]] <- rbind(
    sum_msm_high[[1]],
    sum_usm_high[[1]],
    sum_naive[[1]],
    sum_OLS[[1]],
    sum_msm_high[[2]],
    sum_usm_high[[2]],
    sum_naive[[2]],
    sum_OLS[[2]],
    sum_msm_high[[3]],
    sum_usm_high[[3]],
    sum_naive[[3]],
    sum_OLS[[3]],
    sum_msm_high[[4]],
    sum_usm_high[[4]],
    sum_naive[[4]],
    sum_OLS[[4]])
  
  sum[[r]] <- sum[[r]] %>% arrange(factor(Model, levels = model_level))       
}

full_sum <- rbind(sum[[1]],sum[[2]],sum[[3]],sum[[4]],sum[[5]])

ggplot(data = full_sum, aes(x = factor(Outcome, level=outcome_names), y = EffectSize, ymin = CI_lower, ymax = CI_upper)) +
  geom_pointrange(position=position_dodge(width=.5), size = .01, aes(color = factor(Model, levels = model_type)), linewidth = .2) + 
  scale_color_manual(name='Model', labels=c('MSM', 'USM', 'Naive', 'OLS'), values = c("MSM" = "#440154FF", "USM" = "#2A788EFF", "Naive" = "#CD3700", "OLS" = "#66CDAA")) +
  labs(color = "Outcome", x = "", y = "Effect Size") + 
  geom_hline(yintercept = 0, linetype = 3, linewidth = .4) +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(factor(Exposure,level=exposure_names)), nrow = 2) #+ 
  # scale_colour_viridis_d()

ggsave(filename = paste0("plots/MSM_USM_forest.png"), width=6, height=4, dpi = 300)

# coefficient by frequency plots 

curr_data <- data.frame(W = W_south, Mean = mean[,1,1], CI_upper = upper[,1,1], CI_lower = lower[,1,1], Exposure = exposure_names[1], Outcome = outcome_names[1], Model = "MSM")

for (r in 2:5){
  curr_data2 <- data.frame(W = W_south, Mean = mean[,1,r], CI_upper = upper[,1,r], CI_lower = lower[,1,r], Exposure = exposure_names[1], Outcome = outcome_names[r], Model = "MSM")
  curr_data <- rbind(curr_data2, curr_data)
}

for (e in 2:E){
  for (r in 1:5){
    curr_data3 <- data.frame(W = W_south, Mean = mean[,e,r], CI_upper = upper[,e,r], CI_lower = lower[,e,r], Exposure = exposure_names[e], Outcome = outcome_names[r], Model = "MSM")
    curr_data <- rbind(curr_data3, curr_data)
  }
}

for (e in 1:E){
  for (r in 1:R){
    curr_data_univ <- data.frame(
      W = W_south, 
      Mean = mean_univ[,e,r], 
      CI_upper = upper_univ[,e,r],
      CI_lower = lower_univ[,e,r], 
      Exposure = exposure_names[e],
      Outcome = outcome_names[r],
      Model = paste0("USM")
    )
    curr_data <- rbind(curr_data_univ, curr_data)
  }
}

for (e in 1:E){
  for (r in 1:R){
    curr_data_naive <- data.frame(
      W = W_south,
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
      W = W_south,
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


curr_data_long <- curr_data %>% pivot_longer(cols = c("Mean","CI_upper","CI_lower"), names_to = "line_type", values_to = "estimate")


ggplot(data = curr_data_long, aes(x=W, y=estimate, linetype = line_type, color = Model)) + 
  geom_line(linewidth = .35, aes(color = factor(Model, levels = model_type))) + 
  scale_color_manual(name='Model', labels=c('MSM', 'USM', 'Naive', 'OLS'), values = c("MSM" = "#440154FF", "USM" = "#2A788EFF", "Naive" = "#CD3700", "OLS" = "#66CDAA")) +
  labs(x="Eigenvalue",y="Effect Size") + 
  scale_linetype_manual(values = c("Mean" = "solid", "CI_upper" = "dotted", "CI_lower" = "dotted")) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = .2) + 
  facet_grid(cols=vars(factor(Outcome,level=outcome_names)),rows=vars(factor(Exposure,level=exposure_names)), scales = "free") +
  theme(strip.text = element_text(size=8),
        axis.text = element_text(size=7),
        legend.position = "right") + 
  guides(linetype = "none") # + 
  # scale_colour_viridis_d()

ggsave(filename = paste0("plots/MSM_USM_freq2.png"), width=6, height=4, dpi=300)

ggplot(data = curr_data, aes(x=W, y=Mean, color = Model)) + 
  geom_line(linewidth = .2, aes(color = factor(Model, levels = model_type))) + 
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), linewidth = .5, alpha = .2) + 
  # scale_color_manual(name='Model', labels=c('MSM', 'USM', 'Naive', 'OLS'), values = c("MSM" = "#440154FF", "USM" = "#2A788EFF", "Naive" = "#CD3700", "OLS" = "#66CDAA")) +
  labs(x="Eigenvalue",y="Effect Size") + 
  scale_linetype_manual(values = c("Mean" = "solid", "CI_upper" = "dotted", "CI_lower" = "dotted")) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = .2) + 
  facet_grid(cols=vars(factor(Outcome,level=outcome_names)),rows=vars(factor(Exposure,level=exposure_names)), scales = "free") +
  theme(strip.text = element_text(size=8),
        axis.text = element_text(size=7),
        legend.position = "right") + 
  guides(linetype = "none") + 
  scale_colour_viridis_d()

ggsave(filename = paste0("plots/MSM_USM_freq.png"), width=6, height=4, dpi=300)



########################################################
############# make maps ################################
########################################################

# the full data is not provided 
# the output below can be compiled by data_generation.R

data_south_shape_urban <- readRDS("data_south_shape_urban.rds")

# make maps 

# four themes
themes <- c("theme1", "theme2", "theme3", "theme4")
svi_names <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
for (i in 1:4){
  ggplot(data = data_south_shape_urban) +
    geom_sf(aes(fill = eval(parse(text = paste0("RPL_",themes[i]))), colour = eval(parse(text = paste0("RPL_",themes[i]))))) +
    borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
    theme(legend.position = "none") +
    xlab("") + 
    ylab("") +
    labs(fill = svi_names[i], color = svi_names[i]) +
    theme(legend.position = "right") +
    scale_fill_viridis_c() +
    scale_color_viridis_c()
  ggsave(filename = paste0("plots/svi_", themes[i], ".png"), width=6, height=4)
  
}

# urbanicity

data_south_shape_urban$Urban <- ifelse(data_south_shape_urban$Urbanicity == 0, "Rural", "Urban")

ggplot(data = data_south_shape_urban) +
  geom_sf(aes(fill = factor(Urban), colour = factor(Urban))) +
  borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
  theme(legend.position = "none") +
  xlab("") + 
  ylab("") +
  labs(fill = "Urbanicity", color = "Urbanicity") +
  theme(legend.position = "right") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()
ggsave(filename = paste0("plots/", "urban" , ".png"), width=6, height=4)

# population
ggplot(data = data_south_shape_urban) +
  geom_sf(aes(fill = log(TOTAL), colour = log(TOTAL))) +
  borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
  theme(legend.position = "none") +
  xlab("") + 
  ylab("") +
  labs(fill = "Population", color = "Population") +
  theme(legend.position = "right") +
  scale_fill_viridis_c() +
  scale_color_viridis_c()
ggsave(filename = paste0("plots/", "log_pop" , ".png"), width=6, height=4)


# health conditions
conditions <- c("HYPERT", "CHRNKIDN", "HYPERL", "CHF", "DIABETES")
condition_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")
for (i in 1:5){
  ggplot(data = data_south_shape_urban) +
    geom_sf(aes(fill = log(eval(parse(text = paste0(conditions[i],"_avg_IncRateYear")))), colour = log(eval(parse(text = paste0(conditions[i],"_avg_IncRateYear")))))) +
    borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
    theme(legend.position = "none") +
    xlab("") + 
    ylab("") +
    labs(fill = condition_names[i], color = condition_names[i]) +
    # borders("state", size = .1) + 
    theme(legend.position = "right") +
    scale_fill_viridis_c(limits = c(-6,-0.6)) +
    scale_color_viridis_c(limits = c(-6,-0.6))
  ggsave(filename = paste0("plots/", condition_names[i] , ".png"), width=6, height=4)
}

for (i in 1:5){
  ggplot(data = data_south_shape_urban) +
    geom_sf(aes(fill = log(eval(parse(text = paste0(conditions[i],"_avg_IncRateYear")))), colour = log(eval(parse(text = paste0(conditions[i],"_avg_IncRateYear")))))) +
    borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
    theme(legend.position = "none") +
    xlab("") + 
    ylab("") +
    labs(fill = condition_names[i], color = condition_names[i]) +
    # borders("state", size = .1) + 
    theme(legend.position = "right") +
    scale_fill_viridis_c(limits = c(-6,-0.6)) +
    scale_color_viridis_c(limits = c(-6,-0.6))
  ggsave(filename = paste0("plots/", condition_names[i] , ".png"), width=6, height=4)
}


# health conditions, raw scale (Figures in the supplementary material)
conditions <- c("HYPERT", "CHRNKIDN", "HYPERL", "CHF", "DIABETES")
condition_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")
for (i in 1:5){
  ggplot(data = data_south_shape_urban) +
    geom_sf(aes(fill = (eval(parse(text = paste0(conditions[i],"_avg_IncRateYear")))), colour = (eval(parse(text = paste0(conditions[i],"_avg_IncRateYear")))))) +
    borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
    theme(legend.position = "none") +
    xlab("") + 
    ylab("") +
    labs(fill = condition_names[i], color = condition_names[i]) +
    # borders("state", size = .1) + 
    theme(legend.position = "right") +
    scale_fill_viridis_c() +
    scale_color_viridis_c()
  ggsave(filename = paste0("plots/", condition_names[i] , ".png"), width=6, height=4)
}

# make plots for eigenvectors (Figure 2)
gamma_matrix_south <- cbind(data_south_shape_urban, Gamma_south)

ind <- "11050"
ggplot(data = gamma_matrix_south) +
  geom_sf(aes(fill = eval(parse(text = paste0("X",ind))), colour = eval(parse(text = paste0("X",ind))))) +
  borders("state", regions = c("alabama", "north carolina", "oklahoma", "arkansas", "district of columbia", "delaware", "florida", "georgia", "kentucky", "louisiana", "maryland", "missouri", "mississippi", "south carolina", "tenessee", "texas", "virginia", "west virginia"), size = .1) +
  theme(legend.position = "none") +
  xlab("") + 
  ylab("") +
  # labs(fill = paste0(ind,"th eigenvector"), color = paste0(ind,"th eigenvector")) +
  # borders("state", size = .1) + 
  theme(legend.position = "none") +
  scale_fill_viridis_c() +
  scale_color_viridis_c()
ggsave(filename = "plots/eigenvector_local.png", width=6, height=4)
