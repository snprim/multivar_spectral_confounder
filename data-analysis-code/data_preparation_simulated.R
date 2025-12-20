
#################################################################
########## If you do have access to Medicare data, ##############
############ use data_preparation.R instead. ####################
############ This file prepares a subsest for NC ################
#################################################################

library(findSVI)
library(sf)
library(tidyverse)
library(spdep)


# set the correct path
setwd("")



# gather SVI data 

svi <- find_svi(year = 2018, geography = "zcta", state = 'US')
svi2 <- svi %>% na.omit()

# SES 

ses <- read_csv("data/outputs/US_ZCTA_ACSvars_2015.csv")
ses2 <- ses[ses$TOTAL>0,] %>% dplyr::select(c("ZCTA", "TOTAL")) %>% na.omit()

# crosswalk between ZCTA, state, region for 48 states 

state_zcta_region2 <- read_csv("data/inputs/zcta_state_region.csv")

# combine svi, health, state_zcta_region2, ses

svi_region <- left_join(svi2, state_zcta_region2, by = c("GEOID" = "ZCTA"))

svi_region_health_ses <- left_join(svi_region, ses2, by = c("GEOID" = "ZCTA")) %>% na.omit()

svi_region_health_ses2 <- svi_region_health_ses %>% dplyr::select(-c("year", "state", "FIPS", "STATE_NAME"))

data_nc <- svi_region_health_ses2[svi_region_health_ses2$STATE_ABBR=="NC",] 

##########################################################
######################### shape file #####################
##########################################################

# the zip file for the shape files can be downloaded at 
# https://www2.census.gov/geo/tiger/TIGER2023/ZCTA520/ 
# This is needed if you want to make maps
# but if maps are not needed, skip to the next section 

shape_zcta <- read_sf("tl_2023_us_zcta520/tl_2023_us_zcta520.shp")
data_nc_shape <- right_join(shape_zcta, data_nc, by = c("ZCTA5CE20" = "GEOID")) %>% na.omit()

# find neighborhood matrix 
# these files are provided in local_files

adj.mat_nc <- nb2mat(poly2nb(data_nc_shape), style = "B", zero.policy = T)
saveRDS(adj.mat_nc, "reproducible_code/data_analysis/adj.mat_nc.rds")
M_nc <- diag(rowSums(adj.mat_nc))
Qmat_nc   <- M_nc - adj.mat_nc
eig_nc <- eigen(Qmat_nc)
Gamma_nc <- eig_nc$vectors
W_nc <- eig_nc$values

####################################################
##### Skip here if no shape files are downkoaded ###
####################################################

data_nc_shape <- readRDS("data/simulated/data_nc_shape.rds")

# urbanicity 

urban <- read_csv("data/inputs/nhgis0015_ds172_2010_zcta.csv")
urban2 <- urban %>% dplyr::select(c("NAME", "H7W001", "H7W002", "H7W003", "H7W004", "H7W005", "H7W006"))
urban3 <- urban2 %>% dplyr::rename(Total=H7W001, 
                                   Urban=H7W002,
                                   UrbanArea=H7W003,
                                   UrbanCluster=H7W004,
                                   Rural=H7W005,
                                   NotDefined=H7W006)
urban3$ZCTA <- str_sub(urban3$NAME,7,11)
urban3$Urbanicity <- ifelse(urban3$Urban/urban3$Total>0.5, 1, 0)
urban4 <- urban3 %>% dplyr::select(c("ZCTA", "Urbanicity"))

# final full data 
data_nc_shape_urban <- left_join(data_nc_shape, urban4, by = c("ZCTA5CE20" = "ZCTA"))

# save NC data as csv
# data_nc <- data_nc_shape_urban %>% dplyr::select(c("ZCTA5CE20", "RPL_theme1", "RPL_theme2", "RPL_theme3", "RPL_theme4", "TOTAL", "Urbanicity"))
# write_csv(data_nc, "reproducible_code/data_analysis/data_nc.csv")

# simulate outcomes from the exposures
# set all coefficients as 0.5
exposure_matrix <- data_nc_shape_urban %>% st_drop_geometry() %>% dplyr::select(c("RPL_theme1", "RPL_theme2", "RPL_theme3", "RPL_theme4"))
sim_beta <- matrix(c(.2, .3, .4, .5, .3, .4, .5, .6, .4, .5, .2, .3, .3, .5, .2, .4, .5, .6, .1, .5), nrow = 4, ncol = 5)

# [,1] [,2] [,3] [,4] [,5]
# [1,]  0.2  0.3  0.4  0.3  0.5
# [2,]  0.3  0.4  0.5  0.5  0.6
# [3,]  0.4  0.5  0.2  0.2  0.1
# [4,]  0.5  0.6  0.3  0.4  0.5

set.seed(919)
X  <- as.matrix(exposure_matrix)
Xb <- X%*%sim_beta
E  <- rnorm(nrow(X)*5, 0, 1)
sim_outcome_matrix  <- 5 + X%*%sim_beta + E

# sim_outcome_matrix <- as.matrix(exposure_matrix) %*% sim_beta +
#   abs(rnorm(nrow(exposure_matrix), 0, 0.1)) 

data_nc_shape_urban$HYPERT_avg_IncRateYear <- sim_outcome_matrix[,1]
data_nc_shape_urban$CHRNKIDN_avg_IncRateYear <- sim_outcome_matrix[,2]
data_nc_shape_urban$HYPERL_avg_IncRateYear <- sim_outcome_matrix[,3]
data_nc_shape_urban$CHF_avg_IncRateYear <- sim_outcome_matrix[,4]
data_nc_shape_urban$DIABETES_avg_IncRateYear <- sim_outcome_matrix[,5]

# create Y, X, Z for analysis

data_nc_shape_urban2 <- st_drop_geometry(data_nc_shape_urban)

exposures <- data_nc_shape_urban2 %>% 
  dplyr::select(c("RPL_theme1", "RPL_theme2", "RPL_theme3", "RPL_theme4")) %>% 
  dplyr::rename(theme1 = RPL_theme1,
                theme2 = RPL_theme2,
                theme3 = RPL_theme3,
                theme4 = RPL_theme4)

outcomes <- data_nc_shape_urban2 %>% dplyr::select(contains("IncRateYear")) %>% log(.)

covariates <- cbind(log(data_nc_shape_urban2$TOTAL), data_nc_shape_urban2$Urbanicity)

Y <- as.matrix(outcomes)
X <- as.matrix(exposures)
Z <- as.matrix(covariates)

Y_standardized <- scale(Y)
Z_standardized <- scale(Z)

Y_star <- t(Gamma_nc) %*% Y_standardized
X_star <- t(Gamma_nc) %*% X
Z_star <- t(Gamma_nc) %*% cbind(1, Z_standardized)


saveRDS(Y_star, "data/simulated/Y_star.rds")
saveRDS(X_star, "data/simulated/X_star.rds")
saveRDS(Z_star, "data/simulated/Z_star.rds")

#################################################################
############## run analysis final_analyses.R ####################
#################################################################


# health conditions
conditions <- c("HYPERT", "CHRNKIDN", "HYPERL", "CHF", "DIABETES")
condition_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")
for (i in 1:5){
  ggplot(data = data_nc_shape_urban) +
    geom_sf(aes(fill = eval(parse(text = paste0(conditions[i],"_avg_IncRateYear"))), colour = eval(parse(text = paste0(conditions[i],"_avg_IncRateYear"))))) +
    borders("state", regions = c("north carolina"), size = .1) +
    theme(legend.position = "none") +
    xlab("") + 
    ylab("") +
    labs(fill = condition_names[i], color = condition_names[i]) +
    # borders("state", size = .1) + 
    theme(legend.position = "right") +
    scale_fill_viridis_c() +
    scale_color_viridis_c()
  ggsave(filename = paste0("data-analysis-code/simulated_plots/", condition_names[i] , ".png"), width=6, height=4)
}

# NC only

Gamma_nc <- cbind(data_nc_shape_urban, Gamma_nc)

ind <- "795"
ggplot(data = Gamma_nc) +
  geom_sf(aes(fill = eval(parse(text = paste0("X",ind))), colour = eval(parse(text = paste0("X",ind))))) +
  borders("state", regions = c("north carolina"), size = .1) +
  theme(legend.position = "none") +
  xlab("") + 
  ylab("") +
  labs(fill = paste0(ind,"th eigenvector"), color = paste0(ind,"th eigenvector")) +
  # borders("state", size = .1) + 
  theme(legend.position = "right") +
  scale_fill_viridis_c() +
  scale_color_viridis_c()

