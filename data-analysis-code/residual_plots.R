##########################################################
#### This file makes residual plots in spatial domain ####
##########################################################

library(tidyverse)
library(sp)
library(gstat)

Gamma_south <- readRDS("data/outputs/Gamma_south.rds")

fit_svi_K10_L10_residual <- readRDS("results/fit_svi_K10_L10.rds")
MSM_residual <- rowMeans(fit_svi_K10_L10_residual$residual, dims=2)
spatial_residual <- Gamma_south%*%MSM_residual

data_residual <- cbind(data_south_shape_urban, spatial_residual)

#### Plot residual in maps ####

condition_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

ind <- 1
ggplot(data = data_residual) +
  geom_sf(aes(fill = eval(parse(text = paste0("X",ind))), colour = eval(parse(text = paste0("X",ind))))) +
  # borders("state", regions = c("north carolina", "south carolina"), size = .1) +
  theme(legend.position = "none") +
  xlab("") + 
  ylab("") +
  labs(fill = condition_names[ind], color = condition_names[ind]) +
  # borders("state", size = .1) + 
  theme(legend.position = "right") +
  scale_fill_viridis_c() +
  scale_color_viridis_c()
ggsave(filename = paste0("plots/residuals", condition_names[ind], ".png"), width=6, height=4)

#### density plots ####

# get dataset ready for overlaid density plot

spatial_residual_df <- data.frame(spatial_residual)
colnames(spatial_residual_df) <- condition_names

spatial_residual_long <- pivot_longer(spatial_residual_df, cols = c(Hypertension:Diabetes), names_to = "outcome", values_to = "residual")

spatial_residual_long$outcome <- factor(spatial_residual_long$outcome, levels = condition_names)
ggplot(data = spatial_residual_long, aes(x=residual, fill=outcome)) +
  geom_density(alpha=0.5, position = "identity") + 
  scale_fill_viridis_d() + 
  labs(x="", y="", fill = "Outcome") +
  theme(legend.position = "top", legend.text = element_text(size = 8), legend.title = element_text(size = 8))
ggsave(filename = paste0("plots/hist_residual.png"), width=5, height=4)


#### plot variogram ####


crs = "+proj=aeqd +INTPTLAT20={latitude} +INTPTLON20={longitude} +ellps=WGS84 +units=km"
sf_km <- st_transform(data_residual, crs = crs) 

var1 <- variogram(X1 ~ 1, data = sf_km)
var2 <- variogram(X2 ~ 1, data = sf_km)
var3 <- variogram(X3 ~ 1, data = sf_km)
var4 <- variogram(X4 ~ 1, data = sf_km)
var5 <- variogram(X5 ~ 1, data = sf_km)

postscript(file="plots/vario_residual2.eps", width=5, height=4, horizontal=FALSE)
op <- par(mfrow = c(1,1), mar=c(4.2, 4.2, 1, 6), xpd = T)
plot(var1$dist[1:11], var1$gamma[1:11], ylim = c(0, 0.08), col="black", xlab="Distance (km)", ylab="Semivariance", type = "l", lwd=2)
lines(var2$dist[1:11], var2$gamma[1:11], lty=2, col="brown", lwd=2)
lines(var3$dist[1:11], var3$gamma[1:11], lty=6, col="magenta", lwd=2)
lines(var4$dist[1:11], var4$gamma[1:11], lty=4, col="blue", lwd=2)
lines(var5$dist[1:11], var5$gamma[1:11], lty=5, col="darkgreen", lwd=2)
legend(1075, 0.08, legend=c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes"), col=c("black", "brown", "magenta", "blue", "darkgreen"), lty=c(1,2,6,4,5), horiz =F, cex=.55, lwd = c(2,2,2,2,2))
par(op)
dev.off()



#############################
# Residual vs Z (in the spatial domain)

Z <- readRDS("data/outputs/Z.rds")
Gamma_south <- readRDS("data/outputs/Gamma_south.rds")
fit_svi_K10_L10_residual <- readRDS("results/fit_svi_K10_L10.rds")
MSM_residual <- rowMeans(fit_svi_K10_L10_residual$residual, dims=2)
spatial_residual <- Gamma_south%*%MSM_residual

dataset <- cbind(Z, spatial_residual)
colnames(dataset) <- c("Population", "Urbanicity", "Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")
dataset <- as.data.frame(dataset)

condition_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

plotz <- list()
for (r in 1:5){
  plotz[[r]] <- ggplot(data = dataset, aes(x = Population, y = eval(parse(text=condition_names[r])))) + 
    geom_point() + 
    geom_smooth(method = "loess", se = FALSE, linewidth = 2) +
    labs(title=paste0(condition_names[r])) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("") + 
    xlab("")
}

library(ggpubr)


pdf(paste0("plots/residual_pop.pdf"), width = 9, height = 6)
result_all <- ggarrange(plotz[[1]], plotz[[2]], plotz[[3]], plotz[[4]], plotz[[5]], nrow = 2, ncol = 3, common.legend = TRUE, hjust = -1)
annotate_figure(result_all, left = text_grob("Residual", rot = 90, size = 16), bottom = text_grob("Log Population", size = 16))
dev.off()
  
