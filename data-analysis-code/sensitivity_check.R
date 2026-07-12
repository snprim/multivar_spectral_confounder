rm(list=ls())

# Load the data

theta_star <- readRDS("data/outputs/theta_star.rds")
x_star     <- readRDS("data/outputs/X_star.rds")

# Extract the high-resolution term

use        <- 100
theta_star <- readRDS("data/outputs/theta_star.rds")[1:use,]
x_star     <- readRDS("data/outputs/X_star.rds")[1:use,]

# Compute B matrices

n      <- nrow(theta_star)
E      <- ncol(x_star)
R      <- ncol(theta_star)
Nu     <- 4

Ex     <- eigen(cov(x_star))
B12x   <- Ex$vec[,1:Nu]%*%diag(sqrt(Ex$val[1:Nu]))
Et     <- eigen(cov(theta_star))
B12t   <- Et$vec[,1:Nu]%*%diag(sqrt(Et$val[1:Nu]))

# Compute sensitivities

BB     <- t(B12t%*%B12x)
Vx     <- apply(x_star,2,var)
rho    <- seq(0,0.5,0.01)
rho2   <- rho^2
# Compute sensitivities

BTX    <- B12t%*%t(B12x)
BXX    <- B12x%*%t(B12x)
rho    <- seq(0,0.5,0.01)


# only the four pairs
Gamma_south <- readRDS("data/outputs/Gamma_south.rds")    
X <- readRDS("data/outputs/X.rds")
Z <- readRDS("data/outputs/Z.rds")
Z_standardized <- scale(Z)
X_star <- t(Gamma_south) %*% X
Z_star <- t(Gamma_south) %*% Z_standardized
confound <- cov(X_star[1:100,], Z_star[1:100,])

exposure_names <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
outcome_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

mean <- readRDS("data/outputs/mean.rds")
lower <- readRDS("data/outputs/lower.rds")
upper <- readRDS("data/outputs/upper.rds")

sumz <- round(mean[1,,],3)
upperz <- round(upper[1,,],3)
lowerz <- round(lower[1,,],3)

rownames(sumz) <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
colnames(sumz) <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")


combz <- list(c(1,2), c(1,3), c(1,5), c(4,4))

postscript(file="plots/bias.eps", width=5, height=4, horizontal=FALSE)

op <- par(pty="m", mfrow = c(2,2), mar=c(2,2,2,2), oma = c(3,3,1,1))

combz <- list(c(1,2), c(1,3), c(1,5), c(4,4))

for (comb in combz){
  
  plot(rho,(rho^2)*BTX[comb[2],comb[1]]/BXX[comb[1],comb[1]],type="l", main = paste0(exposure_names[comb[1]], " and ", outcome_names[comb[2]]), cex.main=.8, ylim = c(0,max(abs(upperz[comb[1],comb[2]]), abs(lowerz[comb[1],comb[2]]))))
  abline(h=abs(sumz[comb[1],comb[2]]),lty=1,col="red")
  abline(h=abs(upperz[comb[1],comb[2]]),lty=2,col="red")
  abline(h=abs(lowerz[comb[1],comb[2]]),lty=2,col="red")
  abline(v=min(confound[comb[1],]), lty=3)
  abline(v=max(confound[comb[1],]), lty=3)
}

mtext(expression(rho), side = 1, line = 2, outer = TRUE)
mtext(expression(alpha), side = 2, line = 2, outer = TRUE)

par(op)
dev.off()


# find exact rho value

rho <- .1
for (comb in combz){
  print(abs(lowerz[comb[1], comb[2]]))
  print((rho^2)*BTX[comb[2],comb[1]]/BXX[comb[1],comb[1]])
}

# first pair (1,2)
rho1 <- 0.069
# lower bound
print(abs(lowerz[1,2]))
# bias
print((rho1^2)*BTX[2,1]/BXX[1,1])

# second pair (1,3)
rho2 <- 0.507
# lower bound
print(abs(upperz[1,3]))
# bias
print((rho2^2)*BTX[3,1]/BXX[1,1])


# third pair (1,5)
rho3 <- 0.368
# lower bound
print(abs(lowerz[1,5]))
# bias
print((rho3^2)*BTX[5,1]/BXX[1,1])

# fourth pair (4,4)
rho4 <- 0.117
# lower bound
print(abs(lowerz[4,4]))
# bias
print((rho4^2)*BTX[4,4]/BXX[4,4])
