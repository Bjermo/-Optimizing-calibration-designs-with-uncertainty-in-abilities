################################################################
# This is the code for producing Table 2 in the article
# Optimizing calibration designs with uncertainty in abilities
# published in British Journal of Mathematical and Statistical
# Psychology
###############################################################

# load the mirt package for item parameter and ability estimation
library(mirt)

# save the R-file "Functions for comparing information matrices.R" on
# your computer and use the path in the source command
# Source reads all functions needed for running this code
# Change to the correct path on your computer
source("C:/Users/jonbj99/OneDrive - Linköpings universitet/Optimizing Calibration Designs with Uncertainty in Abilities/Functions/Functions for Comparing information matrices.R")


# Read the saved operational item parameters that was estimated from 
# the second version of the 2018 Swedish SAT. These are needed for 
# genereting the response matrix
# Change to the correct path on your computer
a_org_op <- read.csv("C:/Users/jonbj99/OneDrive - Linköpings universitet/Optimizing Calibration Designs with Uncertainty in Abilities/Results/a_org_op.csv",sep=",",header = T)[,]
b_org_op <- read.csv("C:/Users/jonbj99/OneDrive - Linköpings universitet/Optimizing Calibration Designs with Uncertainty in Abilities/Results/b_org_op.csv",sep=",",header = T)[,]
c_org_op <- read.csv("C:/Users/jonbj99/OneDrive - Linköpings universitet/Optimizing Calibration Designs with Uncertainty in Abilities/Results/c_org_op.csv",sep=",",header = T)[,]


# d is for parametrization in the mirt package
d <- -a_org_op*b_org_op

# sample the abilities
thetas <- as.matrix(rnorm(1000,0,1))


# simulate a response matrix using the estimated item parameters
response <- simdata(a=a_org_op, d, N=1000, itemtype = "2PL",guess=0,Theta= thetas)

#Estimate the model using mirt
mod <- mirt(response,1,"2PL",method="EM")

coef <- coef(mod, simplify=TRUE, IRT=T)$`items`
a_hat = coef[,1]
b_hat = coef[,2]

#Generate the plausible draws
np <- 10000 #nr of plausible draws
thetahats <- fscores(mod,plausible.draws = np,plausible.type="MH")
thetahats_mean <-fscores(mod,full.scores=T,full.scores.SE=T,method="EAP")
thetahats <- matrix(unlist(thetahats),nrow = 1000, ncol = np)

#thetahats is now a matrix with np plausible ability draws for every examine (row)
#thetahats mean is the point estimate and SE of the ability of every examine


# Calculation if the information (M) matrix


k <- 1 #defines the position in the item parameter vectors below
a_new <- c(1,1.5)
b_new <- c(-1,1)


# When no assumption is made about g

z <- sort(thetahats_mean[,1])
x <- rep(0,998)

for (i in 1:(length(z)-1)) {
  
  x[i] = (z[i] + z[i+1])/2
  
}

# S is the information when using plausible draws
S <- matrix(0,ncol = 2, nrow = 2)
for (i in 1:998) {
  S <- S + M_i(i,k,a_new,b_new,thetahats)*dnorm(x[i])*(z[i+1]-z[i])
}


# when g is assumed to be normally distributed

step <- 0.01
x_ <- seq(-7,7,step)
thetahat <-  vector()

for (i in 1:(length(x_)-1)) {
  thetahat[i] = (x_[i+1] + x_[i])/2
}

p_ <- vector()

for(i in 1:length(thetahat)){
  
  p_[i] <- integrate(p_tilde,a_new[k],b_new[k],thetahat[i],I_inv(a_hat,b_hat,thetahat[i]),lower = -7, upper = 7)$value
  
}


eta_11_ <- vector()
eta_12_ <- vector()
eta_21_ <- vector()
eta_22_ <- vector()

for (i in 1:length(thetahat)) {
  
  eta_11_[i] <- eta_11(thetahat[i],a_new,b_new)
  
}

for (i in 1:length(thetahat)) {
  
  eta_12_[i] <- eta_12(thetahat[i],a_new,b_new)
  
}

for (i in 1:length(thetahat)) {
  
  eta_21_[i] <- eta_21(thetahat[i],a_new,b_new)
  
}

for (i in 1:length(thetahat)) {
  
  eta_22_[i] <- eta_22(thetahat[i],a_new,b_new)
  
}

h <- dnorm(thetahat)

M_11 <- sum(((1/(p_*(1-p_)))*eta_11_*h))*step
M_12 <- sum(((1/(p_*(1-p_)))*eta_12_*h))*step
M_21 <- sum(((1/(p_*(1-p_)))*eta_21_*h))*step
M_22 <- sum(((1/(p_*(1-p_)))*eta_22_*h))*step

# M is the information matrix using normal distribution assumption
M <- matrix(c(M_11,M_12,M_21,M_22),nrow = 2, ncol = 2)


# Calculating the information matrix using "true" ability
M_real_11 <- integrate(f11, a_new[k], b_new[k], lower=-Inf, upper=Inf)$value
M_real_12 <- integrate(f12, a_new[k], b_new[k], lower=-Inf, upper=Inf)$value
M_real_22 <- integrate(f22, a_new[k], b_new[k], lower=-Inf, upper=Inf)$value
M_real_21 <- M_real_12

M_real <- matrix(c(M_real_11, M_real_12, M_real_21, M_real_22), nrow = 2, ncol = 2)

# Comparing the results
S
M
M_real

det(S)
det(M)
det(M_real)



