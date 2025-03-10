#######################################################################################################
### Simulations for project OPTIC - OPImal Item Calibration
### 2024-12-31
### Supplementary material for publication: 
### Bjermo J, Fackle-Fornius E, Miller F (2025).
### Optimizing calibration designs with uncertainty in abilities.
### British Journal of Mathematical and Statistical Psychology. 
#######################################################################################################

library(mirt)
library(optical)

calitems <- function(yyy, abil){
  ci   <- NULL
  for (i in 1:length(yyy$ht)){
    for (j in 1:dim(yyy$ht[[i]])[1]){
      if (yyy$ht[[i]][j, 1] <= abil && abil < yyy$ht[[i]][j, 2]) 
        ci <- c(ci, yyy$ht[[i]][j, 3])
    }
  }
  ci
}

N   <-  500  # Number of examinees in simulated tests
sim <- 1000  # Number of repetitions for simulation study

##############################################################################################################
### choose one of the following three item parameters (for 2-, 3-, or 4-item-block, respectively)
##############################################################################################################
ip <- cbind(c(1.6, 1.6),
              c(-1, 1))
ip <- cbind(c(1, 2, 2.5),
            c(-1.5, 0.5, 2))
ip <- cbind(c(1.5, 1, 1, 1.5),
              c(-1.5, -0.25, 0.25, 1.5))

##################
### simulations 
##################
sop  <- 30    # size of operational test
dop  <- 1     # discriminaton of items in operational test
scal <- dim(ip)[1]
ipop <- cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop))

yyyu <- optical(ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)))
yyyk <- optical(ip, uncert=FALSE)

set.seed(2025)
allresau <- allresbu <- allresak <- allresbk <- allresar <- allresbr <- NULL
for (s in 1:sim){

  # Generate true abilities for N examinees and estimate them based on operational test part with mirt
  theta    <- rnorm(N)
  accuracy <- matrix(rep(NA, N*sop), nrow=N, dimnames=list(NULL, 1:sop))
  for (i in 1:sop){
    accuracy[, i] <- rbinom(N, size=1, prob=1/(1+exp(-ipop[i, 1]*(theta-ipop[i, 2]))))
  }
  mod <- mirt(accuracy, model=1, itemtype="2PL")
  thetaest  <- fscores(mod)
  dmat      <- oumat <- okmat <- rdmat <- matrix(rep(NA, N*scal), ncol=scal)

  # Generate first complete datamatrix and apply then ODu, ODk, RD
  for (i in 1:scal){
    dmat[, i] <- rbinom(N, size=1, prob=1/(1+exp(-ip[i, 1]*(theta-ip[i, 2]))))
  }
  yyyr <- sample(scal, size=N, replace=TRUE)
  for (i in 1:scal)
    for (j in 1:N){
      if (calitems(yyyu, thetaest[j])==i)  oumat[j, i] <- dmat[j, i]
      if (calitems(yyyk, thetaest[j])==i)  okmat[j, i] <- dmat[j, i]
      if (yyyr[j]==i)  rdmat[j, i] <- dmat[j, i]
  }

  # estimate item parameters, based on estimated abilities
  resau <- resbu <- resak <- resbk <- resar <- resbr <- NULL
  for (i in 1:scal){
    betahatu <- glm(oumat[, i]~thetaest, family=binomial)$coef
    betahatk <- glm(okmat[, i]~thetaest, family=binomial)$coef
    betahatr <- glm(rdmat[, i]~thetaest, family=binomial)$coef
    names(betahatu) <- NULL
    names(betahatk) <- NULL
    names(betahatr) <- NULL
    ahatu    <- betahatu[2]
    bhatu    <- -betahatu[1]/betahatu[2]
    ahatk    <- betahatk[2]
    bhatk    <- -betahatk[1]/betahatk[2]
    ahatr    <- betahatr[2]
    bhatr    <- -betahatr[1]/betahatr[2]
    resau    <- c(resau, ahatu)
    resbu    <- c(resbu, bhatu)
    resak    <- c(resak, ahatk)
    resbk    <- c(resbk, bhatk)
    resar    <- c(resar, ahatr)
    resbr    <- c(resbr, bhatr)
  }
  allresau <- rbind(allresau, resau)
  allresbu <- rbind(allresbu, resbu)
  allresak <- rbind(allresak, resak)
  allresbk <- rbind(allresbk, resbk)
  allresar <- rbind(allresar, resar)
  allresbr <- rbind(allresbr, resbr)
}

biasau <- biasbu <- biasak <- biasbk <- biasar <- biasbr <- NULL
for (i0 in 1:scal){
  biasau <- c(biasau, mean(allresau[, i0] - ip[i0, 1]))
  biasbu <- c(biasbu, mean(allresbu[, i0] - ip[i0, 2]))
  biasak <- c(biasak, mean(allresak[, i0] - ip[i0, 1]))
  biasbk <- c(biasbk, mean(allresbk[, i0] - ip[i0, 2]))
  biasar <- c(biasar, mean(allresar[, i0] - ip[i0, 1]))
  biasbr <- c(biasbr, mean(allresbr[, i0] - ip[i0, 2]))
}
mseau <- msebu <- mseak <- msebk <- msear <- msebr <- NULL
for (i0 in 1:scal){
  mseau <- c(mseau, mean((allresau[, i0] - ip[i0, 1])^2))
  msebu <- c(msebu, mean((allresbu[, i0] - ip[i0, 2])^2))
  mseak <- c(mseak, mean((allresak[, i0] - ip[i0, 1])^2))
  msebk <- c(msebk, mean((allresbk[, i0] - ip[i0, 2])^2))
  msear <- c(msear, mean((allresar[, i0] - ip[i0, 1])^2))
  msebr <- c(msebr, mean((allresbr[, i0] - ip[i0, 2])^2))
}

# Results in Table for bias for the chosen item-block:
round(biasau, 2)
round(biasak, 2)
round(biasar, 2)

round(biasbu, 2)
round(biasbk, 2)
round(biasbr, 2)

# Results in Table for MSE for the chosen item-block:
round(mseau, 3)
round(mseak, 3)
round(msear, 3)

round(msebu, 3)
round(msebk, 3)
round(msebr, 3)



