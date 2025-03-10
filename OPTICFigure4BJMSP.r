########################################################################################################################################
### Computations for Figure 4, project OPTIC - OPTimal Item Calibration
### 2024-12-31
### Supplementary material for publication: 
### Bjermo J, Fackle-Fornius E, Miller F (2025).
### Optimizing calibration designs with uncertainty in abilities.
### British Journal of Mathematical and Statistical Psychology. 
########################################################################################################################################
  

################# 2-block example
ip <- cbind(c(1.6, 1.6),
              c(-1, 1))
dop  <- 1
eu <- ek <- rep(NA, 20)
for (i in 1:20){
  sop   <- i*6
  yyyu  <- optical(ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)))
  yyyk  <- optical(ip, uncert=FALSE)
  eu[i] <- efficiency(yyyu, ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)), oc="D")
  ek[i] <- efficiency(yyyk, ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)), oc="D")
print(c(sop, ek[[i]]/eu[[i]]))
}
etot2 <- cbind(eu, ek, ek/eu)

################# 3-block example
ip <- cbind(c(1, 2, 2.5),
            c(-1.5, 0.5, 2))
dop  <- 1
eu <- ek <- rep(NA, 20)
for (i in 1:20){
  sop   <- i*6
  yyyu  <- optical(ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)))
  yyyk  <- optical(ip, uncert=FALSE)
  eu[i] <- efficiency(yyyu, ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)), oc="D")
  ek[i] <- efficiency(yyyk, ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)), oc="D")
print(c(sop, ek[i]/eu[i]))
}
etot3 <- cbind(eu, ek, ek/eu)

################# 4-block example
ip <- cbind(c(1.5, 1, 1, 1.5),
              c(-1.5, -0.25, 0.25, 1.5))
dop  <- 1
eu <- ek <- rep(NA, 20)
for (i in 1:20){
  sop   <- i*6
  yyyu  <- optical(ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)))
  yyyk  <- optical(ip, uncert=FALSE)
  eu[i] <- efficiency(yyyu, ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)), oc="D")
  ek[i] <- efficiency(yyyk, ip, uncert=TRUE, ipop=cbind(rep(dop, sop), seq(-1.5, 1.5, length.out=sop)), oc="D")
print(c(sop, ek[i]/eu[i]))
}
etot4 <- cbind(eu, ek, ek/eu)

################# figure
plotarow <- function(etot, ylim1=c(1, 1.45), ylim2=c(0.88, 1), xdata=6*1:20, xlabel=FALSE){
  if (xlabel) xlab0 <- expression(m)  else xlab0 <- ""
  ylim1[2] <- max(c(max(etot[, 1:2]), ylim1[2]))
  ylim2[1] <- min(c(min(etot[, 3]), ylim2[1]))
  plot(c(6, 120), ylim1, type="n", las=1, xlab=xlab0, ylab="Relative efficiency", log="x")
  lines(xdata, etot[, 2], col=4, lwd=2)
  lines(xdata, etot[, 1], col=2, lwd=2)
  plot(c(6, 120), ylim2, type="n", las=1, xlab=xlab0, ylab="", log="x")
  lines(xdata, etot[, 3], col=3, lwd=2)
}
oldpar <- par(mfrow=c(3, 2), oma=c(0.5, 0, 0, 1), mar=c(3.8, 4.5, 0.5, 0))
plotarow(etot2)
plotarow(etot3)
plotarow(etot4, xlabel=T)
par(oldpar)

