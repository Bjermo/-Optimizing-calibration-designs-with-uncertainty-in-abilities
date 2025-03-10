# information and inverse information functions
item_prob<-function(a,b,c,x){
  
  
  item_prob   <- c+((1-c)/(1+exp(-a*(x-b))))
  
  
  return(item_prob)
  
}

pq <- function(a,b,x) {
  
  pq <- item_prob(a,b,0,x)*(1-item_prob(a,b,0,x))
  
  return(pq)
}

I <- function(a,b,x) {
  
  I <- (a^2*pq(a,b,x))
  
  return(I)
}

I_inv <- function(a,b,x){
  
  I_inv <- 1/sqrt(sum(I(a,b,x)))
  
  return(I_inv)
}



# define the integrated function
p_tilde <- function(x,a,b,mu,sigma) {(1/(sigma*sqrt(2*pi)))*exp(-0.5*((x-mu)/sigma)^2)*(1/(1+exp(-a*(x-b))))}

p_tilde_ <- function(thetahat,a,b) {
  integrate(p_tilde,a[k],b[k],thetahat,I_inv(a,b,thetahat), lower = -Inf, upper = Inf)$value
}


# probability of a correct response 2PL model
p <- function(x,a,b) {
  
  p <- (1/(1+exp(-a*(x-b))))
  
  return(p)
}

# Calculating the standardised information matrix M
# when no assumption is made about g

M_i <- function(l,k,a,b,thetahats) {
  
  eta_1 <- mean(-(b_new[k] - thetahats[l,])*p(thetahats[l,],a_new[k],b_new[k])*(1 - p(thetahats[l,],a_new[k],b_new[k])))
  
  eta_2 <- mean(-a_new[k]*p(thetahats[l,],a_new[k],b_new[k])*(1 - p(thetahats[l,],a_new[k],b_new[k])))
  
  eta_11 <- eta_1^2
  eta_12 <- eta_21 <- eta_1*eta_2
  eta_22 <- eta_2^2
  
  p_mc <- mean(p(thetahats[l,],a_new[k],b_new[k]))
  M <- (1/(p_mc*(1-p_mc)))*(matrix(c(eta_11,eta_12,eta_21,eta_22),nrow = 2,ncol=2))
  
  return(M)
  
}

# functions for calculating standardized information matrix M
# when g is assumed to be normally distributed

eta_1 <- function(x,a,b,thetahat) {
  ((-(b-x)*exp(-a*(x-b)))/(1 + exp(-a*(x-b)))^2)*dnorm(x,thetahat,I_inv(a_hat,b_hat,thetahat))
}

eta_2 <- function(x,a,b,thetahat) {
  ((-a*exp(-a*(x-b)))/(1 + exp(-a*(x-b)))^2)*dnorm(x,thetahat,I_inv(a_hat,b_hat,thetahat))
}

eta_11<- function(thetahat,a,b){
  
  integrate(eta_1,a[k],b[k],thetahat,lower = -7, upper = 7)$value*integrate(eta_1,a[k],b[k],thetahat,lower = -7, upper = 7)$value
}


eta_12<- function(thetahat,a,b){
  
  integrate(eta_1,a[k],b[k],thetahat,lower = -7, upper = 7)$value*integrate(eta_2,a[k],b[k],thetahat,lower = -7, upper = 7)$value
}


eta_21<- function(thetahat,a,b){
  
  integrate(eta_2,a[k],b[k],thetahat,lower = -7, upper = 7)$value*integrate(eta_1,a[k],b[k],thetahat,lower = -7, upper = 7)$value
}


eta_22<- function(thetahat,a,b){
  
  integrate(eta_2,a[k],b[k],thetahat,lower = -7, upper = 7)$value*integrate(eta_2,a[k],b[k],thetahat,lower = -7, upper = 7)$value
}

# functions used for calculating the standardized information matrix M
# when real theta is used

f11 <- function(x, a, b) { (1/(1+exp(-a*(x-b)))) * (1-(1/(1+exp(-a*(x-b))))) * ((x-b)^2) * dnorm(x) }
f12 <- function(x, a, b) { (1/(1+exp(-a*(x-b)))) * (1-(1/(1+exp(-a*(x-b))))) * (-a*(x-b)) * dnorm(x) }
f22 <- function(x, a, b) { (1/(1+exp(-a*(x-b)))) * (1-(1/(1+exp(-a*(x-b))))) * a^2 * dnorm(x) }
