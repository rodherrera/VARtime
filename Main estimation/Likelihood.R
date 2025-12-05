
## This file contains the code of likelihood function for the parameter estimation and a code of a simple hessian approximation
## for the standard error estimation 


## Tools for the likelihood function 
matrix.imf<-function(x,d){
  mat    <-matrix(x,d,d)
  return(mat)
}

matrix.diag<-function(x,d){
  mat    <-matrix(0,d,d)
  diag(mat)<-x
  return(mat)
}

# spec  Different options for the matrix A (full,diag)
# rv    vector to include realized measures
# rvf   regular varying function (pareto,weibull,burr,mass,shen)

## Likelihood function 
VARtime.nlik<-function(param,zt,rv=NA,jump=NA,umbral,I1,I0,dyn=F,rho=0.9,spec="diag",rvf="pareto", opt=T) #Log-likelihood 
{
  d    <-dim(zt)[2]
  n    <-dim(zt)[1]
  
  
  if(spec=="diag"){
    A  <-matrix.diag(param[1:d],d)
    param.<- param[-c(1:d)]
  }
  
  if(spec=="full"){
    
    A    <-matrix.imf(param[1:(d*d)],d)
    param.<- param[-c(1:(d*d))]
  }
  
  a    <-matrix(param.,length(param.)/d,d)
  
  
  
  dt.  <-a[1,]
  B  <- t(a[2,])
  C  <- t(a[3,])
  S <-matrix(rep(rho,d),1,d)
  
  
  ####
  
  if(dyn & !rvf=="mass"){
    Ak  <- matrix.diag(a[5,],d)
    Bk  <- t(a[4,])
    Ck  <- t(a[6,])
    Dk <- matrix(rep(rho,d),1,d)
    zt.<-EWMA(Dk,zt)
    
    
    k  <- rcppVAR0(Ak,Bk,Ck,zt.)
    k. <-  exp(k)
    if(rvf=="burr" || rvf=="weibull")  tau<-a[7,]
  }
  
  ## non-dynamic
  else{
    if(rvf=="mass")    k.<-1 
    if(rvf=="pareto")  k.<-a[4,] 
    if(rvf=="burr" || rvf=="weibull") {
      k. <-a[4,] 
      tau<-a[5,]
    }
  }
  
  ####
  
  if(rvf=="burr")    Ai   <-k./(k.+umbral^tau)
  if(rvf=="weibull") Ai   <-exp(-k.*umbral^tau) 
  if(rvf=="pareto")  Ai   <-k./(k.+umbral) 
  if(rvf=="mass")    Ai   <-1/(1+umbral) 

  Ai.<-matrix(Ai,n,d,byrow=T)
  
  xt.<- I1-0.1
  xt.<-EWMA(S,xt.)
  
  if(!is.na(sum(rv))){
    m<-dim(a)[1]
    if(!is.na(sum(jump))){
      D   <- t(a[(m-1),])
      E   <- t(a[m,])
      xit  <- rcppVAR3(A,B,C,D,E,log(sqrt(rv)),jump,xt.) 
      
    }
    
    else{
      D   <- t(a[m,])
      xit  <- rcppVAR2(A,B,C,D,log(sqrt(rv)),xt.)
    }
    xit. <-  exp(xit)
  }
  
  else{
    xit  <- rcppVAR1(A,B,C,xt.)
    xit. <-  exp(xit)
  }
  
  phit<-Ai.^(xit.)
  
  #Log-likelihood
  
  logl<-I0*log(1-phit) +I1*log(phit*(xit./dt.)*(1 + zt/dt.)^(-xit.-1)) ## original massaci
  
  if(opt) {
    
    out<- -sum(logl)
  }
  
  else {
    
    
    out<-list( xit=xit., dt=dt., phit=phit,Ai=Ai, k=k.)
    
    
  }
  
  return(out)
}


#simple hessian approximation
hessb     <- function(f, xx, ep = 0.0001, ...){ 
  eps <- ep * xx
  n <- length(xx)
  m <- matrix(0., ncol = n, nrow = n)
  for(i in 1.:n) {
    for(j in 1.:n) {
      x1 <- xx
      x1[i] <- x1[i] + eps[i]
      x1[j] <- x1[j] + eps[j]
      x2 <- xx
      x2[i] <- x2[i] + eps[i]
      x2[j] <- x2[j] - eps[j]
      x3 <- xx
      x3[i] <- x3[i] - eps[i]
      x3[j] <- x3[j] + eps[j]
      x4 <- xx
      x4[i] <- x4[i] - eps[i]
      x4[j] <- x4[j] - eps[j]
      m[i, j] <- (f(x1, ...) - f(x2, ...) - f(x3, ...) + f(
        x4, ...))/(4. * eps[i] * eps[j])
    }
  }
  m
}
