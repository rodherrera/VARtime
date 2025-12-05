
require(statmod)
require(esback)

#--------------------------------------------------------------------
## Risk measure estimation 

VaR <-  function(p, int,xi, delta, u,model=NA){
  if(model=="MPOT"){
    out<-u+pmax(delta*(((1-p)/int)^(-1/xi)-1),0)
  }
  return(pmax(out,u))
}

ES<-  function(p, int,xi, delta, u,model=NA){
  
  VaRt<-VaR(p, int,xi, delta, u,model= model)
  
  if(model=="MPOT"){
    out<-VaRt*xi/(xi-1)  + (delta-u)/(xi-1)}
  return(out)
}

#---------------------------------------------------------------------
## Risk measure backtesting

## Accuracy tests for VaR

DQVaR.test<-function(hit,VaR){
  
  defaultW <- getOption("warn") 
  
  options(warn = -1)
  
  if(sum(hit)>0){
    hit<-c(0,hit)
    tt<-length(hit)
    m_tt<-tt-1
    hit1 <- hit[1:m_tt]
    hit2 <- hit[2:tt]
    var2 <- VaR[1:m_tt]
    mylogit<- glm(hit2~hit1+var2, family=binomial(link="logit"), na.action=na.pass)
    logLik(mylogit)
    alpha <- -log(length(hit)/sum(hit)-1)
    loglik1 <- -sum(1-hit2)*alpha-(tt-1)*log(1+exp(-alpha))
    emv <- mylogit$coefficients
    emv1 <- emv[1]
    emv2 <- emv[2]
    emv3 <- emv[3]
    loglik2 <- -sum((1-hit2)*(emv1+emv2*hit1+emv3*var2))-sum(log(1+exp(-emv1-emv2*hit1-emv3*var2)))
    out<- -2*(loglik1-loglik2)
    pvalue<-pchisq(out,2,lower.tail=F)
  }
  else pvalue<-1
  options(warn = defaultW)
  return(pvalue)
}

DQhit.test<-function(hit){
  defaultW <- getOption("warn") 
  
  options(warn = -1)
  
  if(sum(hit)>0){
    hit<-c(0,hit)
    tt<-length(hit)
    m_tt<-tt-1
    hit1 <- hit[1:m_tt]
    hit2 <- hit[2:tt]
    #var2 <- VaR[2:tt]
    mylogit<- glm(hit2~hit1, family=binomial(link="logit"), na.action=na.pass)
    logLik(mylogit)
    alpha <- -log(length(hit)/sum(hit)-1)
    loglik1 <- -sum(1-hit2)*alpha-(tt-1)*log(1+exp(-alpha))
    emv <- mylogit$coefficients
    emv1 <- emv[1]
    emv2 <- emv[2]
    
    loglik2 <- -sum((1-hit2)*(emv1+emv2*hit1))-sum(log(1+exp(-emv1-emv2*hit1)))
    out<- -2*(loglik1-loglik2)
    pvalue<-pchisq(out,1,lower.tail=F)
  }
  
  else pvalue<-1
  options(warn = defaultW)
  return(pvalue)
  
}

Markov.test<-function(hit){
  zz <- 0
  umz <- 0
  zum <- 0
  umum <- 0
  tt<-length(hit)
  m_tt <- tt-1
  for(k in 1:m_tt) {
    i<-k+1
    if (hit[k]==0 & hit[i]==0){
      zz <- zz +1
    }
    else if (hit[k]==0 & hit[i]==1){
      zum <- zum +1
    }
    else if (hit[k]==1 & hit[i]==1){
      umum <- umum +1
    }
    else{
      umz <- umz +1
    }
  }
  p00 <- zz/(zz+zum)
  p01 <- zum/(zz+zum)
  p10 <- umz/(umz+umum)
  p11 <- umum/(umz+umum)
  llp <- (zum+umum)/(zz+umz+zum+umum)
  
  if(sum(hit)<70){
    ll2 <- ((1-llp)^(zz+umz))*(llp^(zum+umum))
    ll1 <- (p00^zz)*(p01^zum)*(p10^umz)*(p11^umum)
    out2<- -2*log(ll2/ll1)}
  
  else{
    ll2 <- (zz+umz)*log(1-llp)+(zum+umum)*log(llp)
    ll1 <- zz*log(p00)+zum*log(p01)+umz*log(p10)+umum*log(p11)
    out2<- -2*(ll2-ll1)}
  
  
  (pchisq(out2,1,lower.tail=F))
}   


MCS.uc <- function(viovec,p,if.cc=F)
  
{
  fcperiod <- length(viovec) #L?nge des Vektors
  nsim <- 10000 #Anzahl Simulationen
  i<-0
  viosumsim <- rep(0,nsim)
  
  if(sum(viovec)<1 && p==0.001){
    pvalue.ut<-1
    pvalue.lt<-1
    pvalue.tt<-1
  }
  
  else{  
    #Berechnungen unter H0
    while (i<nsim)
    {
      i <- i+1
      viovecsim <- sample(x=c(0,1), size = fcperiod, replace = TRUE, prob = c((1-p),p))
      viosumsim[i] <- sum(viovecsim) + rnorm(1,0,0.001)
    }
    
    
    options(digits=10)
    
    #Berechnung Teststatistik
    viosum <- sum(viovec) + rnorm(1,0,0.001)
    
    
    pvalue.ut <- length(which(viosumsim >= viosum))/nsim #upper-tail Test
    pvalue.lt <- length(which(viosumsim <= viosum))/nsim #lower-tail Test
    pvalue.tt <- min(pvalue.ut,pvalue.lt)*2 #two-tailed Test
  }
  if(if.cc==F) return(list("pvalue.ut"=pvalue.ut, "pvalue.lt"=pvalue.lt,"pvalue.tt"=pvalue.tt))
  
  else return(list("viosumsim"=viosumsim,"viosum"=viosum))
}


MCS.iid <- function(viovec,p,if.cc=F)
  
{
  fcperiod <- length(viovec) #n
  
  viopoints <- which(viovec==1) # (iii)
  targetsum <- sum(viovec)      #(ii)
  nsim <- 10000 
  i<-0
  viodurvek <- c() 
  viodur <- (viopoints[1])^2
  
  
  
  
  #Berechnung Teststatistik (v)
  
  if(targetsum>=2){ 
    for (j in c(2:targetsum))
    {
      viodur <- viodur + (viopoints[j]-viopoints[j-1])^2         #Berechnung der Quadratsumme
    }
    
    if (fcperiod>viopoints[targetsum])
    {
      viodur <- viodur + (fcperiod-viopoints[targetsum])^2
    }
  }
  
  viodur <- viodur + rnorm(1,0,0.001)
  
  
  #Berechnungen unter H0 (vi)
  
  
  while (i<nsim)
  {
    i <- i+1
    
    viopointscv <- sort(sample(x=c(0:fcperiod), size = targetsum, replace = FALSE)) #(vii)
    #Generierung von Zufallszahlen
    viodurcv <- (viopointscv[1]-1)^2 #Berechnung der ersten Quadratsumme
    if(targetsum>=2){ 
      for (j in c(2:targetsum))
      {
        viodurcv <- viodurcv + (viopointscv[j]-viopointscv[j-1])^2         #Berechnung der
        
      }
      if (fcperiod>viopointscv[targetsum]) #Berechnung der letzen Quadratsumme
      {
        viodurcv <- viodurcv + (fcperiod-viopointscv[targetsum])^2
      }
    }
    viodurvek[i] <- viodurcv + rnorm(1,0,0.001)                                                                #Abspeichern der Quadratsumme
    
  }
  
  
  
  viodursort <- sort(viodurvek) #(viii +vii+ix)
  
  if(sum(viovec)<1 && p==0.001) pvalue<-1
  else  pvalue <- length(which(viodursort >= viodur))/nsim #Berechnung p-Wert
  
  
  if(if.cc==F) return(list("pvalue"=pvalue))
  else return(list("viodurvek"=viodurvek,"viodur"=viodur))
}


MCS.cc <- function(viovec,p,alpha=1/2)
  
{
  if(sum(viovec)<1 && p==0.001){
    pvalue.ut<-1
    pvalue.lt<-1
    pvalue.tt<-1
  }
  
  else{
    fcperiod<-length(viovec)
    x1<-MCS.uc(viovec,p,if.cc=T)
    x2<-MCS.iid(viovec,p,if.cc=T)
    
    viosumsim<-x1$viosumsim
    viosum<-x1$viosum
    viodurvek<-x2$viodurvek
    viodur<-x2$viodur
    rk<-mean(viosumsim)
    nsim<-length(viodurvek)
    
    y11<-abs((viosum/fcperiod-p)/p)
    y12<-ifelse(viodur>=rk,(viodur-rk)/rk,0)
    test1<-alpha*y11+(1-alpha)*y12
    
    y21<-abs((viosumsim/fcperiod-p)/p)
    y22<-ifelse(viodurvek>=rk,(viodurvek-rk)/rk,0)
    test2<-alpha*y21+(1-alpha)*y22 
    
    
    
    pvalue.ut <- length(which(test2 >= test1))/nsim #upper-tail Test
    pvalue.lt <- length(which(test2 <= test1))/nsim #lower-tail Test
    pvalue.tt <- min(pvalue.ut,pvalue.lt)*2 #two-tailed Test
     #cat(pvalue.ut,"\n")
  }
  return(list("pvalue.ut"=pvalue.ut, "pvalue.lt"=pvalue.lt,"pvalue.tt"=pvalue.tt)) 

}


## Additional tests for VaR

LRuc.test<-function(hit,p){
  n<-length(hit)
  T1<-sum(hit>0)
  T0<-n-T1
  out<--2*log((((1-p)*n/T0)**T0)*((p*n/T1)**T1))
  # (pchisq(out,1,lower.tail=F))
  (pchisq(out,1,lower.tail=F)+pchisq(-out,1,lower.tail=T))
}

LRcc.test<-function(hit,p){
  if(sum(hit)>0){
    n<-length(hit)
    T1<-sum(hit>0)
    T0<-n-T1
    out1<- -2*log((((1-p)*n/T0)**T0)*((p*n/T1)**T1))

    zz <- 0
    umz <- 0
    zum <- 0
    umum <- 0
    tt<-length(hit)
    m_tt <- tt-1
    for(k in 1:m_tt) {
      i<-k+1
      if (hit[k]==0 & hit[i]==0){
        zz <- zz +1
      }
      else if (hit[k]==0 & hit[i]==1){
        zum <- zum +1
      }
      else if (hit[k]==1 & hit[i]==1){
        umum <- umum +1
      }
      else{
        umz <- umz +1
      }
    }
    p00 <- zz/(zz+zum)
    p01 <- zum/(zz+zum)
    p10 <- umz/(umz+umum)
    p11 <- umum/(umz+umum)
    llp <- (zum+umum)/(zz+umz+zum+umum)
    
    if(sum(hit)<100){
      ll2 <- ((1-llp)^(zz+umz))*(llp^(zum+umum))
      ll1 <- (p00^zz)*(p01^zum)*(p10^umz)*(p11^umum)
      out2<- -2*log(ll2/ll1)}
    
    else{
      ll2 <- (zz+umz)*log(1-llp)+(zum+umum)*log(llp)
      ll1 <- zz*log(p00)+zum*log(p01)+umz*log(p10)+umum*log(p11)
      out2<- -2*(ll2-ll1)}
    
    pvalue<- (pchisq(out1+out2,2,lower.tail=F))
    
  }
  else pvalue<-1
  return(pvalue)
}



# Test VaR

backtestVaR<-function(VaR,x){
  
  VaR.05<-VaR[,1]
  VaR.025<-VaR[,2]
  VaR.01<-VaR[,3]
  VaR.001<-VaR[,4]
  n<-length(VaR.05)

  suma.05<-sum(VaR.05[1:(n-1)]<x[2:n])
  suma.025<-sum(VaR.025[1:(n-1)]<x[2:n])
  suma.01<-sum(VaR.01[1:(n-1)]<x[2:n])
  suma.001<-sum(VaR.001[1:(n-1)]<x[2:n])
  
  hit.05=as.numeric(VaR.05[1:(n-1)]<x[2:n])
  hit.025=as.numeric(VaR.025[1:(n-1)]<x[2:n])
  hit.01=as.numeric(VaR.01[1:(n-1)]<x[2:n])
  hit.001=as.numeric(VaR.001[1:(n-1)]<x[2:n])
  
  bin1<-LRuc.test(hit.05,0.05)
  in1<-LRcc.test(hit.05,0.05)
  Ma1<-Markov.test(hit.05)
  Na1<-DQhit.test(hit.05*0.95)
  La1<-DQVaR.test(hit.05*0.95,-VaR.05[1:(n-1)])
  MCSiid1<-MCS.iid(hit.05,0.05)$pvalue  
  MCSuc1<-MCS.uc(hit.05,0.05)
  MCScc1<-MCS.cc(hit.05,0.05)
 
  bin2<-LRuc.test(hit.025,0.025)
  in2<-LRcc.test(hit.025,0.025)
  Ma2<-Markov.test(hit.025)
  Na2<-DQhit.test(hit.025*0.975)
  La2<-DQVaR.test(hit.025*0.975,-VaR.025[1:(n-1)])
  MCSiid2<-MCS.iid(hit.025,0.025)$pvalue  
  MCSuc2<-MCS.uc(hit.025,0.025)
  MCScc2<-MCS.cc(hit.025,0.025)
  
  bin11<-LRuc.test(hit.01,0.01)
  in11<-LRcc.test(hit.01,0.01)
  Ma11<-Markov.test(hit.01)
  Na11<-DQhit.test(hit.01*0.99)
  La11<-DQVaR.test(hit.01*0.99,-VaR.01[1:(n-1)])
  MCSiid11<-MCS.iid(hit.01,0.01)$pvalue  
  MCSuc11<-MCS.uc(hit.01,0.01)
  MCScc11<-MCS.cc(hit.01,0.01)
  
  bin111<-LRuc.test( hit.001,0.001)
  in111<-LRcc.test( hit.001,0.001)
  Ma111<-Markov.test( hit.001)
  Na111<-DQhit.test( hit.001*0.999)
  La111<-DQVaR.test(hit.001*0.999,-VaR.001[1:(n-1)])
  MCSiid111<-MCS.iid(hit.001,0.001)$pvalue  
  MCSuc111<-MCS.uc(hit.001,0.001)
  MCScc111<-MCS.cc(hit.001,0.001)
  
  options(warn=-1)
  
  
 #cmat<- cbind(rbind(suma.05,bin1,Ma1,in1,Na1,La1,MCSuc1$pvalue.tt,MCSiid1,MCScc1$pvalue.tt),
            # rbind(suma.025,bin2,Ma2,in2,Na2,La2,MCSuc2$pvalue.tt,MCSiid2,MCScc2$pvalue.tt),
            # rbind(suma.01,bin11,Ma11,in11,Na11,La11,MCSuc11$pvalue.tt,MCSiid11,MCScc11$pvalue.tt),
             # rbind(suma.001,bin111,Ma111,in111,Na111,La111,MCSuc111$pvalue.tt,MCSiid111,MCScc111$pvalue.tt))
  
cmat<- cbind(rbind(suma.05,bin1,Ma1,in1,Na1,La1,MCSuc1$pvalue.ut,MCSiid1,MCScc1$pvalue.ut),
        rbind(suma.025,bin2,Ma2,in2,Na2,La2,MCSuc2$pvalue.ut,MCSiid2,MCScc2$pvalue.ut),
        rbind(suma.01,bin11,Ma11,in11,Na11,La11,MCSuc11$pvalue.ut,MCSiid11,MCScc11$pvalue.ut),
        rbind(suma.001,bin111,Ma111,in111,Na111,La111,MCSuc111$pvalue.ut,MCSiid111,MCScc111$pvalue.ut))
  
  
  
  colnames(cmat)<-c("0.95","0.975","0.99","0.999")
  cmat<-t(round(cmat,2))
  colnames(cmat)<-c("Exceptions","LRuc","LRind","LRcc","DQhit","DQVaR","MCSuc","MCSiid","MCScc")


  cat("________________________________________________________________________________________________________________________ ","\n")
  printCoefmat(t(cmat), digits=5)
  cat("________________________________________________________________________________________________________________________ ","\n")
  cat("Total Observations ",n,"\n")
  #cmat
}


# Test ES
# ES test are based on the low tail of the return distribution  

backtestES<-function(VaR,ES,x){
  
  n<-length(x)
  # the entries are positive
  r<-  -x

  sESR.05 <-  esr_backtest(r = r, q = -VaR[,1], e = -ES[,1], alpha = 0.05, version = 1)[[1]]
  sESR.025<-  esr_backtest(r = r, q = -VaR[,2], e = -ES[,2], alpha = 0.025, version = 1)[[1]]
  sESR.01 <-  esr_backtest(r = r, q = -VaR[,3], e = -ES[,3], alpha = 0.01, version = 1)[[1]]
  sESR.001<-  tryCatch(esr_backtest(r = r, q = -VaR[,4], e = -ES[,4], alpha = 0.001, version = 1)[[1]],error=function (e) NA)

  aESR.05 <-  esr_backtest(r = r, q = -VaR[,1], e = -ES[,1], alpha = 0.05, version = 2)[[1]]
  aESR.025<-  esr_backtest(r = r, q = -VaR[,2], e = -ES[,2], alpha = 0.025, version = 2)[[1]]
  aESR.01 <-  esr_backtest(r = r, q = -VaR[,3], e = -ES[,3], alpha = 0.01, version = 2)[[1]]
  aESR.001<- tryCatch(esr_backtest(r = r, q = -VaR[,4], e = -ES[,4], alpha = 0.001, version = 2)[[1]],error=function (e) NA)

  siESR.05 <-  esr_backtest(r = r, q = -VaR[,1], e = -ES[,1], alpha = 0.05, version = 3)[[1]]
  siESR.025<-  esr_backtest(r = r, q = -VaR[,2], e = -ES[,2], alpha = 0.025, version = 3)[[1]]
  siESR.01 <-  esr_backtest(r = r, q = -VaR[,3], e = -ES[,3], alpha = 0.01, version = 3)[[1]]
  siESR.001<-  tryCatch(esr_backtest(r = r, q = -VaR[,4], e = -ES[,4], alpha = 0.001, version = 3)[[1]],error=function (e) NA)
 
  ER.05 <-  er_backtest(r = r, q = -VaR[,1], e = -ES[,1])[[1]]
  ER.025<-  er_backtest(r = r, q = -VaR[,2], e = -ES[,2])[[2]]
  ER.01 <-  er_backtest(r = r, q = -VaR[,3], e = -ES[,3])[[1]]
  ER.001<-  er_backtest(r = r, q = -VaR[,4], e = -ES[,4])[[1]]
  
  sCC.05 <-  cc_backtest(r = r, q = -VaR[,1], e = -ES[,1], alpha = 0.05, hommel=T) [[2]]
  sCC.025<-  cc_backtest(r = r, q = -VaR[,2], e = -ES[,2], alpha = 0.025, hommel=T)[[2]]
  sCC.01 <-  cc_backtest(r = r, q = -VaR[,3], e = -ES[,3], alpha = 0.01, hommel=T) [[2]]
  sCC.001<-  cc_backtest(r = r, q = -VaR[,4], e = -ES[,4], alpha = 0.001, hommel=T)[[2]]
  
  
  gCC.05 <-  cc_backtest(r = r, q = -VaR[,1], e = -ES[,1], s=1,alpha = 0.05, hommel=T) [[4]]
  gCC.025<-  cc_backtest(r = r, q = -VaR[,2], e = -ES[,2], s=1,alpha = 0.025, hommel=T)[[4]]
  gCC.01 <-  cc_backtest(r = r, q = -VaR[,3], e = -ES[,3], s=1,alpha = 0.01, hommel=T) [[4]]
  gCC.001<-  cc_backtest(r = r, q = -VaR[,4], e = -ES[,4], s=1,alpha = 0.001, hommel=T)[[4]]
  
  options(warn=-1)
  
  
  cmat<- cbind(rbind(ER.05,sCC.05,gCC.05,sESR.05,aESR.05,siESR.05),
               rbind(ER.025,sCC.025,gCC.025,sESR.025,aESR.025,siESR.025),
               rbind(ER.01,sCC.01,gCC.01,sESR.01,aESR.01,siESR.01),
               rbind(ER.001,sCC.001,gCC.001,sESR.001,aESR.001,siESR.001))
  
  colnames(cmat)<-c("0.95","0.975","0.99","0.999")
  cmat<-t(round(cmat,2))
  colnames(cmat)<-c("ER","sCC","gCC","sERS","aERS","siERS")
  
  
  cat("________________________________________________________________________________________________________________________ ","\n")
  printCoefmat(t(cmat), digits=5)
  cat("________________________________________________________________________________________________________________________ ","\n")
  cat("Total Observations ",n,"\n")
  #cmat
}



## Risk measure backtesting
### Comparative backtesting based on strictly consistent scoring functions (loss functions)

nol2<- function(y, VaR, alfa){
  logy<-log(y)
  logy[is.nan(logy)]<-0
  (1-alfa-(y>VaR))*log(VaR)+(y>VaR)*logy
  }


nol6<-function(y, VaR, ES, alfa){
  (y>VaR)*((y-VaR)/ES)+(1-alfa)*(VaR/ES-1+log(ES))
  }










