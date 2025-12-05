
setwd("C:Folder/Main estimation") # Choose Directory

source("Likelihood.R")
source("Accuracy Tests.R")
source("filter.arma.R")
library(zoo)
library(HelpersMG)
library(e1071)
library(evd)
library(MCS)
library("ggplot2")


###################################################################################################################################
# Data
###################################################################################################################################

# Getting and ordering the data
# The stocks are AAPL, AMZN, CSCO, IBM, MSFT

ret<-read.zoo("rets_tech_stocks.csv",header=TRUE, format="%d/%m/%Y",sep=",")
rv<-read.zoo("rv_tech_stocks.csv", header=TRUE,  format="%d/%m/%Y",sep=",")
rsv_n<-read.zoo("rsv_n_tech_stocks.csv", header=TRUE,  format="%d/%m/%Y",sep=",")
rsv_p<-read.zoo("rsv_p_tech_stocks.csv", header=TRUE,  format="%d/%m/%Y",sep=",")
bpv<-read.zoo("bpv_tech_stocks.csv",header=TRUE,  format="%d/%m/%Y",sep=",")
jumps<-read.zoo("jumps_tech_stocks.csv",header=TRUE,  format="%d/%m/%Y",sep=",")


# Definition of in- and out-of-sample data
#  in 01.08.2002 - 31.12.2018
#  out 02.01.2019 - 29.12.2023

## negative returns
yt.in<-  -ret[8:4084,]
yt.out<- -ret[4085:5333,]
head(yt.in)
tail(yt.out)

## realized measures
## standard realized volatility
rv.in<- rv[8:4084,]
rv.out<-rv[4085:5333,]
## negative realized semivariance
rsv_n.in<-rsv_n[8:4084,]
rsv_n.out<-rsv_n[4085:5333,]
## positive realized semivariance
rsv_p.in<-rsv_p[8:4084,]
rsv_p.out<-rsv_p[4085:5333,]
## bipower variation and jumps
bpv.in<-bpv[8:4084,]
bpv.out<-bpv[4085:5333,]
jumps.in<-jumps[8:4084,]
jumps.out<-jumps[4085:5333,]


## in-sample analysis
dim(yt.in)
d<-dim(yt.in)[2]
d
p<-0.9
umbral.in   <-apply(yt.in,2,quantile,p)
zt.in <-pmax(t(t(yt.in)-umbral.in),0)
I0.in<- (zt.in==0)*1
I1.in<- (zt.in>0)*1
n.in<-dim(zt.in)[1]


## out-of-sample analysis
dim(yt.out)
d<-dim(yt.out)[2]
d
p<-0.9
umbral.out   <-apply(yt.out,2,quantile,p)
zt.out <-pmax(t(t(yt.out)-umbral.out),0)
I0.out<- (zt.out==0)*1
I1.out<- (zt.out>0)*1
n.out<-dim(zt.out)[1]


########################################################################################################################################
# Estimation 
########################################################################################################################################

# p stands for Pareto
# b stands for Burr
# w stands for Weibull 

## estimation of diagonal,non-dynamic, without realized measures
fit.diag.p<- nlminb(runif(4*d + d ,0,0.1), rv=NA,jump=NA,rvf="pareto",spec="diag",dyn=F, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.b<- nlminb(runif(5*d + d ,0,0.1), rv=NA,jump=NA,rvf="burr",spec="diag", dyn=F,VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.w<- nlminb(runif(5*d + d ,0,0.1), rv=NA,jump=NA,rvf="weibull",spec="diag", dyn=F,VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

## maximization of the likelihood
fit.diag.p$objective
fit.diag.b$objective
fit.diag.w$objective

## hessian and standard errors
hess.diag.p<-hessb(VARtime.nlik,fit.diag.p$par, rv=NA,jump=NA,rvf="pareto", spec="diag", dyn=F,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.p<-SEfromHessian(hess.diag.p)

hess.diag.b<-hessb(VARtime.nlik,fit.diag.b$par, rv=NA,jump=NA,rvf="burr", spec="diag", dyn=F,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-06)
se.diag.b<-SEfromHessian(hess.diag.b)

hess.diag.w<-hessb(VARtime.nlik,fit.diag.w$par, rv=NA,jump=NA,rvf="weibull", spec="diag", dyn=F,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-06)
se.diag.w<-SEfromHessian(hess.diag.w)

## out-of-sample estimation
out.diag.p<-VARtime.nlik(fit.diag.p$par,rv=NA,jump=NA,rvf="pareto", spec="diag", dyn=F,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.b<-VARtime.nlik(fit.diag.b$par,rv=NA,jump=NA,rvf="burr", spec="diag", dyn=F,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.w<-VARtime.nlik(fit.diag.w$par,rv=NA,jump=NA,rvf="weibull", spec="diag", dyn=F,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


#-------------------------------------------------------------------------------
## estimation of diagonal, dynamic, without realized measures
fit.diag.pd<- nlminb(runif(6*d + d,0,0.1), rv=NA,jump=NA,rvf="pareto",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.bd<- nlminb(runif(7*d + d ,0,0.1), rv=NA,jump=NA,rvf="burr",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.wd<- nlminb(runif(7*d + d ,0,0.1), rv=NA,jump=NA,rvf="weibull",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.diag.pd$objective
fit.diag.bd$objective
fit.diag.wd$objective

hess.diag.pd<-hessb(VARtime.nlik,fit.diag.pd$par, rv=NA,jump=NA,rvf="pareto", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-06)
se.diag.pd<-SEfromHessian(hess.diag.pd)

hess.diag.bd<-hessb(VARtime.nlik,fit.diag.bd$par, rv=NA,jump=NA,rvf="burr", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.bd<-SEfromHessian(hess.diag.bd)

hess.diag.wd<-hessb(VARtime.nlik,fit.diag.wd$par, rv=NA,jump=NA,rvf="weibull", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-010)
se.diag.wd<-SEfromHessian(hess.diag.wd)

out.diag.pd<-VARtime.nlik(fit.diag.pd$par,rv=NA,jump=NA,rvf="pareto", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.bd<-VARtime.nlik(fit.diag.bd$par,rv=NA,jump=NA,rvf="burr", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.wd<-VARtime.nlik(fit.diag.wd$par,rv=NA,jump=NA,rvf="weibull", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


#-------------------------------------------------------------------------------
## estimation of full (non-diagonal), non-dynamic, without realized measures
fit.full.p<- nlminb(runif(4*d + d*d,0,0.1), rv=NA,jump=NA,rvf="pareto",spec="full", dyn=F, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.b<- nlminb(runif(5*d + d*d,0,0.1), rv=NA,jump=NA,rvf="burr",spec="full", dyn=F, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.w<- nlminb(runif(5*d + d*d,0,0.1), rv=NA,jump=NA,rvf="weibull",spec="full", dyn=F, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.full.p$objective
fit.full.b$objective
fit.full.w$objective

hess.full.p<-hessb(VARtime.nlik,fit.full.p$par, rv=NA,jump=NA,rvf="pareto", spec="full", dyn=F,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.p<-SEfromHessian(hess.full.p)

hess.full.b<-hessb(VARtime.nlik,fit.full.b$par, rv=NA,jump=NA,rvf="burr", spec="full", dyn=F,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.b<-SEfromHessian(hess.full.b)

hess.full.w<-hessb(VARtime.nlik,fit.full.w$par, rv=NA,jump=NA,rvf="weibull", spec="full", dyn=F,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.w<-SEfromHessian(hess.full.w)

out.full.p<-VARtime.nlik(fit.full.p$par,rv=NA,jump=NA,rvf="pareto", spec="full", dyn=F,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.b<-VARtime.nlik(fit.full.b$par,rv=NA,jump=NA,rvf="burr", spec="full", dyn=F,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.w<-VARtime.nlik(fit.full.w$par,rv=NA,jump=NA,rvf="weibull", spec="full", dyn=F,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


#-------------------------------------------------------------------------------
## estimation of full (non-diagonal), dynamic, without realized measures

fit.full.pd<- nlminb(runif(6*d + d*d,0,0.1), rv=NA,jump=NA,rvf="pareto",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.bd<- nlminb(runif(7*d + d*d,0,0.1), rv=NA,jump=NA,rvf="burr",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.wd<- nlminb(runif(7*d + d*d,0,0.1), rv=NA,jump=NA,rvf="weibull",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.full.pd$objective
fit.full.bd$objective
fit.full.wd$objective

hess.full.pd<-hessb(VARtime.nlik,fit.full.pd$par, rv=NA,jump=NA,rvf="pareto", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.pd<-SEfromHessian(hess.full.pd)

hess.full.bd<-hessb(VARtime.nlik,fit.full.bd$par, rv=NA,jump=NA,rvf="burr", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.bd<-SEfromHessian(hess.full.bd)

hess.full.wd<-hessb(VARtime.nlik,fit.full.wd$par, rv=NA,jump=NA,rvf="weibull", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.wd<-SEfromHessian(hess.full.wd)

out.full.pd<-VARtime.nlik(fit.full.pd$par,rv=NA,jump=NA,rvf="pareto", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.bd<-VARtime.nlik(fit.full.bd$par,rv=NA,jump=NA,rvf="burr", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.wd<-VARtime.nlik(fit.full.wd$par,rv=NA,jump=NA,rvf="weibull", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


#-------------------------------------------------------------------------------
## estimation of diagonal, dynamic, with realized measures

fit.diag.pdr<- nlminb(runif(7*d + d,0,0.1), rv=rv.in,jump=NA,rvf="pareto",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.bdr<- nlminb(runif(8*d + d ,0,0.1), rv=rv.in,jump=NA,rvf="burr",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.wdr<- nlminb(runif(8*d + d ,0,0.1), rv=rv.in,jump=NA,rvf="weibull",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.diag.pdr$objective
fit.diag.bdr$objective
fit.diag.wdr$objective

hess.diag.pdr<-hessb(VARtime.nlik,fit.diag.pdr$par, rv=rv.in,jump=NA,rvf="pareto", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.pdr<-SEfromHessian(hess.diag.pdr)

hess.diag.bdr<-hessb(VARtime.nlik,fit.diag.bdr$par, rv=rv.in,jump=NA,rvf="burr", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-06)
se.diag.bdr<-SEfromHessian(hess.diag.bdr)

hess.diag.wdr<-hessb(VARtime.nlik,fit.diag.wdr$par, rv=rv.in,jump=NA,rvf="weibull", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.wdr<-SEfromHessian(hess.diag.wdr)

out.diag.pdr<-VARtime.nlik(fit.diag.pdr$par,rv=rv.out,jump=NA,rvf="pareto", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.bdr<-VARtime.nlik(fit.diag.bdr$par,rv=rv.out,jump=NA,rvf="burr", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.wdr<-VARtime.nlik(fit.diag.wdr$par,rv=rv.out,jump=NA,rvf="weibull", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


#-------------------------------------------------------------------------------
## estimation of diagonal, dynamic, with realized measures (negative RSV)

fit.diag.pdrs<- nlminb(runif(7*d + d,0,0.1), rv=rsv_n.in,jump=NA,rvf="pareto",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.bdrs<- nlminb(runif(8*d + d ,0,0.1), rv=rsv_n.in,jump=NA,rvf="burr",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.wdrs<- nlminb(runif(8*d + d ,0,0.1), rv=rsv_n.in,jump=NA,rvf="weibull",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.diag.pdrs$objective
fit.diag.bdrs$objective
fit.diag.wdrs$objective

hess.diag.pdrs<-hessb(VARtime.nlik,fit.diag.pdrs$par, rv=rsv_n.in,jump=NA,rvf="pareto", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-08)
se.diag.pdrs<-SEfromHessian(hess.diag.pdrs)

hess.diag.bdrs<-hessb(VARtime.nlik,fit.diag.bdrs$par, rv=rsv_n.in,jump=NA,rvf="burr", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-08)
se.diag.bdrs<-SEfromHessian(hess.diag.bdrs)

hess.diag.wdrs<-hessb(VARtime.nlik,fit.diag.wdrs$par, rv=rsv_n.in,jump=NA,rvf="weibull", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.wdrs<-SEfromHessian(hess.diag.wdrs)

out.diag.pdrs<-VARtime.nlik(fit.diag.pdrs$par,rv=rsv_n.out,jump=NA,rvf="pareto", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.bdrs<-VARtime.nlik(fit.diag.bdrs$par,rv=rsv_n.out,jump=NA,rvf="burr", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.wdrs<-VARtime.nlik(fit.diag.wdrs$par,rv=rsv_n.out,jump=NA,rvf="weibull", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)

#-------------------------------------------------------------------------------
## estimation of diagonal, dynamic, with realized measures (negative and positive RSV)

fit.diag.pdpn<- nlminb(runif(8*d + d,0,0.1), rv=rsv_n.in,jump=rsv_p.in,rvf="pareto",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.bdpn<- nlminb(runif(9*d + d,0,0.1), rv=rsv_n.in,jump=rsv_p.in,rvf="burr",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.wdpn<- nlminb(runif(9*d + d,0,0.1), rv=rsv_n.in,jump=rsv_p.in,rvf="weibull",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.diag.pdpn$objective
fit.diag.bdpn$objective
fit.diag.wdpn$objective

hess.diag.pdpn<-hessb(VARtime.nlik,fit.diag.pdpn$par, rv=rsv_n.in,jump=rsv_p.in,rvf="pareto", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-07)
se.diag.pdpn<-SEfromHessian(hess.diag.pdpn)

hess.diag.bdpn<-hessb(VARtime.nlik,fit.diag.bdpn$par, rv=rsv_n.in,jump=rsv_p.in,rvf="burr", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.bdpn<-SEfromHessian(hess.diag.bdpn)

hess.diag.wdpn<-hessb(VARtime.nlik,fit.diag.wdpn$par, rv=rsv_n.in,jump=rsv_p.in,rvf="weibull", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.wdpn<-SEfromHessian(hess.diag.wdpn)

out.diag.pdpn<-VARtime.nlik(fit.diag.pdpn$par,rv=rsv_n.out,jump=rsv_p.out,rvf="pareto", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.bdpn<-VARtime.nlik(fit.diag.bdpn$par,rv=rsv_n.out,jump=rsv_p.out,rvf="burr", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.wdpn<-VARtime.nlik(fit.diag.wdpn$par,rv=rsv_n.out,jump=rsv_p.out,rvf="weibull", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


#-------------------------------------------------------------------------------
## estimation of diagonal, dynamic, with realized measures (Bipower and jumps)

fit.diag.pdbj<- nlminb(runif(8*d + d,0,0.1), rv=bpv.in,jump=jumps.in,rvf="pareto",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.bdbj<- nlminb(runif(9*d + d ,0,0.1), rv=bpv.in,jump=jumps.in,rvf="burr",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.diag.wdbj<- nlminb(runif(9*d + d ,0,0.1), rv=bpv.in,jump=jumps.in,rvf="weibull",spec="diag", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.diag.pdbj$objective
fit.diag.bdbj$objective
fit.diag.wdbj$objective

hess.diag.pdbj<-hessb(VARtime.nlik,fit.diag.pdbj$par,  rv=bpv.in,jump=jumps.in,rvf="pareto", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-06)
se.diag.pdbj<-SEfromHessian(hess.diag.pdbj)

hess.diag.bdbj<-hessb(VARtime.nlik,fit.diag.bdbj$par,  rv=bpv.in,jump=jumps.in,rvf="burr", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.bdbj<-SEfromHessian(hess.diag.bdbj)

hess.diag.wdbj<-hessb(VARtime.nlik,fit.diag.wdbj$par,  rv=bpv.in,jump=jumps.in,rvf="weibull", spec="diag", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.diag.wdbj<-SEfromHessian(hess.diag.wdbj)

out.diag.pdbj<-VARtime.nlik(fit.diag.pdbj$par,rv=bpv.out,jump=jumps.out,rvf="pareto", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.bdbj<-VARtime.nlik(fit.diag.bdbj$par,rv=bpv.out,jump=jumps.out,rvf="burr", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.diag.wdbj<-VARtime.nlik(fit.diag.wdbj$par,rv=bpv.out,jump=jumps.out,rvf="weibull", spec="diag", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)

#-------------------------------------------------------------------------------
## estimation of full (non-diagonal), dynamic, with realized measures (standard measure)
fit.full.pdr<- nlminb(runif(7*d + d*d,0,0.1), rv=rv.in,jump=NA, rvf="pareto",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.bdr<- nlminb(runif(8*d + d*d,0,0.1), rv=rv.in,jump=NA,rvf="burr",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.wdr<- nlminb(runif(8*d + d*d,0,0.1), rv=rv.in,jump=NA,rvf="weibull",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.full.pdr$objective
fit.full.bdr$objective
fit.full.wdr$objective

fit.full.pdr$par
fit.full.bdr$par
fit.full.wdr$par

hess.full.pdr<-hessb(VARtime.nlik,fit.full.pdr$par,  rv=rv.in,jump=NA,rvf="pareto", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.pdr<-SEfromHessian(hess.full.pdr)

hess.full.bdr<-hessb(VARtime.nlik,fit.full.bdr$par,  rv=rv.in,jump=NA,rvf="burr", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.bdr<-SEfromHessian(hess.full.bdr)

hess.full.wdr<-hessb(VARtime.nlik,fit.full.wdr$par,  rv=rv.in,jump=NA,rvf="weibull", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.wdr<-SEfromHessian(hess.full.wdr)

out.full.pdr<-VARtime.nlik(fit.full.pdr$par,rv=rv.out,jump=NA,rvf="pareto", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.bdr<-VARtime.nlik(fit.full.bdr$par,rv=rv.out,jump=NA,rvf="burr", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.wdr<-VARtime.nlik(fit.full.wdr$par,rv=rv.out,jump=NA,rvf="weibull", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)

#-------------------------------------------------------------------------------
## estimation of full (non-diagonal), dynamic, with realized measures (negative RSV)

fit.full.pdrs<- nlminb(runif(7*d + d*d,0,0.1), rv=rsv_n.in,jump=NA,rvf="pareto",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.bdrs<- nlminb(runif(8*d + d*d,0,0.1), rv=rsv_n.in,jump=NA,rvf="burr",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.wdrs<- nlminb(runif(8*d + d*d,0,0.1), rv=rsv_n.in,jump=NA,rvf="weibull",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.full.pdrs$objective
fit.full.bdrs$objective
fit.full.wdrs$objective

hess.full.pdrs<-hessb(VARtime.nlik,fit.full.pdrs$par,  rv=rsv_n.in,jump=NA,rvf="pareto", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-08)
se.full.pdrs<-SEfromHessian(hess.full.pdrs)

hess.full.bdrs<-hessb(VARtime.nlik,fit.full.bdrs$par,  rv=rsv_n.in,jump=NA,rvf="burr", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.bdrs<-SEfromHessian(hess.full.bdrs)

hess.full.wdrs<-hessb(VARtime.nlik,fit.full.wdrs$par,  rv=rsv_n.in,jump=NA,rvf="weibull", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.wdrs<-SEfromHessian(hess.full.wdrs)

out.full.pdrs<-VARtime.nlik(fit.full.pdrs$par,rv=rsv_n.out,jump=NA,rvf="pareto", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.bdrs<-VARtime.nlik(fit.full.bdrs$par,rv=rsv_n.out,jump=NA,rvf="burr", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.wdrs<-VARtime.nlik(fit.full.wdrs$par,rv=rsv_n.out,jump=NA,rvf="weibull", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)

#-------------------------------------------------------------------------------
## estimation of full (non-diagonal), dynamic, with realized measures (positive RSV )

fit.full.pdpn<- nlminb(runif(8*d + d*d,0,0.1), rv=rsv_n.in,jump=rsv_p.in,rvf="pareto",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.bdpn<- nlminb(runif(9*d + d*d,0,0.1), rv=rsv_n.in,jump=rsv_p.in,rvf="burr",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.wdpn<- nlminb(runif(9*d + d*d,0,0.1), rv=rsv_n.in,jump=rsv_p.in,rvf="weibull",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.full.pdpn$objective
fit.full.bdpn$objective
fit.full.wdpn$objective

hess.full.pdpn<-hessb(VARtime.nlik,fit.full.pdpn$par,  rv=rsv_n.in,jump=rsv_p.in,rvf="pareto", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.pdpn<-SEfromHessian(hess.full.pdpn)

hess.full.bdpn<-hessb(VARtime.nlik,fit.full.bdpn$par,  rv=rsv_n.in,jump=rsv_p.in,rvf="burr", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.bdpn<-SEfromHessian(hess.full.bdpn)

hess.full.wdpn<-hessb(VARtime.nlik,fit.full.wdpn$par,  rv=rsv_n.in,jump=rsv_p.in,rvf="weibull", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.wdpn<-SEfromHessian(hess.full.wdpn)

out.full.pdpn<-VARtime.nlik(fit.full.pdpn$par,rv=rsv_n.out,jump=rsv_p.out,rvf="pareto", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.bdpn<-VARtime.nlik(fit.full.bdpn$par,rv=rsv_n.out,jump=rsv_p.out,rvf="burr", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1.out,I0=I0.out,opt=F)
out.full.wdpn<-VARtime.nlik(fit.full.wdpn$par,rv=rsv_n.out,jump=rsv_p.out,rvf="weibull", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)

#-------------------------------------------------------------------------------
## estimation of full (non-diagonal), dynamic, with realized measures (Bipower and jumps)
fit.full.pdbj<- nlminb(runif(8*d + d*d,0,0.1), rv=bpv.in,jump=jumps.in,rvf="pareto",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.bdbj<- nlminb(runif(9*d + d*d,0,0.1), rv=bpv.in,jump=jumps.in,rvf="burr",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))
fit.full.wdbj<- nlminb(runif(9*d + d*d,0,0.1), rv=bpv.in,jump=jumps.in,rvf="weibull",spec="full", dyn=T, VARtime.nlik,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,opt=T,control=list(eval.max=20000,iter.max=20000,trace=100))

fit.full.pdbj$objective
fit.full.bdbj$objective
fit.full.wdbj$objective

hess.full.pdbj<-hessb(VARtime.nlik,fit.full.pdbj$par,  rv=bpv.in,jump=jumps.in,rvf="pareto", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.pdbj<-SEfromHessian(hess.full.pdbj)

hess.full.bdbj<-hessb(VARtime.nlik,fit.full.bdbj$par,  rv=bpv.in,jump=jumps.in,rvf="burr", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.bdbj<-SEfromHessian(hess.full.bdbj)

hess.full.wdbj<-hessb(VARtime.nlik,fit.full.wdbj$par,  rv=bpv.in,jump=jumps.in,rvf="weibull", spec="full", dyn=T,zt=zt.in,umbral=umbral.in, I1=I1.in,I0=I0.in,ep = 1e-09)
se.full.wdbj<-SEfromHessian(hess.full.wdbj)

out.full.pdbj<-VARtime.nlik(fit.full.pdbj$par,rv=bpv.out,jump=jumps.out,rvf="pareto", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.bdbj<-VARtime.nlik(fit.full.bdbj$par,rv=bpv.out,jump=jumps.out,rvf="burr", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)
out.full.wdbj<-VARtime.nlik(fit.full.wdbj$par,rv=bpv.out,jump=jumps.out,rvf="weibull", spec="full", dyn=T,zt=zt.out,umbral=umbral.out, I1=I1.out,I0=I0.out,opt=F)


# Getting parameter estimates and standar errors
#pareto
Ap<-matrix(fit.full.pdbj$par[1:(d*d)],d,d)
ap<-matrix(fit.full.pdbj$par[(d*d+1):length(fit.full.pdbj$par)],8,d)
length(Ap)
length(ap)
Ap
ap
#pareto s.e
se.Ap<-matrix(se.full.pdbj[1:(d*d)],d,d)
se.ap<-matrix(se.full.pdbj[(d*d+1):length(se.full.pdbj)],8,d)
length(se.Ap)
length(se.ap)
se.Ap
se.ap

#burr
#Ab<-matrix(fit.full.bdbj$par[1:(d*d)],d,d)
#ab<-matrix(fit.full.bdbj$par[(d*d+1):length(fit.full.bdbj$par)],9,d)
#length(Ab)
#length(ab)
#Ab
#ab
#burr s.e
#se.Ab<-matrix(se.full.bdbj[1:(d*d)],d,d)
#se.ab<-matrix(se.full.bdbj[(d*d+1):length(se.full.bdbj)],9,d)
#length(se.Ab)
#length(se.ab)
#se.Ab
#se.ab

#weibull
#Aw<-matrix(fit.full.wdbj$par[1:(d*d)],d,d)
#aw<-matrix(fit.full.wdbj$par[(d*d+1):length(fit.full.wdbj$par)],9,d)
#length(Aw)
#length(aw)
#Aw
#aw
#weibull s.e
#se.Aw<-matrix(se.full.wdbj[1:(d*d)],d,d)
#se.aw<-matrix(se.full.wdbj[(d*d+1):length(se.full.wdbj)],9,d)
#length(se.Aw)
#length(se.aw)
#se.Aw
#se.aw


###################################################################################################################################
# Akaike information criteria
###################################################################################################################################

AIC.VARtime<-function(fit){
  
  2*length(fit$par) +2*fit$objective
}
AIC.VARtime(fit.diag.p)


# Diagonal and non-dynamic (p pareto, b burr, w weibull) models
sapply(list(fit.diag.p,fit.diag.b,fit.diag.w),AIC.VARtime)

# Diagonal and dynamic ( p pareto, b burr, w weibull) models
sapply(list(fit.diag.pd,fit.diag.bd,fit.diag.wd),AIC.VARtime)

# Diagonal and dynamic models( p pareto, b burr, w weibull) with standard realized measure 
sapply(list(fit.diag.pdr,fit.diag.bdr,fit.diag.wdr),AIC.VARtime)

# Diagonal and dynamic models( p pareto, b burr, w weibull) with negative RSV 
sapply(list(fit.diag.pdrs,fit.diag.bdrs,fit.diag.wdrs),AIC.VARtime)

# Diagonal and dynamic models( p pareto, b burr, w weibull) with negative and positive RSV 
sapply(list(fit.diag.pdpn,fit.diag.bdpn,fit.diag.wdpn),AIC.VARtime)

# Diagonal and dynamic models( p pareto, b burr, w weibull) with Bipower variation and jumps
sapply(list(fit.diag.pdbj,fit.diag.bdbj,fit.diag.wdbj),AIC.VARtime)

# Full (non-diagonal) and non-dynamic (m massaci, p pareto, b burr, w weibull) models
sapply(list(fit.full.p,fit.full.b,fit.full.w),AIC.VARtime)

# Full (non-diagonal) and dynamic (p pareto, b burr, w weibull) models
sapply(list(fit.full.pd,fit.full.bd,fit.full.wd),AIC.VARtime)

# Full (non-diagonal) and dynamic (p pareto, b burr, w weibull) with standard realized measure
sapply(list(fit.full.pdr,fit.full.bdr,fit.full.wdr),AIC.VARtime)

# Full (non-diagonal) and dynamic (p pareto, b burr, w weibull) with negative RSV
sapply(list(fit.full.pdrs,fit.full.bdrs,fit.full.wdrs),AIC.VARtime)

# Full (non-diagonal) and dynamic (p pareto, b burr, w weibull) with negative and positive RSV 
sapply(list(fit.full.pdpn,fit.full.bdpn,fit.full.wdpn),AIC.VARtime)

# Full (non-diagonal) and dynamic (p pareto, b burr, w weibull) with Bipower variation and jumps
sapply(list(fit.full.pdbj,fit.full.bdbj,fit.full.wdbj),AIC.VARtime)




###################################################################################################################################
# Accuracy test of Value-at-risk (VaR) and Expected shortfall (ES) and lost function for full dynamic models 
###################################################################################################################################

n.out<-dim(zt.out)[1]
p<-c(0.95,0.975,0.99,0.999)


#sink("./Accuracy_resultados01.csv", append = T, split= T)
#sink("./Accuracy_resultados01.txt", append = T, split= T)

#------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic PARETO model

#VaR full pd
VAR.full.pd.out<-matrix(0,n.out,d*length(p))
VAR.full.pd.out
dim(VAR.full.pd.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.pd.out[,(i-1)*d+j]<-VaR(p[i], out.full.pd$phit[,j], out.full.pd$xit[,j],  out.full.pd$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.pd.out.0.95 <-VAR.full.pd.out[,1:d]
VAR.full.pd.out.0.975<-VAR.full.pd.out[,(d+1):(2*d)]
VAR.full.pd.out.0.99 <-VAR.full.pd.out[,(2*d+1):(3*d)]
VAR.full.pd.out.0.999<-VAR.full.pd.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.pd.out.0.95[100:dim(VAR.full.pd.out.0.95)[1],1:d])
plot.ts(VAR.full.pd.out.0.975[100:dim(VAR.full.pd.out.0.975)[1],1:d])
plot.ts(VAR.full.pd.out.0.99[100:dim(VAR.full.pd.out.0.99)[1],1:d])
plot.ts(VAR.full.pd.out.0.999[100:dim(VAR.full.pd.out.0.999)[1],1:d])

#ES full pd
ES.full.pd.out<-matrix(0,n.out,d*length(p))
ES.full.pd.out
dim(ES.full.pd.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.pd.out[,(i-1)*d+j]<-ES(p[i], out.full.pd$phit[,j], out.full.pd$xit[,j],  out.full.pd$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.pd.out.0.95 <-ES.full.pd.out[,1:d]
ES.full.pd.out.0.975<-ES.full.pd.out[,(d+1):(2*d)]
ES.full.pd.out.0.99 <-ES.full.pd.out[,(2*d+1):(3*d)]
ES.full.pd.out.0.999<-ES.full.pd.out[,(3*d+1):(4*d)]
plot.ts(ES.full.pd.out.0.95[100:dim(ES.full.pd.out.0.95)[1],1:d])
plot.ts(ES.full.pd.out.0.975[100:dim(ES.full.pd.out.0.975)[1],1:d])
plot.ts(ES.full.pd.out.0.99[100:dim(ES.full.pd.out.0.99)[1],1:d])
plot.ts(ES.full.pd.out.0.999[100:dim(ES.full.pd.out.0.999)[1],1:d])


## Accuracy test full pd
VAR.test.full.pd<-matrix(0,n.out,(length(p))*d)
ES.test.full.pd<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.pd[, (j-1)*length(p)+i]<-VAR.full.pd.out[,(i-1)*d+j]
    ES.test.full.pd[, (j-1)*length(p)+i]<-ES.full.pd.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.pd[2:dim(VAR.test.full.pd)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.pd[2:dim(ES.test.full.pd)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.pd[2:dim(ES.test.full.pd)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full pd
full.pd.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.pd.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.pd.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.pd.out[,(i-1)*d+j],p[i])
    full.pd.out.nol2[,(i-1)*d+j][is.na(full.pd.out.nol2[,(i-1)*d+j])]<-0
    full.pd.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.pd.out[,(i-1)*d+j],ES.full.pd.out[,(i-1)*d+j],p[i])
    full.pd.out.nol6[,(i-1)*d+j][is.na(full.pd.out.nol6[,(i-1)*d+j])]<-0
    full.pd.out.nol6[,(i-1)*d+j][is.infinite(full.pd.out.nol6[,(i-1)*d+j])]<-0
  }
}

#--------------------------------------------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic PARETO model with standard realized measure

#VaR full pdr
VAR.full.pdr.out<-matrix(0,n.out,d*length(p))
VAR.full.pdr.out
dim(VAR.full.pdr.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.pdr.out[,(i-1)*d+j]<-VaR(p[i], out.full.pdr$phit[,j], out.full.pdr$xit[,j],  out.full.pdr$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.pdr.out.0.95 <-VAR.full.pdr.out[,1:d]
VAR.full.pdr.out.0.975<-VAR.full.pdr.out[,(d+1):(2*d)]
VAR.full.pdr.out.0.99 <-VAR.full.pdr.out[,(2*d+1):(3*d)]
VAR.full.pdr.out.0.999<-VAR.full.pdr.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.pdr.out.0.95[100:dim(VAR.full.pdr.out.0.95)[1],1:d])
plot.ts(VAR.full.pdr.out.0.975[100:dim(VAR.full.pdr.out.0.975)[1],1:d])
plot.ts(VAR.full.pdr.out.0.99[100:dim(VAR.full.pdr.out.0.99)[1],1:d])
plot.ts(VAR.full.pdr.out.0.999[100:dim(VAR.full.pdr.out.0.999)[1],1:d])

#ES full pdr
ES.full.pdr.out<-matrix(0,n.out,d*length(p))
ES.full.pdr.out
dim(ES.full.pdr.out)

for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.pdr.out[,(i-1)*d+j]<-ES(p[i], out.full.pdr$phit[,j], out.full.pdr$xit[,j],  out.full.pdr$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.pdr.out.0.95 <-ES.full.pdr.out[,1:d]
ES.full.pdr.out.0.975<-ES.full.pdr.out[,(d+1):(2*d)]
ES.full.pdr.out.0.99 <-ES.full.pdr.out[,(2*d+1):(3*d)]
ES.full.pdr.out.0.999<-ES.full.pdr.out[,(3*d+1):(4*d)]
plot.ts(ES.full.pdr.out.0.95[100:dim(ES.full.pdr.out.0.95)[1],1:d])
plot.ts(ES.full.pdr.out.0.975[100:dim(ES.full.pdr.out.0.975)[1],1:d])
plot.ts(ES.full.pdr.out.0.99[100:dim(ES.full.pdr.out.0.99)[1],1:d])
plot.ts(ES.full.pdr.out.0.999[100:dim(ES.full.pdr.out.0.999)[1],1:d])


## Accuracy test full pdr
VAR.test.full.pdr<-matrix(0,n.out,(length(p))*d)
ES.test.full.pdr<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.pdr[, (j-1)*length(p)+i]<-VAR.full.pdr.out[,(i-1)*d+j]
    ES.test.full.pdr[, (j-1)*length(p)+i]<-ES.full.pdr.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.pdr[2:dim(VAR.test.full.pdr)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.pdr[2:dim(ES.test.full.pdr)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.pdr[2:dim(ES.test.full.pdr)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full pdr
full.pdr.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.pdr.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.pdr.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.pdr.out[,(i-1)*d+j],p[i])
    full.pdr.out.nol2[,(i-1)*d+j][is.na(full.pdr.out.nol2[,(i-1)*d+j])]<-0
    full.pdr.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.pdr.out[,(i-1)*d+j],ES.full.pdr.out[,(i-1)*d+j],p[i])
    full.pdr.out.nol6[,(i-1)*d+j][is.na(full.pdr.out.nol6[,(i-1)*d+j])]<-0
    full.pdr.out.nol6[,(i-1)*d+j][is.infinite(full.pdr.out.nol6[,(i-1)*d+j])]<-0
  }
}
#------------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic PARETO model with negative RSV

#VaR full pdrs
VAR.full.pdrs.out<-matrix(0,n.out,d*length(p))
VAR.full.pdrs.out
dim(VAR.full.pdrs.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.pdrs.out[,(i-1)*d+j]<-VaR(p[i], out.full.pdrs$phit[,j], out.full.pdrs$xit[,j],  out.full.pdrs$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.pdrs.out.0.95 <-VAR.full.pdrs.out[,1:d]
VAR.full.pdrs.out.0.975<-VAR.full.pdrs.out[,(d+1):(2*d)]
VAR.full.pdrs.out.0.99 <-VAR.full.pdrs.out[,(2*d+1):(3*d)]
VAR.full.pdrs.out.0.999<-VAR.full.pdrs.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.pdrs.out.0.95[100:dim(VAR.full.pdrs.out.0.95)[1],1:d])
plot.ts(VAR.full.pdrs.out.0.975[100:dim(VAR.full.pdrs.out.0.975)[1],1:d])
plot.ts(VAR.full.pdrs.out.0.99[100:dim(VAR.full.pdrs.out.0.99)[1],1:d])
plot.ts(VAR.full.pdrs.out.0.999[100:dim(VAR.full.pdrs.out.0.999)[1],1:d])

#ES full pdrs
ES.full.pdrs.out<-matrix(0,n.out,d*length(p))
ES.full.pdrs.out
dim(ES.full.pdrs.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.pdrs.out[,(i-1)*d+j]<-ES(p[i], out.full.pdrs$phit[,j], out.full.pdrs$xit[,j],  out.full.pdrs$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.pdrs.out.0.95 <-ES.full.pdrs.out[,1:d]
ES.full.pdrs.out.0.975<-ES.full.pdrs.out[,(d+1):(2*d)]
ES.full.pdrs.out.0.99 <-ES.full.pdrs.out[,(2*d+1):(3*d)]
ES.full.pdrs.out.0.999<-ES.full.pdrs.out[,(3*d+1):(4*d)]
plot.ts(ES.full.pdrs.out.0.95[100:dim(ES.full.pdrs.out.0.95)[1],1:d])
plot.ts(ES.full.pdrs.out.0.975[100:dim(ES.full.pdrs.out.0.975)[1],1:d])
plot.ts(ES.full.pdrs.out.0.99[100:dim(ES.full.pdrs.out.0.99)[1],1:d])
plot.ts(ES.full.pdrs.out.0.999[100:dim(ES.full.pdrs.out.0.999)[1],1:d])


## Accuracy test full pdrs
VAR.test.full.pdrs<-matrix(0,n.out,(length(p))*d)
ES.test.full.pdrs<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.pdrs[, (j-1)*length(p)+i]<-VAR.full.pdrs.out[,(i-1)*d+j]
    ES.test.full.pdrs[, (j-1)*length(p)+i]<-ES.full.pdrs.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.pdrs[2:dim(VAR.test.full.pdrs)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.pdrs[2:dim(ES.test.full.pdrs)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.pdrs[2:dim(ES.test.full.pdrs)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full pdrs
full.pdrs.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.pdrs.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.pdrs.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.pdrs.out[,(i-1)*d+j],p[i])
    full.pdrs.out.nol2[,(i-1)*d+j][is.na(full.pdrs.out.nol2[,(i-1)*d+j])]<-0
    full.pdrs.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.pdrs.out[,(i-1)*d+j],ES.full.pdrs.out[,(i-1)*d+j],p[i])
    full.pdrs.out.nol6[,(i-1)*d+j][is.na(full.pdrs.out.nol6[,(i-1)*d+j])]<-0
    full.pdrs.out.nol6[,(i-1)*d+j][is.infinite(full.pdrs.out.nol6[,(i-1)*d+j])]<-0
  }
}

#-----------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic PARETO model with negative and positive RSV


#VaR full pdpn
VAR.full.pdpn.out<-matrix(0,n.out,d*length(p))
VAR.full.pdpn.out
dim(VAR.full.pdpn.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.pdpn.out[,(i-1)*d+j]<-VaR(p[i], out.full.pdpn$phit[,j], out.full.pdpn$xit[,j],  out.full.pdpn$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.pdpn.out.0.95 <-VAR.full.pdpn.out[,1:d]
VAR.full.pdpn.out.0.975<-VAR.full.pdpn.out[,(d+1):(2*d)]
VAR.full.pdpn.out.0.99 <-VAR.full.pdpn.out[,(2*d+1):(3*d)]
VAR.full.pdpn.out.0.999<-VAR.full.pdpn.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.pdpn.out.0.95[100:dim(VAR.full.pdpn.out.0.95)[1],1:d])
plot.ts(VAR.full.pdpn.out.0.975[100:dim(VAR.full.pdpn.out.0.975)[1],1:d])
plot.ts(VAR.full.pdpn.out.0.99[100:dim(VAR.full.pdpn.out.0.99)[1],1:d])
plot.ts(VAR.full.pdpn.out.0.999[100:dim(VAR.full.pdpn.out.0.999)[1],1:d])

#ES full pdpn
ES.full.pdpn.out<-matrix(0,n.out,d*length(p))
ES.full.pdpn.out
dim(ES.full.pdpn.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.pdpn.out[,(i-1)*d+j]<-ES(p[i], out.full.pdpn$phit[,j], out.full.pdpn$xit[,j],  out.full.pdpn$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.pdpn.out.0.95 <-ES.full.pdpn.out[,1:d]
ES.full.pdpn.out.0.975<-ES.full.pdpn.out[,(d+1):(2*d)]
ES.full.pdpn.out.0.99 <-ES.full.pdpn.out[,(2*d+1):(3*d)]
ES.full.pdpn.out.0.999<-ES.full.pdpn.out[,(3*d+1):(4*d)]
plot.ts(ES.full.pdpn.out.0.95[100:dim(ES.full.pdpn.out.0.95)[1],1:d])
plot.ts(ES.full.pdpn.out.0.975[100:dim(ES.full.pdpn.out.0.975)[1],1:d])
plot.ts(ES.full.pdpn.out.0.99[100:dim(ES.full.pdpn.out.0.99)[1],1:d])
plot.ts(ES.full.pdpn.out.0.999[100:dim(ES.full.pdpn.out.0.999)[1],1:d])


## Accuracy test full pdpn
VAR.test.full.pdpn<-matrix(0,n.out,(length(p))*d)
ES.test.full.pdpn<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.pdpn[, (j-1)*length(p)+i]<-VAR.full.pdpn.out[,(i-1)*d+j]
    ES.test.full.pdpn[, (j-1)*length(p)+i]<-ES.full.pdpn.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.pdpn[2:dim(VAR.test.full.pdpn)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.pdpn[2:dim(ES.test.full.pdpn)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.pdpn[2:dim(ES.test.full.pdpn)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full pdpn
full.pdpn.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.pdpn.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.pdpn.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.pdpn.out[,(i-1)*d+j],p[i])
    full.pdpn.out.nol2[,(i-1)*d+j][is.na(full.pdpn.out.nol2[,(i-1)*d+j])]<-0
    full.pdpn.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.pdpn.out[,(i-1)*d+j],ES.full.pdpn.out[,(i-1)*d+j],p[i])
    full.pdpn.out.nol6[,(i-1)*d+j][is.na(full.pdpn.out.nol6[,(i-1)*d+j])]<-0
    full.pdpn.out.nol6[,(i-1)*d+j][is.infinite(full.pdpn.out.nol6[,(i-1)*d+j])]<-0
  }
}

#---------------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic PARETO model with bipower and jumps

#VaR full pdbj
VAR.full.pdbj.out<-matrix(0,n.out,d*length(p))
VAR.full.pdbj.out
dim(VAR.full.pdbj.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.pdbj.out[,(i-1)*d+j]<-VaR(p[i], out.full.pdbj$phit[,j], out.full.pdbj$xit[,j],  out.full.pdbj$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.pdbj.out.0.95 <-VAR.full.pdbj.out[,1:d]
VAR.full.pdbj.out.0.975<-VAR.full.pdbj.out[,(d+1):(2*d)]
VAR.full.pdbj.out.0.99 <-VAR.full.pdbj.out[,(2*d+1):(3*d)]
VAR.full.pdbj.out.0.999<-VAR.full.pdbj.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.pdbj.out.0.95[100:dim(VAR.full.pdbj.out.0.95)[1],1:d])
plot.ts(VAR.full.pdbj.out.0.975[100:dim(VAR.full.pdbj.out.0.975)[1],1:d])
plot.ts(VAR.full.pdbj.out.0.99[100:dim(VAR.full.pdbj.out.0.99)[1],1:d])
plot.ts(VAR.full.pdbj.out.0.999[100:dim(VAR.full.pdbj.out.0.999)[1],1:d])

#ES full pdbj
ES.full.pdbj.out<-matrix(0,n.out,d*length(p))
ES.full.pdbj.out
dim(ES.full.pdbj.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.pdbj.out[,(i-1)*d+j]<-ES(p[i], out.full.pdbj$phit[,j], out.full.pdbj$xit[,j],  out.full.pdbj$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.pdbj.out.0.95 <-ES.full.pdbj.out[,1:d]
ES.full.pdbj.out.0.975<-ES.full.pdbj.out[,(d+1):(2*d)]
ES.full.pdbj.out.0.99 <-ES.full.pdbj.out[,(2*d+1):(3*d)]
ES.full.pdbj.out.0.999<-ES.full.pdbj.out[,(3*d+1):(4*d)]
plot.ts(ES.full.pdbj.out.0.95[100:dim(ES.full.pdbj.out.0.95)[1],1:d])
plot.ts(ES.full.pdbj.out.0.975[100:dim(ES.full.pdbj.out.0.975)[1],1:d])
plot.ts(ES.full.pdbj.out.0.99[100:dim(ES.full.pdbj.out.0.99)[1],1:d])
plot.ts(ES.full.pdbj.out.0.999[100:dim(ES.full.pdbj.out.0.999)[1],1:d])


## Accuracy test full pdbj
VAR.test.full.pdbj<-matrix(0,n.out,(length(p))*d)
ES.test.full.pdbj<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.pdbj[, (j-1)*length(p)+i]<-VAR.full.pdbj.out[,(i-1)*d+j]
    ES.test.full.pdbj[, (j-1)*length(p)+i]<-ES.full.pdbj.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.pdbj[2:dim(VAR.test.full.pdbj)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.pdbj[2:dim(ES.test.full.pdbj)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.pdbj[2:dim(ES.test.full.pdbj)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full pdbj
full.pdbj.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.pdbj.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.pdbj.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.pdbj.out[,(i-1)*d+j],p[i])
    full.pdbj.out.nol2[,(i-1)*d+j][is.na(full.pdbj.out.nol2[,(i-1)*d+j])]<-0
    full.pdbj.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.pdbj.out[,(i-1)*d+j],ES.full.pdbj.out[,(i-1)*d+j],p[i])
    full.pdbj.out.nol6[,(i-1)*d+j][is.na(full.pdbj.out.nol6[,(i-1)*d+j])]<-0
    full.pdbj.out.nol6[,(i-1)*d+j][is.infinite(full.pdbj.out.nol6[,(i-1)*d+j])]<-0
  }
}

#-------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic BURR model

#VaR full bd
VAR.full.bd.out<-matrix(0,n.out,d*length(p))
VAR.full.bd.out
dim(VAR.full.bd.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.bd.out[,(i-1)*d+j]<-VaR(p[i], out.full.bd$phit[,j], out.full.bd$xit[,j],  out.full.bd$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.bd.out.0.95 <-VAR.full.bd.out[,1:d]
VAR.full.bd.out.0.975<-VAR.full.bd.out[,(d+1):(2*d)]
VAR.full.bd.out.0.99 <-VAR.full.bd.out[,(2*d+1):(3*d)]
VAR.full.bd.out.0.999<-VAR.full.bd.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.bd.out.0.95[100:dim(VAR.full.bd.out.0.95)[1],1:d])
plot.ts(VAR.full.bd.out.0.975[100:dim(VAR.full.bd.out.0.975)[1],1:d])
plot.ts(VAR.full.bd.out.0.99[100:dim(VAR.full.bd.out.0.99)[1],1:d])
plot.ts(VAR.full.bd.out.0.999[100:dim(VAR.full.bd.out.0.999)[1],1:d])

#ES full bd
ES.full.bd.out<-matrix(0,n.out,d*length(p))
ES.full.bd.out
dim(ES.full.bd.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.bd.out[,(i-1)*d+j]<-ES(p[i], out.full.bd$phit[,j], out.full.bd$xit[,j],  out.full.bd$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.bd.out.0.95 <-ES.full.bd.out[,1:d]
ES.full.bd.out.0.975<-ES.full.bd.out[,(d+1):(2*d)]
ES.full.bd.out.0.99 <-ES.full.bd.out[,(2*d+1):(3*d)]
ES.full.bd.out.0.999<-ES.full.bd.out[,(3*d+1):(4*d)]
plot.ts(ES.full.bd.out.0.95[100:dim(ES.full.bd.out.0.95)[1],1:d])
plot.ts(ES.full.bd.out.0.975[100:dim(ES.full.bd.out.0.975)[1],1:d])
plot.ts(ES.full.bd.out.0.99[100:dim(ES.full.bd.out.0.99)[1],1:d])
plot.ts(ES.full.bd.out.0.999[100:dim(ES.full.bd.out.0.999)[1],1:d])


## Accuracy test full bd
VAR.test.full.bd<-matrix(0,n.out,(length(p))*d)
ES.test.full.bd<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.bd[, (j-1)*length(p)+i]<-VAR.full.bd.out[,(i-1)*d+j]
    ES.test.full.bd[, (j-1)*length(p)+i]<-ES.full.bd.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.bd[2:dim(VAR.test.full.bd)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.bd[2:dim(ES.test.full.bd)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.bd[2:dim(ES.test.full.bd)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full bd
full.bd.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.bd.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.bd.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.bd.out[,(i-1)*d+j],p[i])
    full.bd.out.nol2[,(i-1)*d+j][is.na(full.bd.out.nol2[,(i-1)*d+j])]<-0
    full.bd.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.bd.out[,(i-1)*d+j],ES.full.bd.out[,(i-1)*d+j],p[i])
    full.bd.out.nol6[,(i-1)*d+j][is.na(full.bd.out.nol6[,(i-1)*d+j])]<-0
    full.bd.out.nol6[,(i-1)*d+j][is.infinite(full.bd.out.nol6[,(i-1)*d+j])]<-0
  }
}


#-------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic BURR model with standard realized measure

#VaR full bdr

VAR.full.bdr.out<-matrix(0,n.out,d*length(p))
VAR.full.bdr.out
dim(VAR.full.bdr.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.bdr.out[,(i-1)*d+j]<-VaR(p[i], out.full.bdr$phit[,j], out.full.bdr$xit[,j],  out.full.bdr$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.bdr.out.0.95 <-VAR.full.bdr.out[,1:d]
VAR.full.bdr.out.0.975<-VAR.full.bdr.out[,(d+1):(2*d)]
VAR.full.bdr.out.0.99 <-VAR.full.bdr.out[,(2*d+1):(3*d)]
VAR.full.bdr.out.0.999<-VAR.full.bdr.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.bdr.out.0.95[100:dim(VAR.full.bdr.out.0.95)[1],1:d])
plot.ts(VAR.full.bdr.out.0.975[100:dim(VAR.full.bdr.out.0.975)[1],1:d])
plot.ts(VAR.full.bdr.out.0.99[100:dim(VAR.full.bdr.out.0.99)[1],1:d])
plot.ts(VAR.full.bdr.out.0.999[100:dim(VAR.full.bdr.out.0.999)[1],1:d])

#ES full bdr
ES.full.bdr.out<-matrix(0,n.out,d*length(p))
ES.full.bdr.out
dim(ES.full.bdr.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.bdr.out[,(i-1)*d+j]<-ES(p[i], out.full.bdr$phit[,j], out.full.bdr$xit[,j],  out.full.bdr$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.bdr.out.0.95 <-ES.full.pdr.out[,1:d]
ES.full.bdr.out.0.975<-ES.full.pdr.out[,(d+1):(2*d)]
ES.full.bdr.out.0.99 <-ES.full.pdr.out[,(2*d+1):(3*d)]
ES.full.bdr.out.0.999<-ES.full.pdr.out[,(3*d+1):(4*d)]
plot.ts(ES.full.bdr.out.0.95[100:dim(ES.full.bdr.out.0.95)[1],1:d])
plot.ts(ES.full.bdr.out.0.975[100:dim(ES.full.bdr.out.0.975)[1],1:d])
plot.ts(ES.full.bdr.out.0.99[100:dim(ES.full.bdr.out.0.99)[1],1:d])
plot.ts(ES.full.bdr.out.0.999[100:dim(ES.full.bdr.out.0.999)[1],1:d])


## Accuracy test full bdr
VAR.test.full.bdr<-matrix(0,n.out,(length(p))*d)
ES.test.full.bdr<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.bdr[, (j-1)*length(p)+i]<-VAR.full.bdr.out[,(i-1)*d+j]
    ES.test.full.bdr[, (j-1)*length(p)+i]<-ES.full.bdr.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.bdr[2:dim(VAR.test.full.bdr)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.bdr[2:dim(ES.test.full.bdr)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.bdr[2:dim(ES.test.full.bdr)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full bdr
full.bdr.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.bdr.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.bdr.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.bdr.out[,(i-1)*d+j],p[i])
    full.bdr.out.nol2[,(i-1)*d+j][is.na(full.bdr.out.nol2[,(i-1)*d+j])]<-0
    full.bdr.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.bdr.out[,(i-1)*d+j],ES.full.bdr.out[,(i-1)*d+j],p[i])
    full.bdr.out.nol6[,(i-1)*d+j][is.na(full.bdr.out.nol6[,(i-1)*d+j])]<-0
    full.bdr.out.nol6[,(i-1)*d+j][is.infinite(full.bdr.out.nol6[,(i-1)*d+j])]<-0
  }
}

#------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic BURR model with negative RSV

#VaR full bdrs
VAR.full.bdrs.out<-matrix(0,n.out,d*length(p))
VAR.full.bdrs.out
dim(VAR.full.bdrs.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.bdrs.out[,(i-1)*d+j]<-VaR(p[i], out.full.bdrs$phit[,j], out.full.bdrs$xit[,j],  out.full.bdrs$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.bdrs.out.0.95 <-VAR.full.bdrs.out[,1:d]
VAR.full.bdrs.out.0.975<-VAR.full.bdrs.out[,(d+1):(2*d)]
VAR.full.bdrs.out.0.99 <-VAR.full.bdrs.out[,(2*d+1):(3*d)]
VAR.full.bdrs.out.0.999<-VAR.full.bdrs.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.bdrs.out.0.95[100:dim(VAR.full.bdrs.out.0.95)[1],1:d])
plot.ts(VAR.full.bdrs.out.0.975[100:dim(VAR.full.bdrs.out.0.975)[1],1:d])
plot.ts(VAR.full.bdrs.out.0.99[100:dim(VAR.full.bdrs.out.0.99)[1],1:d])
plot.ts(VAR.full.bdrs.out.0.999[100:dim(VAR.full.bdrs.out.0.999)[1],1:d])

#ES full bdrs
ES.full.bdrs.out<-matrix(0,n.out,d*length(p))
ES.full.bdrs.out
dim(ES.full.bdrs.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.bdrs.out[,(i-1)*d+j]<-ES(p[i], out.full.bdrs$phit[,j], out.full.bdrs$xit[,j],  out.full.bdrs$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.bdrs.out.0.95 <-ES.full.bdrs.out[,1:d]
ES.full.bdrs.out.0.975<-ES.full.bdrs.out[,(d+1):(2*d)]
ES.full.bdrs.out.0.99 <-ES.full.bdrs.out[,(2*d+1):(3*d)]
ES.full.bdrs.out.0.999<-ES.full.bdrs.out[,(3*d+1):(4*d)]
plot.ts(ES.full.bdrs.out.0.95[100:dim(ES.full.bdrs.out.0.95)[1],1:d])
plot.ts(ES.full.bdrs.out.0.975[100:dim(ES.full.bdrs.out.0.975)[1],1:d])
plot.ts(ES.full.bdrs.out.0.99[100:dim(ES.full.bdrs.out.0.99)[1],1:d])
plot.ts(ES.full.bdrs.out.0.999[100:dim(ES.full.bdrs.out.0.999)[1],1:d])


## Accuracy test full bdrs
VAR.test.full.bdrs<-matrix(0,n.out,(length(p))*d)
ES.test.full.bdrs<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.bdrs[, (j-1)*length(p)+i]<-VAR.full.bdrs.out[,(i-1)*d+j]
    ES.test.full.bdrs[, (j-1)*length(p)+i]<-ES.full.bdrs.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.bdrs[2:dim(VAR.test.full.bdrs)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.bdrs[2:dim(ES.test.full.bdrs)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.bdrs[2:dim(ES.test.full.bdrs)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full bdrs
full.bdrs.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.bdrs.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.bdrs.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.bdrs.out[,(i-1)*d+j],p[i])
    full.bdrs.out.nol2[,(i-1)*d+j][is.na(full.bdrs.out.nol2[,(i-1)*d+j])]<-0
    full.bdrs.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.bdrs.out[,(i-1)*d+j],ES.full.bdrs.out[,(i-1)*d+j],p[i])
    full.bdrs.out.nol6[,(i-1)*d+j][is.na(full.bdrs.out.nol6[,(i-1)*d+j])]<-0
    full.bdrs.out.nol6[,(i-1)*d+j][is.infinite(full.bdrs.out.nol6[,(i-1)*d+j])]<-0
  }
}

#-----------------------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic BURR model with positive and negative RSV

#VaR full bdpn
VAR.full.bdpn.out<-matrix(0,n.out,d*length(p))
VAR.full.bdpn.out
dim(VAR.full.bdpn.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.bdpn.out[,(i-1)*d+j]<-VaR(p[i], out.full.bdpn$phit[,j], out.full.bdpn$xit[,j],  out.full.bdpn$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.bdpn.out.0.95 <-VAR.full.bdpn.out[,1:d]
VAR.full.bdpn.out.0.975<-VAR.full.bdpn.out[,(d+1):(2*d)]
VAR.full.bdpn.out.0.99 <-VAR.full.bdpn.out[,(2*d+1):(3*d)]
VAR.full.bdpn.out.0.999<-VAR.full.bdpn.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.bdpn.out.0.95[100:dim(VAR.full.bdpn.out.0.95)[1],1:d])
plot.ts(VAR.full.bdpn.out.0.975[100:dim(VAR.full.bdpn.out.0.975)[1],1:d])
plot.ts(VAR.full.bdpn.out.0.99[100:dim(VAR.full.bdpn.out.0.99)[1],1:d])
plot.ts(VAR.full.bdpn.out.0.999[100:dim(VAR.full.bdpn.out.0.999)[1],1:d])

#ES full bdpn
ES.full.bdpn.out<-matrix(0,n.out,d*length(p))
ES.full.bdpn.out
dim(ES.full.bdpn.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.bdpn.out[,(i-1)*d+j]<-ES(p[i], out.full.bdpn$phit[,j], out.full.bdpn$xit[,j],  out.full.bdpn$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.bdpn.out.0.95 <-ES.full.bdpn.out[,1:d]
ES.full.bdpn.out.0.975<-ES.full.bdpn.out[,(d+1):(2*d)]
ES.full.bdpn.out.0.99 <-ES.full.bdpn.out[,(2*d+1):(3*d)]
ES.full.bdpn.out.0.999<-ES.full.bdpn.out[,(3*d+1):(4*d)]
plot.ts(ES.full.bdpn.out.0.95[100:dim(ES.full.bdpn.out.0.95)[1],1:d])
plot.ts(ES.full.bdpn.out.0.975[100:dim(ES.full.bdpn.out.0.975)[1],1:d])
plot.ts(ES.full.bdpn.out.0.99[100:dim(ES.full.bdpn.out.0.99)[1],1:d])
plot.ts(ES.full.bdpn.out.0.999[100:dim(ES.full.bdpn.out.0.999)[1],1:d])


## Accuracy test full bdpn
VAR.test.full.bdpn<-matrix(0,n.out,(length(p))*d)
ES.test.full.bdpn<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.bdpn[, (j-1)*length(p)+i]<-VAR.full.bdpn.out[,(i-1)*d+j]
    ES.test.full.bdpn[, (j-1)*length(p)+i]<-ES.full.bdpn.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.bdpn[2:dim(VAR.test.full.bdpn)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.bdpn[2:dim(ES.test.full.bdpn)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.bdpn[2:dim(ES.test.full.bdpn)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full bdpn
full.bdpn.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.bdpn.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.bdpn.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.bdpn.out[,(i-1)*d+j],p[i])
    full.bdpn.out.nol2[,(i-1)*d+j][is.na(full.bdpn.out.nol2[,(i-1)*d+j])]<-0
    full.bdpn.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.bdpn.out[,(i-1)*d+j],ES.full.bdpn.out[,(i-1)*d+j],p[i])
    full.bdpn.out.nol6[,(i-1)*d+j][is.na(full.bdpn.out.nol6[,(i-1)*d+j])]<-0
    full.bdpn.out.nol6[,(i-1)*d+j][is.infinite(full.bdpn.out.nol6[,(i-1)*d+j])]<-0
  }
}

#---------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic BURR model with bipower and jumps

#VaR full bdbj
VAR.full.bdbj.out<-matrix(0,n.out,d*length(p))
VAR.full.bdbj.out
dim(VAR.full.bdbj.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.bdbj.out[,(i-1)*d+j]<-VaR(p[i], out.full.bdbj$phit[,j], out.full.bdbj$xit[,j],  out.full.bdbj$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.bdbj.out.0.95 <-VAR.full.bdbj.out[,1:d]
VAR.full.bdbj.out.0.975<-VAR.full.bdbj.out[,(d+1):(2*d)]
VAR.full.bdbj.out.0.99 <-VAR.full.bdbj.out[,(2*d+1):(3*d)]
VAR.full.bdbj.out.0.999<-VAR.full.bdbj.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.bdbj.out.0.95[100:dim(VAR.full.bdbj.out.0.95)[1],1:d])
plot.ts(VAR.full.bdbj.out.0.975[100:dim(VAR.full.bdbj.out.0.975)[1],1:d])
plot.ts(VAR.full.bdbj.out.0.99[100:dim(VAR.full.bdbj.out.0.99)[1],1:d])
plot.ts(VAR.full.bdbj.out.0.999[100:dim(VAR.full.bdbj.out.0.999)[1],1:d])

#ES full bdbj
ES.full.bdbj.out<-matrix(0,n.out,d*length(p))
ES.full.bdbj.out
dim(ES.full.bdbj.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.bdbj.out[,(i-1)*d+j]<-ES(p[i], out.full.bdbj$phit[,j], out.full.bdbj$xit[,j],  out.full.bdbj$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.bdbj.out.0.95 <-ES.full.bdbj.out[,1:d]
ES.full.bdbj.out.0.975<-ES.full.bdbj.out[,(d+1):(2*d)]
ES.full.bdbj.out.0.99 <-ES.full.bdbj.out[,(2*d+1):(3*d)]
ES.full.bdbj.out.0.999<-ES.full.bdbj.out[,(3*d+1):(4*d)]
plot.ts(ES.full.bdbj.out.0.95[100:dim(ES.full.bdbj.out.0.95)[1],1:d])
plot.ts(ES.full.bdbj.out.0.975[100:dim(ES.full.bdbj.out.0.975)[1],1:d])
plot.ts(ES.full.bdbj.out.0.99[100:dim(ES.full.bdbj.out.0.99)[1],1:d])
plot.ts(ES.full.bdbj.out.0.999[100:dim(ES.full.bdbj.out.0.999)[1],1:d])


## Accuracy test full bdbj
VAR.test.full.bdbj<-matrix(0,n.out,(length(p))*d)
ES.test.full.bdbj<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.bdbj[, (j-1)*length(p)+i]<-VAR.full.bdbj.out[,(i-1)*d+j]
    ES.test.full.bdbj[, (j-1)*length(p)+i]<-ES.full.bdbj.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.bdbj[2:dim(VAR.test.full.bdbj)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.bdbj[2:dim(ES.test.full.bdbj)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.bdbj[2:dim(ES.test.full.bdbj)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full bdbj
full.bdbj.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.bdbj.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.bdbj.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.bdbj.out[,(i-1)*d+j],p[i])
    full.bdbj.out.nol2[,(i-1)*d+j][is.na(full.bdbj.out.nol2[,(i-1)*d+j])]<-0
    full.bdbj.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.bdbj.out[,(i-1)*d+j],ES.full.bdbj.out[,(i-1)*d+j],p[i])
    full.bdbj.out.nol6[,(i-1)*d+j][is.na(full.bdbj.out.nol6[,(i-1)*d+j])]<-0
    full.bdbj.out.nol6[,(i-1)*d+j][is.infinite(full.bdbj.out.nol6[,(i-1)*d+j])]<-0
  }
}

#---------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic WEIBULL model 

#VaR full wd
VAR.full.wd.out<-matrix(0,n.out,d*length(p))
VAR.full.wd.out
dim(VAR.full.wd.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.wd.out[,(i-1)*d+j]<-VaR(p[i], out.full.wd$phit[,j], out.full.wd$xit[,j],  out.full.wd$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.wd.out.0.95 <-VAR.full.wd.out[,1:d]
VAR.full.wd.out.0.975<-VAR.full.wd.out[,(d+1):(2*d)]
VAR.full.wd.out.0.99 <-VAR.full.wd.out[,(2*d+1):(3*d)]
VAR.full.wd.out.0.999<-VAR.full.wd.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.wd.out.0.95[100:dim(VAR.full.wd.out.0.95)[1],1:d])
plot.ts(VAR.full.wd.out.0.975[100:dim(VAR.full.wd.out.0.975)[1],1:d])
plot.ts(VAR.full.wd.out.0.99[100:dim(VAR.full.wd.out.0.99)[1],1:d])
plot.ts(VAR.full.wd.out.0.999[100:dim(VAR.full.wd.out.0.999)[1],1:d])

#ES full wd
ES.full.wd.out<-matrix(0,n.out,d*length(p))
ES.full.wd.out
dim(ES.full.wd.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.wd.out[,(i-1)*d+j]<-ES(p[i], out.full.wd$phit[,j], out.full.wd$xit[,j],  out.full.wd$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.wd.out.0.95 <-ES.full.wd.out[,1:d]
ES.full.wd.out.0.975<-ES.full.wd.out[,(d+1):(2*d)]
ES.full.wd.out.0.99 <-ES.full.wd.out[,(2*d+1):(3*d)]
ES.full.wd.out.0.999<-ES.full.wd.out[,(3*d+1):(4*d)]
plot.ts(ES.full.wd.out.0.95[100:dim(ES.full.wd.out.0.95)[1],1:d])
plot.ts(ES.full.wd.out.0.975[100:dim(ES.full.wd.out.0.975)[1],1:d])
plot.ts(ES.full.wd.out.0.99[100:dim(ES.full.wd.out.0.99)[1],1:d])
plot.ts(ES.full.wd.out.0.999[100:dim(ES.full.wd.out.0.999)[1],1:d])

## Accuracy test full wd
VAR.test.full.wd<-matrix(0,n.out,(length(p))*d)
ES.test.full.wd<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.wd[, (j-1)*length(p)+i]<-VAR.full.wd.out[,(i-1)*d+j]
    ES.test.full.wd[, (j-1)*length(p)+i]<-ES.full.wd.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.wd[2:dim(VAR.test.full.wd)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.wd[2:dim(ES.test.full.wd)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.wd[2:dim(ES.test.full.wd)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full wd
full.wd.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.wd.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.wd.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.wd.out[,(i-1)*d+j],p[i])
    full.wd.out.nol2[,(i-1)*d+j][is.na(full.wd.out.nol2[,(i-1)*d+j])]<-0
    full.wd.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.wd.out[,(i-1)*d+j],ES.full.wd.out[,(i-1)*d+j],p[i])
    full.wd.out.nol6[,(i-1)*d+j][is.na(full.wd.out.nol6[,(i-1)*d+j])]<-0
    full.wd.out.nol6[,(i-1)*d+j][is.infinite(full.wd.out.nol6[,(i-1)*d+j])]<-0
  }
}


#-------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic WEIBULL model with standard realized measure

#VaR full wdr
VAR.full.wdr.out<-matrix(0,n.out,d*length(p))
VAR.full.wdr.out
dim(VAR.full.wdr.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.wdr.out[,(i-1)*d+j]<-VaR(p[i], out.full.wdr$phit[,j], out.full.wdr$xit[,j],  out.full.wdr$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.wdr.out.0.95 <-VAR.full.wdr.out[,1:d]
VAR.full.wdr.out.0.975<-VAR.full.wdr.out[,(d+1):(2*d)]
VAR.full.wdr.out.0.99 <-VAR.full.wdr.out[,(2*d+1):(3*d)]
VAR.full.wdr.out.0.999<-VAR.full.wdr.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.wdr.out.0.95[100:dim(VAR.full.wdr.out.0.95)[1],1:d])
plot.ts(VAR.full.wdr.out.0.975[100:dim(VAR.full.wdr.out.0.975)[1],1:d])
plot.ts(VAR.full.wdr.out.0.99[100:dim(VAR.full.wdr.out.0.99)[1],1:d])
plot.ts(VAR.full.wdr.out.0.999[100:dim(VAR.full.wdr.out.0.999)[1],1:d])

#ES full wdr
ES.full.wdr.out<-matrix(0,n.out,d*length(p))
ES.full.wdr.out
dim(ES.full.wdr.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.wdr.out[,(i-1)*d+j]<-ES(p[i], out.full.wdr$phit[,j], out.full.wdr$xit[,j],  out.full.wdr$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.wdr.out.0.95 <-ES.full.pdr.out[,1:d]
ES.full.wdr.out.0.975<-ES.full.pdr.out[,(d+1):(2*d)]
ES.full.wdr.out.0.99 <-ES.full.pdr.out[,(2*d+1):(3*d)]
ES.full.wdr.out.0.999<-ES.full.pdr.out[,(3*d+1):(4*d)]
plot.ts(ES.full.wdr.out.0.95[100:dim(ES.full.wdr.out.0.95)[1],1:d])
plot.ts(ES.full.wdr.out.0.975[100:dim(ES.full.wdr.out.0.975)[1],1:d])
plot.ts(ES.full.wdr.out.0.99[100:dim(ES.full.wdr.out.0.99)[1],1:d])
plot.ts(ES.full.wdr.out.0.999[100:dim(ES.full.wdr.out.0.999)[1],1:d])

## Accuracy test full wdr
VAR.test.full.wdr<-matrix(0,n.out,(length(p))*d)
ES.test.full.wdr<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.wdr[, (j-1)*length(p)+i]<-VAR.full.wdr.out[,(i-1)*d+j]
    ES.test.full.wdr[, (j-1)*length(p)+i]<-ES.full.wdr.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.wdr[2:dim(VAR.test.full.wdr)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.wdr[2:dim(ES.test.full.wdr)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.wdr[2:dim(ES.test.full.wdr)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full wdr
full.wdr.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.wdr.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.wdr.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.wdr.out[,(i-1)*d+j],p[i])
    full.wdr.out.nol2[,(i-1)*d+j][is.na(full.wdr.out.nol2[,(i-1)*d+j])]<-0
    full.wdr.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.wdr.out[,(i-1)*d+j],ES.full.wdr.out[,(i-1)*d+j],p[i])
    full.wdr.out.nol6[,(i-1)*d+j][is.na(full.wdr.out.nol6[,(i-1)*d+j])]<-0
    full.wdr.out.nol6[,(i-1)*d+j][is.infinite(full.wdr.out.nol6[,(i-1)*d+j])]<-0
  }
}

#---------------------------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic WEIBULL model with negative RSV

#VaR full wdrs
VAR.full.wdrs.out<-matrix(0,n.out,d*length(p))
VAR.full.wdrs.out
dim(VAR.full.wdrs.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.wdrs.out[,(i-1)*d+j]<-VaR(p[i], out.full.wdrs$phit[,j], out.full.wdrs$xit[,j],  out.full.wdrs$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.wdrs.out.0.95 <-VAR.full.wdrs.out[,1:d]
VAR.full.wdrs.out.0.975<-VAR.full.wdrs.out[,(d+1):(2*d)]
VAR.full.wdrs.out.0.99 <-VAR.full.wdrs.out[,(2*d+1):(3*d)]
VAR.full.wdrs.out.0.999<-VAR.full.wdrs.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.wdrs.out.0.95[100:dim(VAR.full.wdrs.out.0.95)[1],1:d])
plot.ts(VAR.full.wdrs.out.0.975[100:dim(VAR.full.wdrs.out.0.975)[1],1:d])
plot.ts(VAR.full.wdrs.out.0.99[100:dim(VAR.full.wdrs.out.0.99)[1],1:d])
plot.ts(VAR.full.wdrs.out.0.999[100:dim(VAR.full.wdrs.out.0.999)[1],1:d])

#ES full wdrs
ES.full.wdrs.out<-matrix(0,n.out,d*length(p))
ES.full.wdrs.out
dim(ES.full.wdrs.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.wdrs.out[,(i-1)*d+j]<-ES(p[i], out.full.wdrs$phit[,j], out.full.wdrs$xit[,j],  out.full.wdrs$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.wdrs.out.0.95 <-ES.full.wdrs.out[,1:d]
ES.full.wdrs.out.0.975<-ES.full.wdrs.out[,(d+1):(2*d)]
ES.full.wdrs.out.0.99 <-ES.full.wdrs.out[,(2*d+1):(3*d)]
ES.full.wdrs.out.0.999<-ES.full.wdrs.out[,(3*d+1):(4*d)]
plot.ts(ES.full.wdrs.out.0.95[100:dim(ES.full.wdrs.out.0.95)[1],1:d])
plot.ts(ES.full.wdrs.out.0.975[100:dim(ES.full.wdrs.out.0.975)[1],1:d])
plot.ts(ES.full.wdrs.out.0.99[100:dim(ES.full.wdrs.out.0.99)[1],1:d])
plot.ts(ES.full.wdrs.out.0.999[100:dim(ES.full.wdrs.out.0.999)[1],1:d])

## Accuracy test full wdrs
VAR.test.full.wdrs<-matrix(0,n.out,(length(p))*d)
ES.test.full.wdrs<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.wdrs[, (j-1)*length(p)+i]<-VAR.full.wdrs.out[,(i-1)*d+j]
    ES.test.full.wdrs[, (j-1)*length(p)+i]<-ES.full.wdrs.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.wdrs[2:dim(VAR.test.full.wdrs)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.wdrs[2:dim(ES.test.full.wdrs)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.wdrs[2:dim(ES.test.full.wdrs)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full wdrs
full.wdrs.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.wdrs.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.wdrs.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.wdrs.out[,(i-1)*d+j],p[i])
    full.wdrs.out.nol2[,(i-1)*d+j][is.na(full.wdrs.out.nol2[,(i-1)*d+j])]<-0
    full.wdrs.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.wdrs.out[,(i-1)*d+j],ES.full.wdrs.out[,(i-1)*d+j],p[i])
    full.wdrs.out.nol6[,(i-1)*d+j][is.na(full.wdrs.out.nol6[,(i-1)*d+j])]<-0
    full.wdrs.out.nol6[,(i-1)*d+j][is.infinite(full.wdrs.out.nol6[,(i-1)*d+j])]<-0
  }
}

#---------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic WEIBULL model with positive and negative RSV 

#VaR full wdpn
VAR.full.wdpn.out<-matrix(0,n.out,d*length(p))
VAR.full.wdpn.out
dim(VAR.full.wdpn.out)

for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.wdpn.out[,(i-1)*d+j]<-VaR(p[i], out.full.wdpn$phit[,j], out.full.wdpn$xit[,j],  out.full.wdpn$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.wdpn.out.0.95 <-VAR.full.wdpn.out[,1:d]
VAR.full.wdpn.out.0.975<-VAR.full.wdpn.out[,(d+1):(2*d)]
VAR.full.wdpn.out.0.99 <-VAR.full.wdpn.out[,(2*d+1):(3*d)]
VAR.full.wdpn.out.0.999<-VAR.full.wdpn.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.wdpn.out.0.95[100:dim(VAR.full.wdpn.out.0.95)[1],1:d])
plot.ts(VAR.full.wdpn.out.0.975[100:dim(VAR.full.wdpn.out.0.975)[1],1:d])
plot.ts(VAR.full.wdpn.out.0.99[100:dim(VAR.full.wdpn.out.0.99)[1],1:d])
plot.ts(VAR.full.wdpn.out.0.999[100:dim(VAR.full.wdpn.out.0.999)[1],1:d])

#ES full wdpn
ES.full.wdpn.out<-matrix(0,n.out,d*length(p))
ES.full.wdpn.out
dim(ES.full.wdpn.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.wdpn.out[,(i-1)*d+j]<-ES(p[i], out.full.wdpn$phit[,j], out.full.wdpn$xit[,j],  out.full.wdpn$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.wdpn.out.0.95 <-ES.full.wdpn.out[,1:d]
ES.full.wdpn.out.0.975<-ES.full.wdpn.out[,(d+1):(2*d)]
ES.full.wdpn.out.0.99 <-ES.full.wdpn.out[,(2*d+1):(3*d)]
ES.full.wdpn.out.0.999<-ES.full.wdpn.out[,(3*d+1):(4*d)]
plot.ts(ES.full.wdpn.out.0.95[100:dim(ES.full.wdpn.out.0.95)[1],1:d])
plot.ts(ES.full.wdpn.out.0.975[100:dim(ES.full.wdpn.out.0.975)[1],1:d])
plot.ts(ES.full.wdpn.out.0.99[100:dim(ES.full.wdpn.out.0.99)[1],1:d])
plot.ts(ES.full.wdpn.out.0.999[100:dim(ES.full.wdpn.out.0.999)[1],1:d])


## Accuracy test full wdpn
VAR.test.full.wdpn<-matrix(0,n.out,(length(p))*d)
ES.test.full.wdpn<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.wdpn[, (j-1)*length(p)+i]<-VAR.full.wdpn.out[,(i-1)*d+j]
    ES.test.full.wdpn[, (j-1)*length(p)+i]<-ES.full.wdpn.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.wdpn[2:dim(VAR.test.full.wdpn)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.wdpn[2:dim(ES.test.full.wdpn)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.wdpn[2:dim(ES.test.full.wdpn)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full wdpn
full.wdpn.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.wdpn.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.wdpn.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.wdpn.out[,(i-1)*d+j],p[i])
    full.wdpn.out.nol2[,(i-1)*d+j][is.na(full.wdpn.out.nol2[,(i-1)*d+j])]<-0
    full.wdpn.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.wdpn.out[,(i-1)*d+j],ES.full.wdpn.out[,(i-1)*d+j],p[i])
    full.wdpn.out.nol6[,(i-1)*d+j][is.na(full.wdpn.out.nol6[,(i-1)*d+j])]<-0
    full.wdpn.out.nol6[,(i-1)*d+j][is.infinite(full.wdpn.out.nol6[,(i-1)*d+j])]<-0
  }
}

#---------------------------------------------------------------------------------
# Full (non-diagonal) and dynamic WEIBULL model with bipower and jumps
#VaR full wdbj
VAR.full.wdbj.out<-matrix(0,n.out,d*length(p))
VAR.full.wdbj.out
dim(VAR.full.wdbj.out)
for(i in 1:length(p)){
  for(j in 1:d){
    VAR.full.wdbj.out[,(i-1)*d+j]<-VaR(p[i], out.full.wdbj$phit[,j], out.full.wdbj$xit[,j],  out.full.wdbj$dt[j] ,umbral.out[j], model="MPOT")
  }
}
VAR.full.wdbj.out.0.95 <-VAR.full.wdbj.out[,1:d]
VAR.full.wdbj.out.0.975<-VAR.full.wdbj.out[,(d+1):(2*d)]
VAR.full.wdbj.out.0.99 <-VAR.full.wdbj.out[,(2*d+1):(3*d)]
VAR.full.wdbj.out.0.999<-VAR.full.wdbj.out[,(3*d+1):(4*d)]
plot.ts(VAR.full.wdbj.out.0.95[100:dim(VAR.full.wdbj.out.0.95)[1],1:d])
plot.ts(VAR.full.wdbj.out.0.975[100:dim(VAR.full.wdbj.out.0.975)[1],1:d])
plot.ts(VAR.full.wdbj.out.0.99[100:dim(VAR.full.wdbj.out.0.99)[1],1:d])
plot.ts(VAR.full.wdbj.out.0.999[100:dim(VAR.full.wdbj.out.0.999)[1],1:d])

#ES full wdbj
ES.full.wdbj.out<-matrix(0,n.out,d*length(p))
ES.full.wdbj.out
dim(ES.full.wdbj.out)
for(i in 1:length(p)){
  for(j in 1:d){
    ES.full.wdbj.out[,(i-1)*d+j]<-ES(p[i], out.full.wdbj$phit[,j], out.full.wdbj$xit[,j],  out.full.wdbj$dt[j],  umbral.out[j], model="MPOT")
  }
}
ES.full.wdbj.out.0.95 <-ES.full.wdbj.out[,1:d]
ES.full.wdbj.out.0.975<-ES.full.wdbj.out[,(d+1):(2*d)]
ES.full.wdbj.out.0.99 <-ES.full.wdbj.out[,(2*d+1):(3*d)]
ES.full.wdbj.out.0.999<-ES.full.wdbj.out[,(3*d+1):(4*d)]
plot.ts(ES.full.wdbj.out.0.95[100:dim(ES.full.wdbj.out.0.95)[1],1:d])
plot.ts(ES.full.wdbj.out.0.975[100:dim(ES.full.wdbj.out.0.975)[1],1:d])
plot.ts(ES.full.wdbj.out.0.99[100:dim(ES.full.wdbj.out.0.99)[1],1:d])
plot.ts(ES.full.wdbj.out.0.999[100:dim(ES.full.wdbj.out.0.999)[1],1:d])


## Accuracy test full wdbj
VAR.test.full.wdbj<-matrix(0,n.out,(length(p))*d)
ES.test.full.wdbj<-matrix(0,n.out,(length(p))*d)
for(j in 1:d){
  for(i in 1:length(p)){
    VAR.test.full.wdbj[, (j-1)*length(p)+i]<-VAR.full.wdbj.out[,(i-1)*d+j]
    ES.test.full.wdbj[, (j-1)*length(p)+i]<-ES.full.wdbj.out[,(i-1)*d+j]
  }
  VAR.test<-backtestVaR(VAR.test.full.wdbj[2:dim(VAR.test.full.wdbj)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
  ES.test<-backtestES(VAR.test.full.wdbj[2:dim(ES.test.full.wdbj)[1],(length(p)*(j-1)+1):(length(p)*j)], ES.test.full.wdbj[2:dim(ES.test.full.wdbj)[1],(length(p)*(j-1)+1):(length(p)*j)],yt.out[2:dim(yt.out)[1],j])
}

## loss functions for full wdbj
full.wdbj.out.nol2<-matrix(0,dim(yt.out)[1],d*length(p))
full.wdbj.out.nol6<-matrix(0,dim(yt.out)[1],d*length(p))
for(i in 1:length(p)){
  for(j in 1:d){
    full.wdbj.out.nol2[,(i-1)*d+j] <-nol2(yt.out[,j],VAR.full.wdbj.out[,(i-1)*d+j],p[i])
    full.wdbj.out.nol2[,(i-1)*d+j][is.na(full.wdbj.out.nol2[,(i-1)*d+j])]<-0
    full.wdbj.out.nol6[,(i-1)*d+j] <-nol6(yt.out[,j],VAR.full.wdbj.out[,(i-1)*d+j],ES.full.wdbj.out[,(i-1)*d+j],p[i])
    full.wdbj.out.nol6[,(i-1)*d+j][is.na(full.wdbj.out.nol6[,(i-1)*d+j])]<-0
    full.wdbj.out.nol6[,(i-1)*d+j][is.infinite(full.wdbj.out.nol6[,(i-1)*d+j])]<-0
  }
}

###############################################################################################
# Model Confidence Set
################################################################################################

#sink("./MCS_resultados01.csv", append = T, split= T)
#sink("./MCS_resultados01.txt", append = T, split= T)


## Zero-degree homogeneous functions nol2 for VaR, and nol6 for VaR and ES



fnol2.M<-rep(0,length(p)*d)
for(i in 1:length(p)){
  for(j in 1:d){
    fnol2.M[(i-1)*d+j]<- MCSprocedure(Loss = cbind(

      full.pd.out.nol2[,(i-1)*d+j],
      full.bd.out.nol2[,(i-1)*d+j],
      full.wd.out.nol2[,(i-1)*d+j],
      
      full.pdr.out.nol2[,(i-1)*d+j],
      full.bdr.out.nol2[,(i-1)*d+j],
      full.wdr.out.nol2[,(i-1)*d+j],
      
      full.pdrs.out.nol2[,(i-1)*d+j],
      full.bdrs.out.nol2[,(i-1)*d+j],
      full.wdrs.out.nol2[,(i-1)*d+j],
      
      full.pdpn.out.nol2[,(i-1)*d+j],
      full.bdpn.out.nol2[,(i-1)*d+j],
      full.wdpn.out.nol2[,(i-1)*d+j],
      
      full.pdbj.out.nol2[,(i-1)*d+j],
      full.bdbj.out.nol2[,(i-1)*d+j],
      full.wdbj.out.nol2[,(i-1)*d+j]
    ),
    alpha = 0.95, B = 5000, statistic = "Tmax")@Info$model.names
  }
}

fnol2.M<-matrix(fnol2.M, d,length(p))
colnames(fnol2.M)<-c("0.95","0.975","0.99","0.999")
rownames(fnol2.M)<-c("D1","D2","D3","D4","D5")

fnol2.M


fnol6.M<-rep(0,length(p)*d)
for(i in 1:length(p)){
  for(j in 1:d){
    fnol6.M[(i-1)*d+j]<- MCSprocedure(Loss = cbind(
    
      full.pd.out.nol6[,(i-1)*d+j],
      full.bd.out.nol6[,(i-1)*d+j],
      full.wd.out.nol6[,(i-1)*d+j],
      
      full.pdr.out.nol6[,(i-1)*d+j],
      full.bdr.out.nol6[,(i-1)*d+j],
      full.wdr.out.nol6[,(i-1)*d+j],
      
     full.pdrs.out.nol6[,(i-1)*d+j],
     full.bdrs.out.nol6[,(i-1)*d+j],
     full.wdrs.out.nol6[,(i-1)*d+j],
      
      full.pdpn.out.nol6[,(i-1)*d+j],
      full.bdpn.out.nol6[,(i-1)*d+j],
      full.wdpn.out.nol6[,(i-1)*d+j],
      
      full.pdbj.out.nol6[,(i-1)*d+j],
      full.bdbj.out.nol6[,(i-1)*d+j],
      full.wdbj.out.nol6[,(i-1)*d+j]
    ),
    alpha = 0.95, B = 5000, statistic = "Tmax")@Info$model.names
  }
}

fnol6.M<-matrix(fnol6.M, d,length(p))
colnames(fnol6.M)<-c("0.95","0.975","0.99","0.999")
rownames(fnol6.M)<-c("D1","D2","D3","D4","D5")

fnol6.M

############################################################################################################################
# Getting data for plots 
#############################################################################################################################

fecha<-read.csv("rets_tech_stocks.csv",header=TRUE,sep=",")
fecha<-fecha[4085:5333,1]
head(fecha)

## getting date
# write.csv(ret[4085:5333,1],file="C:\\Users\\Claudio Candia C\\Dropbox\\Claudio Candia\\Paper 3s.1\\VAR time X\\datum.csv")

d1<-data.frame(yt.out[,1],VAR.full.pd.out.0.99[,1], ES.full.pd.out.0.975[,1], rep("D1", length(yt.out[,1])))
colnames(d1)<-c("retornos", "var","es", "market")

d2<-data.frame(yt.out[,2],VAR.full.pd.out.0.99[,2], ES.full.pd.out.0.975[,2], rep("D2", length(yt.out[,2])))
colnames(d2)<-c("retornos", "var","es", "market")

d3<-data.frame(yt.out[,3],VAR.full.pd.out.0.99[,3], ES.full.pd.out.0.975[,3], rep("D3", length(yt.out[,3])))
colnames(d3)<-c("retornos", "var","es", "market")

d4<-data.frame(yt.out[,4],VAR.full.pd.out.0.99[,4], ES.full.pd.out.0.975[,4], rep("D4", length(yt.out[,4])))
colnames(d4)<-c("retornos", "var","es", "market")

d5<-data.frame(yt.out[,5],VAR.full.pd.out.0.99[,5], ES.full.pd.out.0.975[,5], rep("D5", length(yt.out[,5])))
colnames(d5)<-c("retornos", "var","es", "market")


RISK.TAIL<-rbind(d1,d2,d3,d4,d5)
head(RISK.TAIL)

## getting VaR (0.99) and ES (0.975)
#write.csv(RISK.TAIL,file="C:\\Users\\Claudio Candia C\\Dropbox\\Claudio Candia\\Paper 3s.1\\VAR time X\\RISK.TAIL.csv")



























