#This file outputs results in Figure 1#
rm(list=ls())
setwd("C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II")
source("Functions.R")
alpha = 0.2;
#=============================#
#The diagnosis of mean shift# 
#=============================#
set.seed(1234)
OC = 'gaussian'
cov.type = 'ar'
del = 0.3
n1 = 400;
n2 = 100;
p = 200;
op = 0.1;
rho = 0.8
det1=-del
det2=del
###
Sig1<-toeplitz(rho^(0:(p-1)))#
dat<-Gen_highmean(n1=n1,n2=n2,p=p,op=op,det1=det1,det2=det2,Sig1,OC=OC,OCmu='II')
X1 <- dat$X1
X2 <- dat$X2
THETA <- dat$THETA
#--------------#
Det <- Func_SDA2(X1, X2,alpha, Omega=NULL, model ='mean',stable = F)                                      
tru_mean <- which(THETA!=0)
W_mean=Det$stat[[1]]
plot(W_mean)
###The diagnosis of the general linear profile##
model='glm';
cov.type='ar';
p=200;
rho=0;
op = 20;
m=1500
n1 <- 200
n2 <- 100
del=0.05
det1=det2=del
# Sig<-Gen_Sigma(p,rho)
# X1<-rmvnorm(m,rep(0,p),Sig)
# system.time(dat<-Gen_linearprofile(X1,n1,n2,m,p,rho=0,op,det1,det2,model))
# Y1<-dat$Y1
# Y2<-dat$Y2
# Ind<-dat$Ind
# system.time(Ome<-Func_omega(X1,X2=NULL,Y1,Y2,model,maxit=maxit))
# Beta1<-Ome$Beta1
# Beta2<-Ome$Beta2
# S1<-Ome$S1
# S2<-Ome$S2
#Det <- Func_SDA2(X1, X2=NULL,Y1,Y2,Beta1=Beta1,Beta2=Beta2,S1=S1,S2=S2,alpha=alpha, Omega=NULL, model =model,stable=F)
#tru <- Ind
#W=Det$stat[[1]]
###
#ind <- rep(0,p)
#ind[tru]=1
#dat_W <- cbind(ind,W)
#write.csv(dat_W,file='glm_W.csv')
dat_W <- read.csv('glm_W.csv')[,-1]
W_glm <- dat_W$W;
plot(W_glm)
tru_glm <- which(dat_W$ind!=0)
