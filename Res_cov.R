rm(list=ls())
library(parallel)
library(foreach)
library(doParallel)
setwd('C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II')
source('Functions.R')
###-------------------------------------###
##We consider only covariance change##   
###-------------------------------------###
num_cl <- detectCores(all=T)  
cl <- makeCluster(num_cl)
registerDoParallel(cl)
rep <- 100
dat_g <- foreach(ii=1:rep,
                .combine = 'rbind',
                .maxcombine = rep,
                .packages = c('mvtnorm','glmnet','MASS',
                              'matrixcalc','Matrix')
  )%dopar%{ 
    set.seed(ii+1234)
    n1 = 1000; 
    n2 = 500; 
    p = 10; 
    det1 = 0.1;
    det2 = 0.3
    alpha = c(0.05,0.1,0.2);
    Sig1 = Gen_Sigma(p,rho=0.5,'ar');
    Sig2=Sig1
    mu1=mu2=rep(0,p)
    diffnumber = 10;#oc elements of Sig1#
    M = upper.tri(Sig1,diag = F)
    ind <- which(M!=0)
    odind<-sample(ind,diffnumber)
    MM=matrix(0,p,p)
    MM[odind]<-runif(length(odind),det1,det2)*(2*rbinom(length(odind),1,1/2)-1)
    MM=MM+t(MM)
    Sig2 <- Sig1+MM
    #
    eigen1<-eigen(Sig1)$values
    eigen2<-eigen(Sig2)$values  
    delta<-abs(min(eigen1,eigen2)) + 0.01
    #
    Sig1 <- Sig1 + delta*diag(p)
    Sig2 <- Sig2 + delta*diag(p)
    #
    dat1<-mvrnorm(n1,mu1,Sig1)
    dat2<-rmvnorm(n2,mu2,Sig2)
    #lower.tri locations#
    low.loc<-lower.tri(Sig1,diag = T)
    lower_S1 <- Sig1[low.loc]
    lower_S2 <- Sig2[low.loc]
    beta1 <- c(lower_S1);
    beta2 <- c(lower_S2);  
    #True ind#
    Ind <- which((beta2-beta1)!=0)
    #And the consistent beta estimation# 
    S1_est <- cov(dat1);
    S2_est <- cov(dat2);
    lower_S1_est <- S1_est[low.loc]
    lower_S2_est <- S2_est[low.loc]
    ##
    beta1 <- c(lower_S1_est);
    beta2 <- c(lower_S2_est);
    beta <- beta1 - beta2;
    #estimate the covariance of correlation# 
    mu1_est <- apply(dat1,2,mean);
    mu2_est <- apply(dat2,2,mean);
    L <- elimination.matrix(p)
    bar1 <- mu1_est;
    bar2 <- mu2_est;
    ###
    sum_i<-0
    for(i in 1:n1){
      xi<-dat1[i,]
      tmp1<-(xi-bar1)%*%t(xi-bar1)
      tmp<-kronecker(tmp1,tmp1)
      sum_i<-sum_i+tmp
    }
    M1 <- sum_i/n1
    ###
    sum_j<-0
    for(j in 1:n2){
      xj<-dat2[j,]
      tmp1<-(xj-bar2)%*%t(xj-bar2)
      tmp<-kronecker(tmp1,tmp1)
      sum_j<-sum_j+tmp
    }
    M2<-sum_j/n2
    #
    V1 <- M1 - vec(S1_est)%*%t(vec(S1_est)); 
    V2 <- M2 - vec(S2_est)%*%t(vec(S2_est));
    ###covariance of sample covariance matrix#
    S1 <- L%*%V1%*%t(L)
    S2 <- L%*%V2%*%t(L)
    #
    S <- 1/n1*S1+1/n2*S2
    Omega <- solve(S)
    Gamma = Sqrt(Omega)
    #
    Det_SDA <- Func_SDA2(dat1,dat2,alpha=alpha,Omega=Omega, model ='cov',K=10)$out
    measure_SDA <- do.call(rbind,lapply(1:length(alpha),function(x){c(Func_Measure(Det_SDA[[x]],Ind,length(beta)),alpha[x])}))
    #
    Det_LEB <- Func_LEB2(dat1, dat2, beta=beta,Omega=Omega)
    measure_LEB <- do.call(rbind,lapply(1:length(alpha),function(x){c(Func_Measure(Det_LEB,Ind,length(beta)),alpha[x])}))
    return(as.data.frame(rbind(measure_LEB,measure_SDA)))
    gc()
  }
  stopCluster(cl)
  ###
  names(dat_g) <- c('MS','FDP','TDP','alpha')
  dat_g$method <- rep(c(rep('LEB',3),rep('SDA',3)))
  ###
  dat_LEB <- dat_g[dat_g$method=='LEB',]
  d_LEB <- apply(dat_LEB[,1:4],2,mean)
  dat_SDA <- dat_g[dat_g$method=='SDA',]
  alp0.05 <- dat_SDA[dat_SDA$alpha==0.05,]
  alp0.1 <- dat_SDA[dat_SDA$alpha==0.1,]
  alp0.2 <- dat_SDA[dat_SDA$alpha==0.2,]
  d_SDA1 <- apply(alp0.05[,1:4],2,mean)
  d_SDA2 <- apply(alp0.1[,1:4],2,mean)
  d_SDA3 <- apply(alp0.2[,1:4],2,mean)
  dat <- rbind(d_LEB,d_SDA1,d_SDA2,d_SDA3)
###
library(xtable)
xtable(dat)
