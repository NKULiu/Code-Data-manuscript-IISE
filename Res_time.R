###This file records the time of MDW,LEB and dSDA#
rm(list=ls())
setwd('C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II')
source('Functions.R')
library(foreach)
library(doParallel)
ns = 100
cod_OC = c('gaussian')
cod_cov = c('ar')
cod_op = c(0.1)
cod_delta = 0.3
cod_n2 = 100
cod_p = c(100,200,400,600,800)
df.temp1 <- data.frame()
#---#
tic = proc.time()
cl = makeCluster(8)
registerDoParallel(cl)
#---#
for(ind.n2 in 1:length(cod_n2)){
  for(ind.cov in 1:length(cod_cov)){
    for(ind.OC in 1:length(cod_OC)){
      for(ind.op in 1:length(cod_op)){
        for(ind.delta in 1:length(cod_delta)){
          for(ind.p in 1:length(cod_p)){
          ###
          n2 <- cod_n2[ind.n2] 
          p <- cod_p[ind.p]
          delta <- cod_delta[ind.delta]
          OC=cod_OC[ind.OC]
          cov.type=cod_cov[ind.cov]
          op=cod_op[ind.op]
          n1=400;n2=n2;del = delta;
          if(cov.type == 'ar'){rho = 0.8};
          if(cov.type == 'bd'){rho = 0.6};
          det1 = -delta
          det2 = delta
          alpha=0.2
          Loop=foreach(l=1:ns,.combine=list,.multicombine = TRUE,.maxcombine = ns,
                       .packages= c("MASS","glmnet",'pracma','ks',
                                    'foreach','doParallel','mvtnorm'))%dopar%{
                                      set.seed(l+1234)
                                      print(l)
                                      if(cov.type=='ar'){
                                        Sig1<-toeplitz(rho^(0:(p-1)))
                                      }
                                      if(cov.type=='bd'){
                                        Sig1<-Gen_Sigma(p,rho,cov.type=cov.type)
                                      }
                                      ## 
                                      op1 = 0.5;sig12 = 0.05; sig2 = sig12^2;
                                      dat<-Gen_highmean(n1=n1,n2=n2,p=p,op=op,det1=det1,det2=det2,Sig1,OC=OC,OCmu='II')
                                      X1<-dat$X1
                                      Sig1<-dat$Sig1
                                      mu1<-dat$mu1
                                      X2<-dat$X2
                                      Sig2<-dat$Sig2
                                      mu2<-dat$mu2
                                      Ind<-dat$Ind
                                      THETA<-dat$THETA
                                      #estimated omega via full data# 
                                      Time_est<-system.time(Omega<-Func_omega(X1,X2,Y1,Y2,model='mean'))[3]
                                      #--------------#
                                      LEB = 1
                                      SDA = 1
                                      MDW = 1  
                                      ###
                                      RES1<-list()
                                      Allres = list(RES1=RES1);k1=0;Method1=NULL;
                                      if(SDA){
                                        k1=k1+1
                                        Tim_SDA<-system.time(Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K=1))[3]
                                        Allres$RES1[[k1]]<-Tim_SDA
                                        Method1<-c(Method1,'SDA')
                                      }
                                      if(MDW){
                                        k1=k1+1
                                        Tim_MDW<-system.time(res<-Func_MDW(X2,lam=0.1,b_value = 5,alpha))[3]
                                        Allres$RES1[[k1]]=Tim_MDW
                                        Method1<-c(Method1,'MDW')
                                      }
                                      if(LEB){
                                        k1=k1+1
                                        Tim_LEB<-system.time(Det <- Func_LEB2(X1, X2, Omega=Omega,model='mean'))[3]+Time_est
                                        Allres$RES1[[k1]]<-Tim_LEB
                                        Method1=c(Method1,'LEB');
                                      }
                                      ###time###
                                      ResSummary1 = lapply(1:length(Allres$RES1),function(x){Allres$RES1[[x]]})
                                      ResSummary1 = do.call(rbind,ResSummary1);
                                      rownames(ResSummary1) = Method1;
                                      colnames(ResSummary1) = 'Time';
                                      Result1=round(ResSummary1,digits = 3)
                                      return(RES=list(Result1=Result1))
                                    }
          #----------------#
          output1<-do.call(rbind,lapply(1:ns,function(x){Loop[[x]]$Result1}))
          nn1=nrow(output1)
          MeanRes1<-do.call(rbind,lapply(1:(nn1/ns),function(x){
            #x=1
            dat<-as.data.frame(mean(output1[seq(x,nn1,nn1/ns),]))
            rownames(dat)=rownames(output1)[x]
            return(dat)
          }))
          MeanRes1 <- round(MeanRes1,3)
          paras=c(p)
          paras_name='p'
          ##
          mat_set1=matrix(0,nn1/ns,length(paras))
          for(i in 1:length(paras)){
            mat_set1[,i]=paras[i]
          }
          dfmat1=cbind(MeanRes1,mat_set1)
          colnames(dfmat1)=c('meanTime',paras_name)
          ##
          filename=c()
          for (i in 1:length(paras)){
            temp=paste0(paras_name[i],'-',paras[i])
            filename=paste0(filename,temp,sep='-')
          }
          print(paste0(filename,Sys.time(),sep='-'))
          df.temp1=rbind(df.temp1,dfmat1)
         }
        }
      }
    }
   }
  }
stopCluster(cl)
proc.time() - tic
###Table of time###
p100 <- as.data.frame(matrix(df.temp1[df.temp1$p==100,1],ncol=1))
p200 <- as.data.frame(matrix(df.temp1[df.temp1$p==200,1],ncol=1))
p400 <- as.data.frame(matrix(df.temp1[df.temp1$p==400,1],ncol=1))
p600 <- as.data.frame(matrix(df.temp1[df.temp1$p==600,1],ncol=1))
p800 <- as.data.frame(matrix(df.temp1[df.temp1$p==800,1],ncol=1))
bb <- cbind(p100,p200,p400,p600,p800)
bb <- round(bb,1)
row.names(bb) = c('dSDA','MDW','LEB')
colnames(bb) = c(100,200,400,600,800)
library(xtable)
print(xtable::xtable(bb))
