rm(list=ls())
setwd("C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II")
source("Functions.R")
library(foreach)
library(doParallel)
ns=2
cod_cov=c('ar','bd')
cod_OC=c('gaussian','t','gamma')
cod_rho=seq(0.1,0.9,0.1)
df.temp2 <- data.frame()
#---#
cl=makeCluster(2)
registerDoParallel(cl)
#---#
for(ind.cov in 1:length(cod_cov)){
 for(ind.OC in 1:length(cod_OC)){
  for(ind.rho in 1:length(cod_rho)){
    cov.type = cod_cov[ind.cov]
    OC=cod_OC[ind.OC]
    rho=cod_rho[ind.rho]
    ###The setting of Figure 1###
    n1=400;n2=100;p=600;
    op=0.1;del=0.3
    ###
    det1=-del
    det2=del
    alpha=c(0.2)
    Loop=foreach(l=1:ns,.combine=list,.multicombine = TRUE,.maxcombine = ns,
                 .packages= c("MASS","glmnet",'ks','pracma',
                              'foreach','doParallel','mvtnorm'))%dopar%{
                                ####
                                set.seed(l+1234)
                                if(cov.type=='ar'){Sig1<-toeplitz(rho^(0:(p-1)))}
                                if(cov.type=='bd'){
                                  Sig1<-Gen_Sigma(p,rho,cov.type=cov.type) 
                                }
                                ##
                                op1=0.5;sig12=0.05;sig2=sig12^2
                                dat<-Gen_highmean(n1=n1,n2=n2,p=p,op=op,det1=det1,det2=det2,Sig1,OC=OC,OCmu='II')
                                X1<-dat$X1
                                Sig1<-dat$Sig1
                                mu1<-dat$mu1
                                X2<-dat$X2
                                Sig2<-dat$Sig2
                                mu2<-dat$mu2
                                Ind<-dat$Ind
                                THETA<-dat$THETA
                                #
                                n1<-nrow(X1);n2<-nrow(X2)
                                Sig<-1/n1*Sig1+1/n2*Sig2
                                Omega0<-solve(Sig)
                                Omega0<-(Omega0+t(Omega0))/2
                                #                              
                                SDA0=0
                                MOW=0
                                #
                                SDA=1
                                MDW=1
                                #
                                reject.ratio=1
                                #--------------#
                                RES2<-list()
                                Allres = list(RES2=RES2);k2=0;Method2=NULL
                                if(SDA0){
                                  k2=k2+1
                                  Det <- Func_SDA2(X1, X2,alpha=alpha, Omega=Omega0, model ='mean',stable = F)
                                  stat=Det$stat[[1]]
                                  RES2<-Func_TPR(THETA,stat,reject.ratio)
                                  Allres$RES2[[k2]]<-RES2
                                  Method2<-c(Method2,do.call(c,lapply(reject.ratio, function(x){paste0('SDA0-',x)})))
                                }
                                if(MOW){
                                  k2=k2+1
                                  res<-Func_MOW(X2,op,op1,det1,det2,sig2,alpha)
                                  stat<-res$stat
                                  RES2<-Func_TPR(THETA,-stat,reject.ratio)
                                  Allres$RES2[[k2]]=RES2
                                  Method2<-c(Method2,do.call(c,lapply(reject.ratio, function(x){paste0('MOW-',x)})))
                                }
                                if(SDA){
                                  k2=k2+1
                                  Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K=1)
                                  stat=Det$stat[[1]]
                                  RES2<-Func_TPR(THETA,stat,reject.ratio)
                                  Allres$RES2[[k2]]<-RES2
                                  Method2<-c(Method2,do.call(c,lapply(reject.ratio, function(x){paste0('SDA-',x)})))
                                }
                                if(MDW){
                                  k2=k2+1
                                  res<-Func_MDW(X2,lam=0.1,b_value = 5,alpha)
                                  stat<-res$stat
                                  RES2<-Func_TPR(THETA,-stat,reject.ratio)
                                  Allres$RES2[[k2]]=RES2
                                  Method2<-c(Method2,do.call(c,lapply(reject.ratio, function(x){paste0('MDW-',x)})))
                                }
                                ######corrected power######
                                ResSummary2 = lapply(1:length(Allres$RES2),function(x){Allres$RES2[[x]]})
                                ResSummary2 = do.call(rbind,ResSummary2)
                                rownames(ResSummary2)=Method2
                          
                                Result2=round(ResSummary2,digits = 3)
                                return(RES=list(Result2=Result2))
                              }
    #----------------#
    output2<-as.matrix(do.call(rbind,lapply(1:ns,function(x){Loop[[x]]$Result2})))
    nn2=nrow(output2)
    ###
    MeanRes2<-do.call(rbind,lapply(1:(nn2/ns),function(x){
      #x=1
      dat<-as.data.frame(t(mean(output2[seq(x,nn2,nn2/ns),])))
      rownames(dat)=rownames(output2)[x]
      return(dat)
    }))
    names(MeanRes2)='TPR'
    paras=c(OC,p,cov.type,rho,op,n1,n2,det1,det2)
    paras_name=c('OC','p','cov.type','rho','op','n1','n2','det1','det2')
    ##
    mat_set2=matrix(0,nn2/ns,length(paras))
    for(i in 1:length(paras)){
      mat_set2[,i]=paras[i]
    }
    dfmat2=cbind(MeanRes2,mat_set2)
    colnames(dfmat2)=c('TPR',paras_name)                  
    ##
    df.temp2=rbind(df.temp2,dfmat2)
  }
 }
}
stopCluster(cl)
print(df.temp2)
#write.csv(df.temp1,file = 'Res_Figure1.csv')





