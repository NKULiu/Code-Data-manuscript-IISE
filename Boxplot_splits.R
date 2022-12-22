rm(list=ls())
setwd('C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II')
source('Functions.R')
setwd('C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II/Box')
library(foreach)
library(doParallel)
tic=proc.time()
ns = 100
#---------------#
num_cl <- detectCores()
cl=makeCluster(num_cl)
registerDoParallel(cl)
cod_OC = c('gaussian','t','gamma')
cod_cov = c('ar','bd')
cod_op = c(0.1)
df.temp <- data.frame()
for(ind.cov in 1:length(cod_cov)){
  for(ind.OC in 1:length(cod_OC)){
    for(ind.op in 1:length(cod_op)){
      ###
      OC=cod_OC[ind.OC]
      cov.type=cod_cov[ind.cov]
      op=cod_op[ind.op]
      n1 = 400;n2 = 100;p = 600; det1 = 0.3;det2 = 0.3
      ###
      if(cov.type%in%c('ar')){rho = 0.8};
      if(cov.type%in%c('bd')){rho = 0.6};
      alpha = c(0.05,0.1,0.2)
      Loop=foreach(l=1:ns,.combine=list,.multicombine = TRUE,.maxcombine = ns,
                   .packages= c("MASS","glmnet",'ks','pracma',
                                'foreach','doParallel','mvtnorm'))%dopar%{
                                  set.seed(l+1234)
                                  if(cov.type=='ar'){Sig1<-toeplitz(rho^(0:(p-1)))}
                                  if(cov.type%in%c('bd')){
                                    Sig1<-Gen_Sigma(p,rho,cov.type=cov.type) 
                                  }
                                  op1 = 0.5;sig12 = 0.05;sig2 = sig12^2;
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
                                  SDA_K1 = 1#estimated omega by the first split
                                  
                                  SDA_K5 = 1
                                  SDA_K10 = 1
                                  SDA_K15 = 1
                                  SDA_K20 = 1
                                  RES1<-list()
                                  Allres = list(RES1=RES1);k1=0;Method1=NULL;
                                  if(SDA_K1){
                                    k1=k1+1
                                    Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K=1)
                                    Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                    Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA-K1',x)})))
                                  }
                                  #
                                  
                                  #
                                  if(SDA_K5){
                                    k1=k1+1
                                    Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K=5)
                                    Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                    Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA-K5',x)})))
                                  }
                                  #
                                  if(SDA_K10){
                                    k1=k1+1
                                    Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K = 10)
                                    Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                    Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA-K10',x)})))
                                  }
                                  #
                                  if(SDA_K15){
                                    k1=k1+1
                                    Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K = 15)
                                    Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                    Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA-K15',x)})))
                                  }
                                  #
                                  if(SDA_K20){
                                    k1=k1+1
                                    Det <- Func_SDA2(X1, X2,alpha=alpha, model ='mean',K=20)
                                    Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                    Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA-K20',x)})))
                                  }
                                  #####FDP and TPP#####
                                  ResSummary1 = lapply(1:length(Allres$RES1),function(x){Allres$RES1[[x]]})
                                  ResSummary1 = do.call(rbind,ResSummary1);
                                  rownames(ResSummary1)=Method1;
                                  colnames(ResSummary1)=c('MS',"FDP","TPP")
                                  Result1=round(ResSummary1,digits = 3)
                                  return(RES=list(Result1=Result1))
                                }
      #----------------#
      output1<-do.call(rbind,lapply(1:ns,function(x){Loop[[x]]$Result1}))
      nn1=nrow(output1)
      MeanRes1<-do.call(rbind,lapply(1:(nn1/ns),function(x){
        #x=1
        dat<-as.data.frame(t(colMeans(output1[seq(x,nn1,nn1/ns),])))
        rownames(dat)=rownames(output1)[x]
        return(dat)
      }))
      ##
      sdRes1<-do.call(rbind,lapply(1:(nn1/ns),function(x){
        #x=1
        dat<-as.data.frame(t(apply(output1[seq(x,nn1,nn1/ns),],2,sd)))
        rownames(dat)=rownames(output1)[x]
        return(dat)
      }))
      ##
      MeanRes1 <- round(MeanRes1,3)
      sdRes1 <- round(sdRes1,3)
      Res1<-cbind(MeanRes1,sdRes1)
      ###
      paras=c(OC,p,cov.type,rho,op,n1,n2,det1,det2)
      paras_name=c('OC','p','cov.type','rho','op','n1','n2','det1','det2')
      ###
      mat_set1=matrix(0,nn1/ns,length(paras))
      for(i in 1:length(paras)){
        mat_set1[,i]=paras[i]
      }
      dfmat1=cbind(Res1,mat_set1)
      colnames(dfmat1)=c('meanMS','meanFDP','meanTPP','sdMS','sdFDP','sdTPP',paras_name)
      nn1=nrow(output1)
      mat_set=matrix(0,nn1,length(paras))
      for(i in 1:length(paras)){
        mat_set[,i]=paras[i]
      }
      dfmat=cbind(output1,mat_set)
      colnames(dfmat)=c('MS','FDP','TPP',paras_name)
      ###
      filename=c()
      for (i in 1:length(paras)){
        temp=paste0(paras_name[i],'-',paras[i])
        filename=paste0(filename,temp,sep='-')
      }
      print(paste0(filename,Sys.time(),sep='-'))
      write.csv(dfmat,file=paste0(filename,'.csv'))
      df.temp=rbind(df.temp,dfmat1)
      }
    }
  }
stopCluster(cl)
proc.time()-tic
print(df.temp)

