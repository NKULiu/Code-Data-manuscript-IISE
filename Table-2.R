rm(list=ls())
setwd("C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II")
source("Functions.R")
library(foreach)
library(doParallel)
ns = 100
cod_p=c(400,600)
cod_op=c(30,60)
cod_ncase=c(1,2)
df.temp <- data.frame()
#----#
num_cl <- detectCores()
cl=makeCluster(num_cl)
registerDoParallel(cl)
for(ind.ncase in 1:length(cod_ncase)){
  for(ind.p in 1:length(cod_p)){
    for(ind.op in 1:length(cod_op)){
      p=cod_p[ind.p]
      op=cod_op[ind.op]
      n_case=cod_ncase[ind.ncase]
      OC=1;Ak=Ck=1;del=0.55;
      #####
      if(n_case==1){
        n1=400;n2=400
      }
      if(n_case==2){
        n1=400;n2=200
      }

      E_n=n1*n2/(n1+n2)
      delta=del
      det1=del
      det2=del
      alpha=c(0.05)
      ###
      a0=0
      sig_epsilon=sig_vk=sig_wk=1
      ###
      multi_cov<-Gen_multicov(p,Ak,Ck,a0,sig_epsilon,sig_vk,sig_wk)
      Gam<-multi_cov$Gam
      invGam<-multi_cov$invGam
      covY<-multi_cov$covY
      covX<-multi_cov$covYY
      Loop=foreach(l=1:ns,.combine=list,.multicombine = TRUE,.maxcombine = ns,
                   .packages= c("MASS","glmnet",
                                'foreach','doParallel','mvtnorm'))%dopar%{
                                  ##############
                                  set.seed(l+1234)
                                  #generate IC sample X1 and OC sample X2#
                                  dat <- Gen_multistage(n1=n1,n2=n2,p=p,op=op,det1,det2,Ak=Ak,Ck=Ck,a0=0,sig_epsilon=1,sig_vk=1,sig_wk=1)
                                  ###
                                  X1<-dat$X1
                                  X2<-dat$X2
                                  e1<-dat$e1
                                  e2<-dat$e2
                                  Ind<-dat$Ind
                                  delta<-dat$delta
                                  beta<-invGam%*%delta##The oracle beta##
                                  S<-1/n1*covX+1/n2*covX
                                  Omega0<-solve(S)
                                  #avoid Omega0 is not symmetric#
                                  Omega0<-(Omega0+t(Omega0))/2
                                  if(!(matrixcalc::is.positive.definite(Omega0))){
                                    Omega0<-as.matrix(Matrix::nearPD(Omega0, conv.tol = 1e-15, conv.norm = "F")$mat)
                                  }
                                  #Omega<-Func_omega(X1,X2,Y1,Y2,model='multistage',
                                  #                 cov.type,invGam,K_hat=1,stars=F)
                                  LEB0=1
                                  SDA0=1
                                  #--------------#
                                  Allres = list();k=0;Method1=NULL;
                                  ###
                                  if(SDA0){
                                    k=k+1
                                    Det <- Func_SDA2(X1, X2, alpha=alpha,Omega=Omega0, invGam = invGam,model ='multistage')
                                    Allres[[k]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                    Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA0',x)})))
                                  }
                                  # 
                                  if(LEB0){
                                    k=k+1
                                    Det <- Func_LEB2(X1, X2,Omega=Omega0,invGam=invGam,model='multistage')
                                    Allres[[k]]<-Func_Measure(Det,Ind,p)
                                    Method1=c(Method1,'LEB0');
                                  }
                                  ResSummary1 = lapply(1:length(Allres),function(x){Allres[[x]]})
                                  ResSummary1 = do.call(rbind,ResSummary1);
                                  rownames(ResSummary1)=Method1;
                                  colnames(ResSummary1)=c('MS','FDP','TPP')
                                  Result1=round(ResSummary1,digits = 3)
                                }
      #----------------#
      output<-do.call(rbind,lapply(1:ns,function(x){Loop[[x]]}))
      paras=c(OC,p,op,n1,n2,delta)
      paras_name=c('OC','p','op','n1','n2','delta')
      n<-nrow(output)
      MeanRes<-do.call(rbind,lapply(1:(n/ns),function(x){
        #x=1
        dat<-as.data.frame(t(colMeans(output[seq(x,n,n/ns),])))
        rownames(dat)=rownames(output)[x]
        return(dat)
      }))
      sdRes<-do.call(rbind,lapply(1:(n/ns),function(x){
        #x=1
        dat<-as.data.frame(t(apply(output[seq(x,n,n/ns),],2,sd)))
        rownames(dat)=rownames(output)[x]
        return(dat)
      }))
      MeanRes <- round(MeanRes,3)
      sdRes <- round(sdRes,3)
      Res <- cbind(MeanRes,sdRes)
      #
      mat_set=matrix(0,n/ns,length(paras))
      for(i in 1:length(paras)){
        mat_set[,i]=paras[i]
      }
      dfmat=cbind(Res,mat_set)
      colnames(dfmat)=c('meanMS','meanFDP','meanTPP','sdMS','sdFDP','sdTPP',paras_name)
      #-----#
      filename=c()
      for (i in 1:length(paras)){
        temp=paste0(paras_name[i],'-',paras[i])
        filename=paste0(filename,temp,sep='-')
      }
      print(paste0(filename,Sys.time(),sep='-'))
      df.temp=rbind(df.temp,dfmat)
    }
  }
}
stopCluster(cl)
#print(df.temp)
#write.csv(df.temp,file = 'Res_Table3.csv')
