rm(list=ls())
setwd("C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II")
source("Functions.R")
library(foreach)
library(doParallel)
ns=100
cod_OC=c('gaussian','t','gamma')
cod_cov=c('ar','bd')
cod_delta=seq(0.1,0.55,0.05)               
df.temp1 <- data.frame()
#---#
num_cl <- detectCores()
cl=makeCluster(num_cl)
registerDoParallel(cl)
#---#
  for(ind.cov in 1:length(cod_cov)){
    for(ind.OC in 1:length(cod_OC)){
            for(ind.del in 1:length(cod_delta)){
              OC=cod_OC[ind.OC]
              cov.type=cod_cov[ind.cov]
              del=cod_delta[ind.del]
              if(cov.type=='ar'){
                rho=0.8;
              }  
              if(cov.type=='bd'){
              del = del +1;
              rho=0.6;
              }            
              n1=400;n2=300;p=600;op=0.1;
              det1=-del
              det2=del
              alpha=c(0.05,0.1,0.2)
              Loop=foreach(l=1:ns,.combine=list,.multicombine = TRUE,.maxcombine = ns,
                           .packages= c("MASS","glmnet",
                                        'foreach','doParallel','mvtnorm'))%dopar%{
                                          ####
                                          set.seed(l+1234)
                                          if(cov.type=='ar'){Sig1<-toeplitz(rho^(0:(p-1)))}
                                          if(cov.type%in%c('bd')){
                                            Sig1<-Gen_Sigma(p,rho,cov.type=cov.type)
                                          }
                                          dat<-Gen_highmean(n1=n1,n2=n2,p=p,op=op,det1=det1,det2=det2,Sig1,OC=OC,OCmu='II')
                                          X1<-dat$X1
                                          Sig1<-dat$Sig1
                                          mu1<-dat$mu1
                                          X2<-dat$X2
                                          Sig2<-dat$Sig2
                                          mu2<-dat$mu2
                                          Ind<-dat$Ind
                                          THETA<-dat$THETA
                                          n1<-nrow(X1);n2<-nrow(X2)
                                          Sig<-1/n1*Sig1+1/n2*Sig2
                                          #oracle omega#
                                          Omega0<-solve(Sig)
                                          Omega0<-(Omega0+t(Omega0))/2
                                          #estimated omega via full data# 
                                          Omega<-Func_omega(X1,X2,Y1,Y2,model='mean')
                                          p=length(THETA)
                                          #--------------#
                                          LEB0=1#oracle omega
                                          LEB=1#estimated omega
                                          SDA0=1#oracle omega, that is Omega0
                                          SDA=1#estimated omega by the first split, that is, Omega==NULL
                                          RES1<-list()
                                          Allres = list(RES1=RES1);k1=0;Method1=NULL;
                                          if(SDA0){#with oracle Omega
                                            k1=k1+1
                                            Det <- Func_SDA2(X1, X2,alpha=alpha,Omega=Omega0, model ='mean')
                                            Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                            Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA0-',x)})))
                                          }
                                          if(SDA){#The first split data estimate omega
                                            k1=k1+1
                                            Det <- Func_SDA2(X1, X2,alpha=alpha, Omega=NULL, model ='mean')
                                            Allres$RES1[[k1]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                            Method1<-c(Method1,do.call(c,lapply(alpha, function(x){paste0('SDA-',x)})))
                                          }
                                          if(LEB0){
                                            k1=k1+1
                                            Det <- Func_LEB2(X1, X2, Omega=Omega0,model='mean')
                                            Allres$RES1[[k1]]<-Func_Measure(Det,Ind,p)
                                            Method1=c(Method1,'LEB0');
                                          }
                                          if(LEB){
                                            k1=k1+1
                                            Det <- Func_LEB2(X1, X2, Omega=Omega,model='mean')
                                            Allres$RES1[[k1]]<-Func_Measure(Det,Ind,p)
                                            Method1=c(Method1,'LEB');
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
              ##
              mat_set1=matrix(0,nn1/ns,length(paras))
              for(i in 1:length(paras)){
                mat_set1[,i]=paras[i]
              }
              dfmat1=cbind(Res1,mat_set1)
              colnames(dfmat1)=c('meanMS', 'meanFDP','meanTPP','sdMS', 'sdFDP','sdTPP',paras_name)
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
stopCluster(cl)
#print(df.temp1)
#write.csv(df.temp1,file = 'Res_Figure2&FigureS1&FigureS2.csv')
 


