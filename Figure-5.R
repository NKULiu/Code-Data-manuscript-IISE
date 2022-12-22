rm(list=ls())
setwd("C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II")
source("Functions.R")
library(foreach)
library(doParallel)
tic=proc.time()
ns = 100
df.temp <- data.frame()
#---------------#
num_cl <- detectCores()
cl=makeCluster(num_cl)
registerDoParallel(cl)
            model='glm';cov.type='ar';p=200;rho=0;op=20;m=1500
            n1<-200
            n2<-200
            del=0.04
            det1=det2=del
            alpha=c(0.05,0.1,0.2)
            ###
            Loop=foreach(l=1:ns,.combine=list,.multicombine = TRUE,.maxcombine = ns,
                         .packages= c("MASS","glmnet",
                                      'foreach','doParallel','mvtnorm'))%dopar%{
                                        maxit=100
                                        set.seed(l+1234)
                                        Sig<-Gen_Sigma(p,rho)
                                        X1<-rmvnorm(m,rep(0,p),Sig)
                                        system.time(dat<-Gen_linearprofile(X1,n1,n2,m,p,rho=0,op,det1,det2,model))
                                        #
                                        Y1<-dat$Y1
                                        Y2<-dat$Y2
                                        Ind<-dat$Ind
                                        system.time(Ome<-Func_omega(X1,X2=NULL,Y1,Y2,model,maxit=maxit))
                                        ###### 
                                        Beta1<-Ome$Beta1### store betas of n1 samples
                                        Beta2<-Ome$Beta2### store betas of n2 samples
                                        S1<-Ome$S1### store Ss of n1 samples
                                        S2<-Ome$S2### store Ss of n2 samples
                                        LEB=1#estimated omega
                                        SDA=1#estimated omega by the first split, that is, Omega==NULL
                                        n_case=c(1:9)
                                        Result<-data.frame()
                                        for(ncase in n_case){
                                          ############
                                          if(ncase==1){n1_r<-200;n2_r<-40}
                                          if(ncase==2){n1_r<-200;n2_r<-60}
                                          if(ncase==3){n1_r<-200;n2_r<-80}
                                          if(ncase==4){n1_r<-200;n2_r<-100}
                                          if(ncase==5){n1_r<-200;n2_r<-120}
                                          if(ncase==6){n1_r<-200;n2_r<-140}
                                          if(ncase==7){n1_r<-200;n2_r<-160}
                                          if(ncase==8){n1_r<-200;n2_r<-180}
                                          if(ncase==9){n1_r<-200;n2_r<-200}
                                          Y1_r<-Y1[,1:(n1_r)]
                                          Y2_r<-Y2[,1:(n2_r)]  
                                          Beta1_r<-Beta1[1:(n1_r),]
                                          Beta2_r<-Beta2[1:(n2_r),]
                                          S1_r<-S1[,,1:(n1_r)]
                                          S2_r<-S2[,,1:(n2_r)]
                                          ###
                                          beta1_r<-apply(Beta1_r,2,mean)
                                          beta2_r<-apply(Beta2_r,2,mean)
                                          beta_r<-beta1_r-beta2_r
                                          s2_r<-s1_r<-matrix(0,p,p)
                                          for(i in 1:n1_r){
                                            s1_r<-s1_r+S1_r[,,i]
                                          }
                                          for(j in 1:n2_r){ 
                                            s2_r<-s2_r+S2_r[,,j]
                                          }
                                          s1_r<-s1_r/n1_r
                                          s2_r<-s2_r/n2_r
                                          Omega_r<-solve(1/n1_r*s1_r+1/n2_r*s2_r)
                                          Allres = list();k=0;Method1=NULL;
                                          ######
                                          if(SDA){
                                            k=k+1
                                            Det <- Func_SDA2(X1, X2=NULL,Y1_r,Y2_r,Beta1=Beta1_r,Beta2=Beta2_r,S1=S1_r,S2=S2_r,alpha=alpha, Omega=NULL, model =model)
                                            Allres[[k]]<-t(sapply(1:length(alpha),function(x){Func_Measure(Det$out[[x]],Ind,p)}))
                                            Method1<-c(Method1,do.call(cbind,lapply(alpha, function(xx){paste0('SDA',xx,'-','n1-',n1_r,'-','n2-',n2_r)})))
                                          }
                                          if(LEB){
                                            k=k+1
                                            system.time(Det <- Func_LEB2(X1, X2=NULL,Y1_r,Y2_r,beta = beta_r, Omega=Omega_r,model=model,maxit=maxit))
                                            Allres[[k]]<-Func_Measure(Det,Ind,p)
                                            Method1=c(Method1,paste0('LEB-','n1-',n1_r,'-','n2-',n2_r));
                                          }
                                          ResSummary1 = lapply(1:length(Allres),function(x){Allres[[x]]})
                                          ResSummary1 = do.call(rbind,ResSummary1);
                                          colnames(ResSummary1)=c('MS',"FDP","TPP")
                                          Result1=as.data.frame(round(ResSummary1,digits = 3))
                                          #######
                                          Result1$Method=Method1
                                          Result=rbind(Result1,Result)
                                        }
                                        return(Result)
                                      }
            #----------------#
            output<-do.call(rbind,lapply(1:ns,function(x){Loop[[x]]}))
            paras=c(model,p,cov.type,rho,op,m,det1,det2)
            paras_name=c('model','p','cov.type','rho','op','m','det1','det2')
            n<-nrow(output)
            MeanRes<-do.call(rbind,lapply(1:(n/ns),function(x){
              #x=1
              dat<-as.data.frame(t(colMeans(output[seq(x,n,n/ns),-ncol(output)]    )))
              rownames(dat)=rownames(output)[x]
              return(dat)
            }))
            sdRes<-do.call(rbind,lapply(1:(n/ns),function(x){
              #x=1
              dat<-as.data.frame(t(apply(output[seq(x,n,n/ns),-ncol(output)],2,sd)))
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
            colnames(dfmat)=c('meanMS','meanFDP','meanTPP','sdMS', 'sdFDP','sdTPP',paras_name)
            dfmat$Method<-unique(output$Method)
            #-----#
            filename=c()
            for (i in 1:length(paras)){
              temp=paste0(paras_name[i],'-',paras[i])
              filename=paste0(filename,temp,sep='-')
            }
            print(paste0(filename,Sys.time(),sep='-'))
            df.temp=rbind(df.temp,dfmat)
  stopCluster(cl)
#############################
aa<-df.temp$Method
name_c<-do.call(cbind,lapply(1:length(aa),function(x){return( strsplit(aa,split='-',fixed=TRUE)[[x]][1]  )}))
n1_c<-do.call(cbind,lapply(1:length(aa),function(x){return( as.numeric(strsplit(aa,split='-',fixed=TRUE)[[x]][3] )  )}))
n2_c<-do.call(cbind,lapply(1:length(aa),function(x){return( as.numeric(strsplit(aa,split='-',fixed=TRUE)[[x]][5] )  )}))
df.temp$Method<-t(name_c)
df.temp$n1<-t(n1_c)
df.temp$n2<-t(n2_c)
proc.time()-tic
# df.temp
# write.csv(df.temp,file = 'Res_Figure4.csv')