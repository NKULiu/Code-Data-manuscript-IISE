rm(list=ls())
setwd('C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II')
library(doParallel)
source('Functions.R')
library(Matrix)
library(RNOmni) #RankNorm for normal transformation#
secom<-read.table("secom.data.txt", sep = " ")
secom_labels<-read.table('secom_labels.data.txt',sep=' ')
#The ratio of null values# 
sum(is.na(secom))/(dim(secom)[1]*dim(secom)[2])
#First, impute na's with mean values#  
for(i in 1:dim(secom)[2]){
  secom[,i][is.na(secom[,i])] = mean(secom[,i], na.rm = TRUE)
}
#Second, remove discrete columns with less than 2 values#
num_Discrete <- apply(secom,2,function(x)length(unique(x)))
sum(num_Discrete<=2)  #116
ind_Discrete<- which(num_Discrete<=2)
secom <- secom[,-ind_Discrete]
dim(secom)
#
IC_index<-which(secom_labels$V1==-1)
OC_index<-which(secom_labels$V1==1)
#
secom_IC <- secom[IC_index,]
secom_OC <- secom[OC_index,]
#Third remove the constant columns#
sd1 <- apply(secom_IC,2,function(x)sd(x,na.rm = T))
sum(sd1<=1e-12)
#
sd2 <- apply(secom_OC,2,function(x)sd(x,na.rm = T))
sum(sd2<=1e-12)
#
sds <- apply(secom,2,function(x)sd(x,na.rm = T))
sum(sds<=1e-12)
##
ind_sd0 <- which(sd2<=1e-12)
if(length(ind_sd0)>0){
  secom <- secom[,-ind_sd0]
}
#
dim(secom) #1567 by 468#
secom_IC <- secom[IC_index,];
secom_OC <- secom[OC_index,];
###normal IC###
dat.IC <- do.call(cbind,lapply(1:dim(secom_IC)[2],function(x){
  return(RankNorm(secom_IC[,x], k = 0.5, ties.method = "average"))
}))
###normal OC###
dat.OC <- do.call(cbind,lapply(1:dim(secom_IC)[2],function(x){
  FF <- (rank(secom_IC[,x], na.last = "keep",
              ties.method ='average') - 0.5) / sum(!is.na(secom_IC[,x]))
  Fx = approx(secom_IC[,x],FF,xout = secom_OC[,x],rule=2,method = 'linear',ties = mean)$y
  return((qnorm(Fx,mean=0,sd=1)))
}))
n.IC <- nrow(dat.IC);
n.OC <- nrow(dat.OC);
p <- ncol(dat.IC);
#estimate covariance S#
library(corpcor)
S.IC <- matrix(as.numeric(corpcor::cov.shrink(dat.IC,verbose=T)), nrow=p )
S.IC <- Matrix::nearPD(S.IC, conv.tol = 1e-15, conv.norm = "F")$mat
S.IC = as.matrix(S.IC)
##
S.OC <- matrix(as.numeric(corpcor::cov.shrink(dat.OC,verbose=T)), nrow=p )
S.OC <- Matrix::nearPD(S.OC, conv.tol = 1e-15, conv.norm = "F")$mat
S.OC = as.matrix(S.OC)
##
S <- 1/n.IC*S.IC +1/n.OC*S.OC
S <- Matrix::nearPD(S, conv.tol = 1e-15, conv.norm = "F")$mat
S = as.matrix(S)
image(S,col = gray.colors(256))
###estimate Omega###
nlambda = 20;
lambda.min.ratio = 1e-3
lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
lambda.min = lambda.min.ratio * lambda.max
lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
out <- rho.glasso(S,n=1/(1/n.IC+1/n.OC),lambda = lambda,nlambda=20,crit = 'AIC',parallel = T)   
Ome1=out$Omega
# image(Ome1,col=gray.colors(256))
# sum(Ome1!=0)
# lambda.opt <- out$lambda.opt
# print(lambda.opt)
alpha <- c(0.1,0.15,0.2,0.25,0.3)
###The SDA###
Det_SDA <- Func_SDA2(dat.IC, dat.OC,alpha=alpha,Omega=Ome1, model ='mean')$out
###The LEB###
Det_LEB <- Func_LEB2(dat.IC, dat.OC, Omega = Ome1,model='mean')
###MDW###
Det_MDW<-Func_MDW(dat.OC,lam = 0.1,b_value = 5,alpha = alpha)$Res
###replicate the dSDA to calculate the average number detected by SDA#
repli <- 20
library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
system.time(Loop<-foreach(l=1:repli,.combine=list,.multicombine = TRUE,.maxcombine = repli,
                          .packages= c("MASS","glmnet",'foreach','doParallel','mvtnorm'
                          ))%dopar%{
                            set.seed(l+1234)
                            Det <- Func_SDA2(dat.IC, dat.OC,alpha=alpha, Omega=Ome1, model ='mean')$out
                            return(Det)
                          })
stopCluster(cl)  
##
number_dsda1 <- matrix(0,nrow= length(alpha),ncol = repli)
for(j in 1:repli){
  #j=1
  Det <- Loop[[j]]
  number_dsda1[,j] = do.call(c,lapply(Det,length))
}
number_dsda = round(apply(number_dsda1,1,mean))
number_LEB = length(Det_LEB)  #Only detect 1 coordinate#
number_MDW = sapply(Det_MDW,length)
numbers <- rbind(number_dsda,number_LEB,number_MDW)
#Table 4#
xtable::xtable(numbers)
###Visualize the oc mean of dSDA and MDW###  
num_max <- apply(number_dsda1,1,which.max)
Det_SDA <- Loop[[num_max[2]]]
###plot the results###
OC_mean = apply(dat.OC,2,mean)
dat <- as.data.frame(matrix(0,p,5))
dat[,1] = 1:p;
dat[,2] = OC_mean;
dat[,-c(1,2)] = 'IC';
#dat[Det_LEB,3]  = 'OC';
dat[Det_SDA[[2]],3] = 'OC'; #alpha = 0.15;
dat[Det_MDW[[2]],4] = 'OC';  #alpha = 0.15;
names(dat) <- c('Index','Value','dSDA','MDW')
##
dat$`MDW` = as.factor(dat$`MDW`);
#dat$LEB = as.factor(dat$LEB);
dat$`dSDA` = as.factor(dat$`dSDA`);
###
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
cbbPalette<-brewer.pal(9,'Set1')
Func_p<-function(df,x,x.title,y,y.title,title,shap,colo,text.xy=19,title.xy=21,
                 legend.key.size=3,legend.text.size=19,legend.title.size=21,
                 strip.text.size=19){
  p_1 <- ggplot(df,aes(x=x,y=y,shape = shap ,colour=colo))+
    geom_point(size=2)+
    labs(x=x.title,y=y.title,title = title)+
    theme(axis.text.x = element_text(size=text.xy),axis.title.x=element_text(size=title.xy))+
    theme(axis.text.y = element_text(size=text.xy),axis.title.y=element_text(size=title.xy))+
    theme(legend.position = "none",legend.key.size =unit(legend.key.size,'line'),
          legend.title=element_text(size=legend.title.size),
          legend.text=element_text(size=legend.text.size))+
    theme(
      strip.text.x = element_text(margin = margin(0, 0, 0, 0))
    )+
    theme(strip.text.x = element_text(size = strip.text.size),strip.text.y = element_text(size=strip.text.size))+
    theme(legend.title = element_blank())+
    theme(legend.position = 'none')+
    theme(plot.title=element_text(hjust=0.5))+
    theme(title = element_text(size = title.xy))+
    scale_color_manual(values = cbbPalette[c(1,3)]) 
  return(p_1)
}
p1<-Func_p(dat,x=dat$Index,'Index',y=dat$Value,'Value','(a) dSDA',shap = dat$dSDA,colo = dat$dSDA)
p2<-Func_p(dat,x=dat$Index,'Index',y=dat$Value,'Value','(b) MDW',shap = dat$MDW,colo = dat$MDW)
###
#grid.arrange(p1,p2,nrow=1,ncol=2)
#Figure 6#
#pdf('Real-Figure.pdf',width = 15,height = 6)
pdf('Figure-6.pdf',width = 15,height = 6)
grid.arrange(p1,p2,nrow=1,ncol=2)
dev.off()

#visualize covariance before and after the change point# 
##Figure S4##
#pdf('covariance.pdf',width = 15,height = 6)
pdf('Figure-S4.pdf',width = 15,height = 6)
par(mfrow=c(1,2))
image(S.IC,main = 'IC covariance matrix')
image(S.OC,main = 'OC covariance matrix')
par(mfrow=c(1,1))
dev.off()
###artificially contaminate the IC data### 
###and do a simulation study###
cod_delta <- seq(0.1,0.5,0.05);
cod_op <- c(0.1,0.2)
alpha = c(0.2)
rep <- 100
dat_measure <- data.frame()
num_cl <- detectCores(all=T)
cl <- makeCluster(num_cl)
registerDoParallel(cl)
for(ind.delta in 1:length(cod_delta)){
  for(ind.op in 1:length(cod_op)){
    #      
    del <- cod_delta[ind.delta]
    op <- cod_op[ind.op]
    N_OC<-p*op
    #
    Loop <- foreach(ii=1:rep,
                    .combine = rbind,
                    .maxcombine = rep,
                    .packages = c('glmnet','ks','corpcor',
                                  'pracma','glassoFast','mvtnorm')
    )%dopar%{
      set.seed(ii+1234)
      #randomly choose OC data streams#
      ind_OC=sample(1:p,N_OC)
      mu2 <- rep(0,p);
      mu2[ind_OC] = rnorm(length(ind_OC),del,0.05)
      X.OC <- rmvnorm(n.OC,mu2,S.IC)
      #assume only mean shift before and after change point#
      S.OC <- S.IC;
      S<-1/n.IC*S.IC+1/n.OC*S.OC
      ##set lambda
      nlambda = 10;
      lambda.min.ratio = 1e-3
      lambda.max = 0.1
      lambda.min = lambda.min.ratio * lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      ##estimate Omega##
      out <- rho.glasso(S,n=1/(1/n.IC+1/n.OC),lambda = lambda,crit = 'AIC') 
      Ome=out$Omega
      Det_SDA <- Func_SDA2(dat.IC, X.OC, alpha=alpha, Omega=Ome, model ='mean',K=10)$out
      Det_LEB <- Func_LEB2(dat.IC, X.OC, Omega=Ome, model='mean')
      Det_MDW <- Func_MDW(X.OC,lam = 0.1,b_value = 5,alpha = alpha)$Res
      ###
      measure_MDW <- t(sapply(1:length(alpha),function(x){
        return(c(Func_Measure(Det_MDW[[x]],ind_OC,p),alpha[x]))
      }))
      measure_SDA <- t(sapply(1:length(alpha),function(x){
        return(c(Func_Measure(Det_SDA[[x]],ind_OC,p),alpha[x]))
      }))
      measure_LEB <- t(sapply(1:length(alpha),function(x){
        return(c(Func_Measure(Det_LEB,ind_OC,p),alpha[x]))
      }))
      measures <- as.data.frame(rbind(measure_LEB,measure_SDA,measure_MDW))
      colnames(measures) = c('MS','FDP','TDP','alpha')
      measures$method = rep(c('LEB','dSDA','MDW'),each=length(alpha))
      measures$delta = del;
      measures$op = op;
      return(measures)
    }
    dat_measure <- rbind(dat_measure,Loop)
  }
}     
#write.csv(dat_measure,file ='real-dat.csv')
#The results are stored in real-dat.csv# 
dat_measure <- read.csv('real-dat.csv')
dat_measure$method[dat_measure$method=='leb']='LEB';
dat_measure$method[dat_measure$method=='sda']='dSDA';
dat_measure$method[dat_measure$method=='mdr']='MDW';
##And we plot the boxplots##
df <- dat_measure[,-c(1,2)]
df_fdp<-df[,-c(2)]
df_fdp$measure<-'FDP'
df_tpp<-df[,-c(1)]
df_tpp$measure<-'TDP'
names(df_fdp)[1]=names(df_tpp)[1]='value'
###
df<-rbind(df_fdp,df_tpp)
df$measure <- as.factor(df$measure);
df$delta <- as.factor(df$delta);
df$op <- as.factor(df$op);
#
Func_Box<-function(df,x,x.title,y,y.title='TDR',facet1,facet2,alpha,linesize=1,text.xy=15,title.xy=17,
                   legend.key.size=2,legend.text.size=15,legend.title.size=17,labeller=label_value,
                   strip.text.size=15,space='free',scales='free'){
  p <- ggplot(df,aes(factor(x),y=y))+ 
    geom_boxplot(aes(fill=Method))+
    scale_fill_manual(values=cbbPalette)+
    labs(x=x.title,y=y.title)+
    theme(axis.text.x = element_text(size=text.xy),axis.title.x=element_text(size=title.xy))+
    theme(axis.text.y = element_text(size=text.xy),axis.title.y=element_text(size=title.xy))+
    theme(legend.position = "top",legend.key.size =unit(legend.key.size,'line'),
          legend.title=element_text(size=legend.title.size), 
          legend.text=element_text(size=legend.text.size))+
    facet_grid(formula(paste(facet2,"~",facet1)),space=space,scales=scales,labeller=labeller)+
    theme(strip.text.x = element_text(size = strip.text.size),strip.text.y = element_text(size=strip.text.size))+
    geom_hline(data=df[df$measure=='FDP',],aes(yintercept = alpha), colour="black", linetype="dashed")+
    theme(legend.title = element_text(size = title.xy))+
    scale_color_manual(values=cbbPalette)
  return(p)
}
##
names(df)[3] = 'Method';
df1 <- df[df$alpha==0.2&df$delta%in%c(0.15,0.2,0.25,0.3,0.35,0.4),]
levels(df1$op) = c('OC ratio: 10%','OC ratio: 20%')
Box_p1<-Func_Box(df1,x=df1$delta,x.title = expression(gamma),y=df1$value,y.title=" ",
                 facet1 = 'measure',facet2 = 'op',linesize=1,text.xy=17,title.xy=19,legend.key.size=2,legend.text.size=17,
                 legend.title.size=19,labeller=label_value,strip.text.size=17,space='free_x',scales='free')
print(Box_p1)
###Figure 7###
pdf('Figure-7.pdf',width = 15,height = 9)
Box_p1
dev.off()
