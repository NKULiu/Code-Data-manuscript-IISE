#This file plots the W and FDP in Figure 1#
rm(list=ls())
setwd("C:/Users/lyh/Desktop/Revision-IISE/Code&Data-II")
source("Functions.R")
library(ggplot2)
library(reshape)
library(reshape2)
Func_plotW <- function(WWj,tru1,lab_x = ' ',lab_y =' ',fig_title=' '){
  Wj =as.vector(WWj) 
  p_dat=as.data.frame(matrix(0,nrow=length(Wj),ncol=3))
  p_dat[,1]=Wj
  p_dat[,2]=1:length(Wj)
  p_dat[tru1,3]=1
  p_dat[-tru1,3]=0
  names(p_dat)=c("Wj","index","ind")
  p_dat$ind=as.factor(p_dat$ind)
  ylim_range <- ceiling(max(abs(Wj)))
  Wjplot=ggplot(p_dat,aes(x = index, y = Wj, colour = ind,shape=ind,size=ind,alpha=ind))+
    geom_point( )+scale_color_manual(values = c("black","red"))+
    scale_size_manual(values = c(2,2))+
    scale_alpha_manual(values = c(0.5,1))+
    labs(x = lab_x,y=lab_y,title=fig_title) + coord_cartesian(ylim=c(-ylim_range,ylim_range))+
    geom_hline(aes(yintercept=0), linetype="dashed",colour="blue",size=1)+
    theme_bw()+theme(legend.position = "none",panel.grid.major =element_blank(), 
                     panel.grid.minor = element_blank())+
    theme(legend.text=element_text(size=16))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15))+
    theme(axis.title.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20)) +
    theme(plot.title=element_text(hjust=0.5,size = 20))
    return(Wjplot)
}
###plot of FDP(t)###
Func_plotFDP <- function(Wj,tru1,lab_x = ' ',lab_y = ' ',fig_title =' ',leg.poi=c(1,1)){
  t=sort(abs(Wj))
  Ta=sapply(t,function(x){(sum(Wj<=(-x)))/max(1,sum(Wj>=x))})
  bestlam = min(t[which(Ta<=0.2)])
  y1=matrix(0,length(t),1)
  y2=matrix(0,length(t),1)
  for(i in 1:length(t)){
    FZ=which(Wj<=(-t[i]))
    FM=which(Wj>=t[i])
    y1[i]=length(FZ)/max(1,length(FM)) #FDP(t)
    y2[i]=length(setdiff(FM,tru1))/max(1,length(FM)) #FDPW(t)
  }
  FDPbox=as.data.frame(cbind(y1,y2))
  FDPbox0=melt(FDPbox)
  FDPbox0$t=rep(t,times=2)
  FDPbox0$lineind=FDPbox0$variable
  lab=c(expression(widehat(FDP)(t)),expression("FDP(t)"))
  FDPplot=ggplot(FDPbox0,aes(x = t, y = value))+
    geom_line(aes(linetype=variable,color=variable),size=1)+
    scale_linetype_manual(values=c("dashed", "solid"),breaks=c("V1", "V2"),labels=lab)+
    scale_color_manual(values = c("red","black"),breaks=c("V1", "V2"),labels=lab)+
    labs(x=lab_x,y=lab_y,title=fig_title)+
    theme_bw()+theme(legend.position=leg.poi, legend.justification=c(1,1),legend.title=element_blank())+
    theme(legend.text = element_text(size=16))+
    theme(legend.background = element_rect(fill = 'white', colour = 'black'))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15))+
    theme(axis.title.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20)) +
    theme(plot.title=element_text(hjust=0.5,size = 20))
  return(FDPplot)
}
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
tru <- which(THETA!=0)
W=Det$stat[[1]]
W_mean <- Func_plotW(WWj=as.vector(W),tru1=tru,lab_y = 'W',lab_x = 'Index',fig_title = 'Scenario (I)')
FDP_mean <- Func_plotFDP(Wj=as.vector(W),lab_y = 'FDP', lab_x = 't',tru1=tru,fig_title = 'Scenario (I)',leg.poi=c(1,1))
print(W_mean)
print(FDP_mean)
FDP_mean <- FDP_mean + 
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks = c(0,25,50,75))
print(FDP_mean)
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
W <- dat_W$W;
tru <- which(dat_W$ind!=0); 
##Plot the W##
W_glm <- Func_plotW(WWj=as.vector(W),tru1=tru,fig_title = 'Scenario (II)')
##Plot the FDP##
FDP_glm <- Func_plotFDP(Wj=as.vector(W),tru1=tru,fig_title = 'Scenario (II)',leg.poi=c(1,1))+
    coord_cartesian(xlim = c(0,15)) 
print(W_glm)
print(FDP_glm)
###
library(grid)   
library(gridExtra)
gridExtra::grid.arrange(W_mean,W_glm,FDP_mean,FDP_glm,nrow=2,ncol=2) 
#
pdf('Figure-1.pdf',width = 15,height = 9)
gridExtra::grid.arrange(W_mean,W_glm,FDP_mean,FDP_glm,nrow=2,ncol=2) 
dev.off()






