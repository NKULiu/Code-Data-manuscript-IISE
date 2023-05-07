rm(list = ls())
library(R.matlab)  #read matlab file into R#
setwd('C:/Users/lyh/Desktop/mixture/iise-real-data')
Det <- readRDS('detected_oga_sda_leb.RData')
aa<-readMat('fabData.mat')
Y=aa$Y;
X=aa$X;
Ysigma = aa$Ysigma;
Xstage = aa$Xstage;
stageTable = aa$stageTable
Ind_zero <- which(apply(X,2,mean)==0)
stageTable <- stageTable[-Ind_zero,]
###print faulty tools for each wafer sample###
###selected by three methods###
for(i in 1:length(Det)){
  S_dSDA <- which(Det[[i]][,2]!=0)
  S_OGA_new <- which(Det[[i]][,1]!=0)  
  S_LEB <- which(Det[[i]][,3]!=0)
  print(i)
  print(c('dSDA'))
  print(S_dSDA)
  print(c('OGA_new'))
  print(S_OGA_new)
  print(c('LEB'))
  print(S_LEB)
}
##
num <- do.call(rbind,lapply(1:length(Det),
                            function(xx){apply(Det[[xx]],2,sum)}))
##
num = as.data.frame(num)
colnames(num) = c('oga_post','sda_fix','leb');
print(num)
apply(num,2,mean)
Det_sum <- matrix(0,220,3)
for(i in 1:length(Det)){
  Det_sum = Det_sum + Det[[i]]
}
##
Det_sum <- cbind(Det_sum,matrix(stageTable),matrix(c(1:220)))
Det_sum <- as.data.frame(Det_sum)
colnames(Det_sum) = c('OGA_new','dSDA','LEB','Stage','Ind');
###How many tools that are selected?#
sum(Det_sum$OGA_new!=0) #32
sum(Det_sum$dSDA!=0)  #46
sum(Det_sum$LEB!=0)#11
##scatter plot by ggplot##
library(ggplot2)
library(grid)
library(gridExtra)
Det_sum1 <- as.data.frame(matrix(0, nrow = 660, ncol =  4)) 
names(Det_sum1) = c('Frequency','Stage','Ind','fac')
Det_sum1$Frequency <- c(Det_sum$dSDA,Det_sum$OGA_new,Det_sum$LEB)
Det_sum1$Stage <- rep(Det_sum$Stage,3)
Det_sum1$Ind <- rep(Det_sum$Ind,3)
Det_sum1$fac <- rep(c('dSDA','OGA_new','LEB'),each = 220)
Det_sum1$fac = as.factor(Det_sum1$fac)
###
pp <- ggplot(data = Det_sum1, aes(x = Ind, y = Frequency)) + 
  geom_bar(position = position_dodge(),stat = 'identity') +
  coord_cartesian(ylim = c(0, max(Det_sum1$Frequency))) +
  facet_wrap(~fac) + theme(legend.position = 'top') +
  theme(
    legend.title = element_text(size = 20))+
  theme(legend.key.size = unit(20,"pt")) +
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20))+
  theme(plot.title=element_text(hjust=0.5)) +                                                                # Change font size
  theme(strip.text.x = element_text(size = 20))
print(pp)
###