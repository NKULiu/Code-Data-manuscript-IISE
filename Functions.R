#############################################
###The following function 'Func_SDA2' is mainly adapted from the functions in the R package 'sdafilter' constructed by Lilun Du, Xu Guo, Wenguang Sun, and Changliang Zou
###which is available from 'https://github.com/dulilun/sdafilter'.
#####################################
library(MASS)
library(mvtnorm)
library(glasso)
library(glassoFast)
library(glmnet)
library(lars)
library(Matrix)###'nearPD'---make the covariance matrix positive definite
library(matrixcalc)#is.positive.definite to check positiveness 
library(fastglm)###fast glm###
library(selectiveInference)
library(ks)
library(pracma)
##
Func_TPR<-function(THETA,stat,reject.ratio){
  tru=which(THETA!=0)
  #rank the stat
  or=order(stat,decreasing=T)
  cellnum=sum(THETA!=0)
  #where to cut
  Res<-as.data.frame(sapply(reject.ratio,function(x){
    ##
    num_reject=cellnum*x
    indrejectedcell=or[1:num_reject]
    det.in.reject=intersect(tru,indrejectedcell)
    TPR=length(det.in.reject)/length(tru)#the true signal/overall signal 
    return(c(TPR))
  }))
  names(Res)=c('TPR')
  return(Res)
}
####
Func_Measure<-function(det,Ind,p){
  Modelsize= length(det)
  FDP= length(setdiff(det,Ind))/max(1,length(det))
  TPP= length(intersect(det,Ind))/max(1,length(Ind))
  return(c(Modelsize,FDP,TPP))
}
######mean and covariance of standard M estimator#####
Func_betaS<-function(X,Y,model='lm',maxit=150){
  func_betaS<-function(X,Y,model='lm',maxit=1){
    if(model%in%c('lm')){
      A <- solve(crossprod(X))
      beta <- as.vector(A %*% t(X) %*% Y)#########as.vector(matrix(0,p*q,p*q))
      residuals <- Y - X %*% beta ###residuals
      p <- ncol(X)##do not consider intercept
      n <- nrow(X)
      residual_var <- as.numeric(crossprod(residuals) / (n - p))
      S<- residual_var * A
    }
    if(model=='rlm2'){
      mod<-MASS::rlm(Y~X-1,maxit=maxit)
      S<-vcov(mod)
      beta<-mod$coefficients
    }    
    if(model=='glm'){
      mod<-fastglm::fastglm(X,Y,family=binomial(),maxit=maxit,method=2)#The iteration is high to guarantee convergence
      beta<-mod$coefficients
      ###write variance-covariance matrix estimation by myself###
      piX = mod$fitted.values
      W = diag(piX * (1 - piX))
      infMatrix = t(X) %*% W %*% X 
      S <- solve(infMatrix)
    }
    return(list(beta=beta,S=S))
  }
  ######
  p<-dim(X)[2]
  Ydim<-dim(Y)
  ndim<-length(Ydim)
  n<-Ydim[ndim]
  ######
  BETAS<-lapply(1:n,function(x){
    if(ndim<3){
      y<-Y[,x]
    }else{
      y<-Y[,,x]
    }
    mod<-func_betaS(X,y,model,maxit = maxit)
    return(mod)
  })
  return(BETAS)
}
#####################help functions###############
# Omega^{1/2}
Sqrt <- function(Omega_1){
  EgSig0 = eigen(Omega_1)
  EgVec = EgSig0$vectors
  Gamma = EgVec%*%diag( sqrt(EgSig0$values) )%*%t(EgVec)
  return( Gamma )
}
##########
chole<-function(sample,band){
  dim=ncol(sample)
  s_size=nrow(sample)
  sample_cen<-sample-matrix(colMeans(sample),ncol=dim,nrow=s_size,byrow=T)
  d<-NULL; d[1]<-var(sample_cen[,1])
  a<-matrix(0,ncol=dim,nrow=dim)
  #estimate a_i and d_i^2 using least square regression
  #band=1
  for(j in 2:dim){
    X_mat<-sample_cen[,(max(j-band,1)):(j-1)]
    y<-sample_cen[,j]
    beta_hat<-solve(crossprod(X_mat),crossprod(X_mat,y))
    a[j,(max(j-band,1)):(j-1)]<-beta_hat
    d[j]<-1/(s_size)*sum((y-sample_cen%*%a[j,])^2)	
  }
  inv_hat<-t(diag(1,dim)-a)%*%solve(diag(d,dim))%*%(diag(1,dim)-a)
  sig_hat<-solve(diag(1,dim)-a)%*%diag(d,dim)%*%solve(t(diag(1,dim)-a))
  return(list(Omega=inv_hat,Sigma=sig_hat))
}
#############################
###output beta and omega when model is lm, rlm, glm###  
#############################
Func_omega<-function(X1=NULL,X2=NULL,
                     Y1=NULL,Y2=NULL,
                     model='mean',maxit=150){
  ##########################
  if(model=='mean'){
    n1=nrow(X1);
    n2<-nrow(X2);
    ######
    S1 <- chole(X1,band=1)$Sigma
    S2 <- chole(X2,band=1)$Sigma
    Omega <- solve(1/n1*S1+1/n2*S2)
   
  }else if(model%in%c('lm','rlm2','glm')){
    ####
    ndim<-length(dim(Y1))
    n1<-dim(Y1)[ndim]
    n2<-dim(Y2)[ndim]
    p<-dim(X1)[2]
    ##########beta#########
    mod1<-Func_betaS(X1,Y1,model,maxit)
    mod2<-Func_betaS(X1,Y2,model,maxit)
    Beta1<-do.call(rbind,lapply(1:n1,function(x){return(mod1[[x]]$beta)}))
    SS1<-array(unlist(lapply(1:n1,function(x){return(mod1[[x]]$S)})),dim=c(p,p,n1))
    ###
    Beta2<-do.call(rbind,lapply(1:n2,function(x){return(mod2[[x]]$beta)}))
    SS2<-array(unlist(lapply(1:n2,function(x){return(mod2[[x]]$S)})),dim=c(p,p,n2))
    ##mean the beta and S##
    beta1<-apply(Beta1,2,mean)
    beta2<-apply(Beta2,2,mean)
    #####
    S1<-matrix(0,p,p)
    for(i in 1:n1){
      S1<-S1+SS1[,,i]
    }
    #####
    S2<-matrix(0,p,p)
    for(j in 1:n2){
      S2<-S2+SS2[,,j]
    }
    S1<-S1/n1
    S2<-S2/n2
    ########
    beta<-beta1-beta2
    Omega <-solve(1/n1*S1+1/n2*S2)
  }
  ###make sure Omega is positive definite###
  if(!(matrixcalc::is.symmetric.matrix(Omega))){
    Omega<-(Omega+t(Omega))/2
  }
  if(!(matrixcalc::is.positive.definite(Omega))){
    Omega<-as.matrix(Matrix::nearPD(Omega, conv.tol = 1e-15, conv.norm = "F")$mat)
  }
  if(model%in%c('mean')){
    return(Omega)
  }else{
    #S1, Beta1 is the overall beta and S of n1 samples;
    #S2, Beta2 is the overall beta and S of n2 samples;
    #beta is the difference of beta
    return(list(Beta1=Beta1,Beta2=Beta2,S1=SS1,S2=SS2,beta=beta,Omega=Omega))
  }
}
###BIC for choosing tuning parameter for glasso###
rho.glasso<-function(S,n,lambda=NULL,nlambda=50,crit="AIC",g=0.5,parallel = F){
  p=ncol(S)
  if(is.null(lambda)){
    lambda=sqrt(log(p)/n)*seq(0.001,0.5,length=nlambda)
  }
  if(parallel){
    num_core <- detectCores(all=T)
    cl <- makeCluster(num_core)
    registerDoParallel(cl)
    library(foreach)
    library(doParallel)
    myscore<-foreach(x=1:length(lambda),.combine = c,
                     .maxcombine = length(lambda), 
                     .packages = c('glassoFast'))%dopar%{
                       #x=1
                       lamb=lambda[x]
                       out<-glassoFast(S,rho=lamb)
                       Ome=out$wi
                       loglike=sum(diag(S%*%Ome))-determinant(Ome,logarithm = TRUE)$modulus
                       df=(sum(Ome!=0)-p)/2
                       if (crit== "AIC") {
                         CV_ERR=loglike+df*2/n
                       }
                       if (crit== "BIC") {
                         CV_ERR=loglike+df*log(n)/n
                       }
                       if(crit=='EBIC'){
                         CV_ERR=loglike+df*(log(n)+4*g*log(p))/n
                       }
                       return(CV_ERR)
                     }
    stopCluster(cl)
  }else{
    myscore<-sapply(1:length(lambda),function(x){
      #x=1
      lamb=lambda[x]
      out<-glassoFast(S,rho=lamb)
      Ome=out$wi
      loglike=sum(diag(S%*%Ome))-determinant(Ome,logarithm = TRUE)$modulus
      df=(sum(Ome!=0)-p)/2
      if (crit== "AIC") {
        CV_ERR=loglike+df*2/n
      }
      if (crit== "BIC") {
        CV_ERR=loglike+df*log(n)/n
      }
      if(crit=='EBIC'){
        CV_ERR=loglike+df*(log(n)+4*g*log(p))/n
      }
      CV_ERR
    })
  }
  opt.i=which.min(myscore)
  lambda.opt=lambda[opt.i]
  out<-glassoFast(S,rho=lambda.opt)
  Omega=out$wi
  return(list(lambda.opt=lambda.opt,Omega=Omega))
}
###mu2 is the oracle mean vector to serve as weight, and THETA is the oracle nonzero index###   
Func_MOW<-function(X,op,op1,det1,det2,sig2,alpha){
  n2<-dim(X)[1]
  p<-dim(X)[2]
  sigma=1/sqrt(n2)
  if(n2>1){X=apply(X,2,mean)}
  g=op1*dnorm((X-det1)/sqrt(sigma^2+sig2))/sqrt(sigma^2+sig2)+
    (1-op1)*dnorm((X-det2)/sqrt(sigma^2+sig2))/sqrt(sigma^2+sig2);#g1
  #
  f=(1-op)*dnorm(X,0,sigma)+op*g;#g
  ##
  #if(zmu){
  partI=op1*(sigma^2*sig2/(sigma^2+sig2)+((sig2*X+sigma^2*det1)/(sigma^2+sig2))^2)*dnorm(X,det1,sqrt(sigma^2+sig2))
  partII=(1-op1)*(sigma^2*sig2/(sigma^2+sig2)+((sig2*X+sigma^2*det2)/(sigma^2+sig2))^2)*dnorm(X,det2,sqrt(sigma^2+sig2))
  HX=partI+partII
  #}else{
   # HX=0  
  #}
  LAMBDA=(1-op)*dnorm(X,0,sigma)/((1-op)*dnorm(X,0,sigma)+op*(HX+g));
  ##
  DX=op*(HX+g)/f#D(X)
  #LAM_sort<-sort(LAMBDA)
  ind<-order(LAMBDA)
  DX_sort=DX[ind]
  result<-rep(0,p)
  for(j in 1:p){
    result[j]=sum(DX_sort[(j+1):p])/sum(DX_sort) 
  }
  Res<-lapply(alpha,function(x){
    cut<-which((result<=x)!=0)[1]
    delta<-rep(0,p) #decision set
    delta[ind[1:cut]]=1
    return(which(delta!=0))
  })
  return(list(stat=LAMBDA,Res=Res))
}
###
gmu2<-function(p,n,lam,X,mu){
  sigma<-1/sqrt(n);
  Time=linspace(-1/lam,1/lam,100);
  resu<-sapply(Time,function(x){
    return(mean(cos(x*(mu-X)))*exp(0.5*sigma^2*x^2)/lam/50)
  })
  result<-sum(resu)
  result=result-2*(1-p)*sin(mu/lam)/mu;
  result=result/(2*pi*p);
  g=max(result,0);
  return(g)  
}
###
epsest.func <- function(x,u,sigma)
{
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt = linspace(0,tmax,100);
  epsest<-sapply(1:length(tt),function(xxx){
    t  = tt[xxx];
    f  = t*xi;
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co<-sapply(1:length(xi),function(yyy){
      return(mean(cos(t*xi[yyy]*z)))
    })
    return(1 - sum(w*f*co)/sum(w))
  })
  return(epsest = max(epsest))
}
###
Func_MDW<-function(X2,lam=0.1,b_value=5,alpha){
  n<-dim(X2)[1]
  p<-dim(X2)[2]
  sigma=1/sqrt(n)
  if(n>1){X=apply(X2,2,mean)}
  phat <- epsest.func(X,0,sigma)
  if(phat==0){phat=1e-8};
  zv.ds <- kde(X,density = T)
  f = predict(zv.ds,x=X)
  bound <- b_value;
  n_bound <- 200
  MU<-linspace(-bound,bound,n_bound)
  GMU <- sapply(MU,function(x){
    return(gmu2(phat,n,lam,X,x))
  })
  wi<-sapply(X,function(x){
    int<-sum(sapply(1:(length(MU)-1),function(ii){
      return(GMU[ii]*(1+MU[ii]^2)*dnorm(x,MU[ii],1/sqrt(n))/n_bound*2*bound)
    }))
    return(int)     
  })
  #plot(X,wi)  
  LAMBDA=(1-phat)*dnorm(X,0,sigma)/((1-phat)*dnorm(X,0,sigma)+phat*wi);
  DX=phat*wi/f;
  ind<-order(LAMBDA)
  DX_sort=DX[ind]
  result <- sapply(1:p,function(x){
    sum(DX_sort[x:p])/sum(DX_sort)
  })
  Res<-lapply(alpha,function(x){
    cut<-which((result<=x)!=0)[1]
    delta<-rep(0,p)
    delta[ind[1:cut]]=1 
    return(which(delta!=0))
  })
  return(list(stat=LAMBDA,Res=Res))
}
####S1 and S2 only useful when need to estimate omega via the first split data#### 
Func_SDA2 <- function(X1=NULL,X2=NULL,Y1=NULL,Y2=NULL,
                      scale=T,beta=NULL,Beta1=Beta1,Beta2=Beta2,S1=NULL,S2=NULL,
                      alpha=0.2,Omega=NULL,invGam=NULL,stable=TRUE,K=10,model='mean',
                      maxit=150){
  ##
  Model_Select <- function(Y, X){
    # Method I: AIC
    fit <- glmnet::glmnet(X, Y, family = "gaussian")
    p<-dim(X)[2]
    k <- fit$df
    AIC <- stats::deviance(fit)+2*k
    i_min = which.min(AIC)
    lambda_select = fit$lambda[i_min]
    fit_AIC= coef(fit,s=lambda_select)[-1] ####remove the intercept###
    w1 = fit_AIC
    sv = as.vector( which(w1!=0) )
    #method II: upper bound p/3
    k= ceiling(p/3)
    if(length(sv)>k){
      wv = which(fit$df==max(fit$df[fit$df<k]))[1]
      sv = which(fit$beta[, wv]!=0)
      w1 = fit$beta[, wv]
      lambda_select<-fit$lambda[wv]
    }
    return( list(index_S=sv, w1 = w1,lambda_select=lambda_select))
  }
  ##
  W_value <- function(X1,X2,Y1,Y2,scale,beta,Beta1,Beta2,S1,S2,Index1,Index2,Omega,invGam,model,maxit){
    #-----------#
    if(model=='mean'){
      X1_1=X1[Index1,];
      X1_2=X1[-Index1,]
      X2_1=X2[Index2,];
      X2_2=X2[-Index2,]
      n1_1 <- nrow(X1_1)
      n1_2 <- nrow(X1_2) 
      n2_1 <- nrow(X2_1)
      n2_2 <- nrow(X2_2) 
      n_1 = n1_1+n2_1
      n_2 = n1_2+n2_2
      ##compute beta####  
      beta1_1<-apply(X1_1,2,mean);
      beta2_1<-apply(X2_1,2,mean)
      beta1_2<-apply(X1_2,2,mean);
      beta2_2<-apply(X2_2,2,mean)
      beta1<-beta1_1-beta2_1;
      beta2<-beta1_2-beta2_2
      ##compute Omega###
      if(is.null(Omega)){
        Omega<-Func_omega(X1_1,X2_1,Y1=NULL,Y2=NULL,model)
        Gamma = Sqrt(Omega)
      }else{
        Gamma<-Sqrt(Omega)
      }
    }else if(model=='cov'){
      X1_1=X1[Index1,];
      X1_2=X1[-Index1,]
      X2_1=X2[Index2,];
      X2_2=X2[-Index2,]
      n1_1 <- nrow(X1_1)
      n1_2 <- nrow(X1_2) 
      n2_1 <- nrow(X2_1)
      n2_2 <- nrow(X2_2) 
      n_1 = n1_1+n2_1
      n_2 = n1_2+n2_2
      #
      S_est1_1 <- cor(X1_1);
      S_est2_1 <- cor(X2_1);
      S_est1_2 <- cor(X1_2);
      S_est2_2 <- cor(X2_2);
      ##lower tri locs##
      lowerloc <- lower.tri(S_est1_1,diag = T)
      ##
      S_est1_1 <- S_est1_1[lowerloc];
      S_est2_1 <- S_est2_1[lowerloc];
      S_est1_2 <- S_est1_2[lowerloc];
      S_est2_2 <- S_est2_2[lowerloc];
      #
      beta1 <- S_est1_1-S_est2_1;
      beta2 <- S_est1_2-S_est2_2;
      ###
      if(is.null(Omega)){
        Omega<-Func_omega(X1_1,X2_1,Y1=NULL,Y2=NULL,model)
        Gamma = Sqrt(Omega)
      }else{
        Gamma<-Sqrt(Omega)
      }
    }else if(model=='multistage'){
      X1_1=X1[Index1,];
      X1_2=X1[-Index1,]
      X2_1=X2[Index2,];
      X2_2=X2[-Index2,]
      n1_1 <- nrow(X1_1)
      n1_2 <- nrow(X1_2) 
      n2_1 <- nrow(X2_1)
      n2_2 <- nrow(X2_2) 
      ##
      n_1 = n1_1+n2_1#nrow(dat_1);
      n_2 = n1_2+n2_2#nrow(dat_2)
      ##compute beta######
      if(is.null(beta)){
        beta1_1<-apply(X1_1,2,mean);
        beta2_1<-apply(X2_1,2,mean)
        beta1_2<-apply(X1_2,2,mean);
        beta2_2<-apply(X2_2,2,mean)
        beta1<-beta1_1-beta2_1;
        beta2<-beta1_2-beta2_2
        beta1 <- invGam%*%beta1
        beta2 <- invGam%*%beta2     
      }else{
        beta1=beta2=beta
      }
      Gamma <- Sqrt(Omega)
      }else if(model%in%c('lm','rlm2','glm')){
      n1_1<-length(Index1)
      n1_2<-n1-n1_1
      n2_1<-length(Index2)
      n2_2<-n2-n2_1
      ####  
      n_1<-n1_1+n2_1
      n_2<-n1_2+n2_2  
      #######compute beta, no need to compute S#######
      Beta1_1<-Beta1[Index1,]
      Beta1_2<-Beta1[-Index1,]
      Beta2_1<-Beta2[Index2,]
      Beta2_2<-Beta2[-Index2,]
      ######obtain average beta#####
      beta1_1<-apply(Beta1_1,2,mean)
      beta1_2<-apply(Beta1_2,2,mean)
      beta2_1<-apply(Beta2_1,2,mean)
      beta2_2<-apply(Beta2_2,2,mean)
      ######
      beta1<-as.vector(beta1_1-beta2_1)#diff of first half
      beta2<-as.vector(beta1_2-beta2_2)
      ######
      if(is.null(Omega)){
        S1_1<-S1[,,Index1]
        S2_1<-S2[,,Index2]
        s2<-s1<-matrix(0,p,p)
        for(i in 1:n1_1){
          s1<-s1+S1_1[,,i]
        }
        for(j in 1:n2_1){ 
          s2<-s2+S2_1[,,j]
        }
        s1<-s1/n1_1
        s2<-s2/n2_1
        Omega<-solve(1/n1_1*s1+1/n2_1*s2)
        Gamma = Sqrt(Omega)
      }else{
        Gamma=Sqrt(Omega)
      }
    }
    X<-Gamma
    Y<-as.vector(X%*%beta1)
    Yt<-as.vector(X%*%beta2)
    MS = Model_Select(Y, X)
    sv = MS$index_S
    w1 = MS$w1
    if ( length(sv)>0&&length(sv)<dim(X)[1]){
          bt2 = stats::lm(Yt~X[,sv]-1)$coefficients
          w2<-rep(0,dim(X)[1])                          #w2 = rep(0,p)
          w2[sv]=bt2
          sigma_w =rep(1,dim(X)[1])                     #rep(1, p)
          DIAG = diag( solve( t(X[,sv])%*%X[,sv] ) )
          sigma_w[sv] = DIAG
      if(scale){
        W = sqrt(n_1*n_2)*w1*w2/sqrt(sigma_w)##0.05鐨勬椂鍊欎繚瀹堢殑鎯呭喌鍙互瑙ｅ喅浜?##
      }else{
        W = sqrt(n_1*n_2)*w1*w2/sigma_w
      }
    }else{
      W = rep(0, dim(X)[1])
    }
    return(W)
  }
  # # FDR control based mirror
  ###alpha can be a vector c(0.05,0.1,0.2)###
  W_det <- function(Wj, alpha, options){
    t = sort(abs(Wj))
    
    if(options=='+'){
      Ta=sapply(t,function(x){(1+sum(Wj<=(-x)))/max(1,sum(Wj>=x))})
    }else{
      Ta=sapply(t,function(x){(0+sum(Wj<=(-x)))/max(1,sum(Wj>=x))})
    }
    #
    out<-do.call(rbind,lapply(alpha, function(x){
      bestlam = min(t[which(Ta<=x)])
      det=which(Wj>=bestlam)
      aa = rep(0, length(Wj))
      aa[det] = 1
      return(aa)
    }))
    return(out)
  }
  ####
  if(model%in%c('lm','rlm2','glm')){
    n1<-ncol(Y1)  
    n2<-ncol(Y2)  
    mp<-dim(X1)  
    m<-mp[1]
    p<-mp[2]
  }
  if(model%in%c('mean','multistage')){
    n1<-nrow(X1) 
    n2<-nrow(X2)
    p<-ncol(X1)
  }
  ######
  m1 = as.integer(2/3*n1)
  m2 = as.integer(2/3*n2)
  # stable option: single or multiple splitting
  if(stable){
    K=K
  }else{K=1}
  #
  Index1_K = lapply(1:K, function(x){sample(1:n1, m1)})
  Index2_K = lapply(1:K, function(x){sample(1:n2, m2)})
  #Wa is of dimension p*K, each column stands for one split result# 
  ####output the statistic Wa####K times#### 
  Wa = sapply(1:K, function(x){
    W_value(X1,X2,Y1,Y2,scale,beta,Beta1,Beta2,S1,S2,Index1_K[[x]],Index2_K[[x]],Omega,invGam,model,maxit)
  })
  deta = lapply(1:K, function(x){ W_det(Wa[, x],alpha,'-')}) #alpha[y], '-') })
  out<-list()
  stat<-list()
  for (i in 1:length(alpha)){
    tmpdeta<-as.matrix(do.call(rbind,lapply(1:K,function(x){deta[[x]][i,]})))
    #sum each column, corresponding to each split  
    mdet = apply(tmpdeta, 2, sum)
    det1 = which(mdet>=as.integer(0.5*K)+1)
    # the majority vote
    aa=rep(0, dim(Wa)[1])#dim(X)[1])#rep(0,p)
    if(length(det1)>0) aa[det1]=1
    # pick the best one
    k_hat = which.min(sapply(1:K,function(x){sum(abs(aa-tmpdeta[x,]))}))
    ##OUTPUT THE statistic##
    stat1=Wa[,k_hat]#The W value corresponding to k_hat
    stat[[i]]=stat1
    ##
    det1=(1:dim(Wa)[1])[tmpdeta[k_hat,]==1]
    if (length(det1)==0){
      out[[i]] = numeric(0)#'no rejection'
    }else{
      out[[i]] = det1
    }
  }
  Res<-list(out=out,stat=stat)
  return(Res)  
}
Func_LEB2 <-function(X1=NULL, X2=NULL, Y1=NULL,Y2=NULL, beta=NULL,Omega=NULL,invGam=NULL,model='mean',maxit=150){
  ##
  Model_Select <- function(Y, X,n,weight){
    fit <- glmnet::glmnet(X, Y, family = "gaussian",nlambda=500,penalty.factor=weight,intercept = F)
    Rss<-stats::deviance(fit)
    k<-fit$df
    c <-(log(n)+2*log(p))
    Err <- Rss+ c*k #Different constants for different criterion 
    i_min=which.min(Err)
    lambda_select = fit$lambda[i_min]
    fit<- glmnet::glmnet(X, Y, family = "gaussian", lambda = lambda_select, penalty.factor =weight,intercept = F)
    w1 = fit$beta[,1]
    sv = as.vector( which(w1!=0) )
    return( list(index_S=sv, w1= w1))
  }
  if(model=='mean'){
    n1<-nrow(X1);
    n2<-nrow(X2);
    #####beta###########
    #beta1<-apply(X1,2,mean);
    #beta2<-apply(X2,2,mean)
    #beta<-beta1-beta2
    if(is.null(beta)){
      beta1 <-apply(X1,2,mean)
      beta2 <-apply(X2,2,mean)
      beta <-beta1-beta2  
    }else{
      beta=beta
    }
    ######omega########
    if ( is.null(Omega) ){
      Omega<-Func_omega(X1,X2,Y1=NULL,Y2=NULL,model)
      Gamma = Sqrt(Omega)
    }else{Gamma<-Sqrt(Omega)}
    ##multistage is similar to high dimensional mean case##
  }else if(model=='multistage'){
    n1<-nrow(X1);
    n2<-nrow(X2);
    ###beta#######
    if(is.null(beta)){
      beta1 <-apply(X1,2,mean)
      beta2 <-apply(X2,2,mean)
      beta1 <-invGam%*%beta1
      beta2 <-invGam%*%beta2
      beta <-beta1-beta2  
    }else{
      beta=beta
    }
    Gamma<-Sqrt(Omega)#}
  }else if(model%in%c('lm','glm','rlm2')){
    n1=ncol(Y1)
    n2=ncol(Y2)
    ##########beta#########
    if(is.null(beta)|is.null(Omega)){
      mod1<-Func_betaS(X1,Y1,model,maxit = maxit)
      mod2<-Func_betaS(X1,Y2,model,maxit = maxit)
      beta1<-mod1$beta
      beta2<-mod2$beta
      beta<-beta1-beta2
      S1<-mod1$S
      S2<-mod2$S
      Omega <-solve(1/n1*S1+1/n2*S2)
      Gamma <-Sqrt(Omega)
    }else{
      beta<-beta
      Gamma<-Sqrt(Omega)
    }
  }
  ##############
  X<-Gamma
  Y<-as.vector(X%*%beta)
  ##effective sample size##
  E_n<-(n1*n2)/(n1+n2)###when n1 is large enough,the effective sample size is ~ n2 
  ##weights are the tilde beta
  MS = Model_Select(Y, X,n=E_n,weight=1/abs(beta))
  det1 = MS$index_S
  if (length(det1)==0){
    out = numeric(0)#'no rejection'
  }else{
    out = det1
  }
  return(out)
}
###LarInf###
Func_Lar <-function(X1=NULL, X2=NULL, Y1=NULL,Y2=NULL, beta=NULL,Omega=NULL,invGam=NULL,model='mean',maxit=150,alpha =alpha,Inf_type='lar'){
  ##Inf_type = 'lar' OR 'fs'##
  ##
  if(model=='mean'){
    n1<-nrow(X1);
    n2<-nrow(X2);
    #####beta###########
    beta1<-apply(X1,2,mean);
    beta2<-apply(X2,2,mean)
    beta<-beta1-beta2
    ######omega########
    if ( is.null(Omega) ){
      Omega<-Func_omega(X1,X2,Y1=NULL,Y2=NULL,model)
      Gamma = Sqrt(Omega)
    }else{Gamma<-Sqrt(Omega)}
    ##multistage is similar to high dimensional mean case##
  }else if(model=='multistage'){
    n1<-nrow(X1);
    n2<-nrow(X2);
    ###beta#######
    if(is.null(beta)){
      beta1 <-apply(X1,2,mean)
      beta2 <-apply(X2,2,mean)
      beta1 <-invGam%*%beta1
      beta2 <-invGam%*%beta2
      beta <-beta1-beta2  
    }else{
      beta=beta
    }
    Gamma<-Sqrt(Omega)#}
  }else if(model%in%c('lm','glm','rlm2')){
    n1=ncol(Y1)
    n2=ncol(Y2)
    ##########beta#########
    if(is.null(beta)|is.null(Omega)){
      mod1<-Func_betaS(X1,Y1,model,maxit = maxit)
      mod2<-Func_betaS(X1,Y2,model,maxit = maxit)
      beta1<-mod1$beta
      beta2<-mod2$beta
      beta<-beta1-beta2
      S1<-mod1$S
      S2<-mod2$S
      Omega <-solve(1/n1*S1+1/n2*S2)
      Gamma <-Sqrt(Omega)
    }else{
      beta<-beta
      Gamma<-Sqrt(Omega)
    }
  }
  ##############
  X <- Gamma
  Y <- as.vector(X%*%beta)
  p = dim(X)[1]
  #
  sigmahat <- estimateSigma(X,Y)$sigmahat
  if(Inf_type=='lar'){
    larfit <- lar(X,Y,maxsteps = ceiling(p/3))
    gc()
    out.aic <- larInf(larfit,type="aic",sigma = sigmahat)
    larfit_khat <- out.aic$vars#which(larfit$beta[,(out.aic$khat)]!=0)
    out <- lapply(alpha,function(x){sort(larfit_khat[which(p.adjust(out.aic$pv,method = 'BH')<x)])}) 
  }
  if(Inf_type=='fs'){
    fsfit <- fs(X,Y,maxsteps = ceiling(p/3))
    gc()
    out.aic <- fsInf(fsfit,type="aic",sigma = sigmahat)
    fsfit_khat <- out.aic$vars#which(fsfit$beta[,(out.aic$khat+1)]!=0)
    #out <- lapply(alpha,function(x){which(p.adjust(out.aic$pv,method = 'BH')<x)})  
    out <- lapply(alpha,function(x){sort(fsfit_khat[which(p.adjust(out.aic$pv,method = 'BH')<x)])}) 
  }
  return(out)
}
#####
Gen_Sigma<-function(p,rho=0.6,cov.type='ar'){
  if(cov.type=="ar"){
    sigma <- toeplitz(rho^(0:(p-1)))
  }
  if(cov.type=="bd"){             
    #(2) block diagonal matrix
    matr1<-diag(1,p,p);
    matr2<-cbind(rep(0,p-1),diag(c(rep(c(rho,0),p-2),rho),p-1,p-1))
    matr3<-rbind(matr2,rep(0,p))
    sigma<-matr1+matr3+t(matr3)
  }
  return(sigma)         
} 
##
Gen_highmean<-function(n1=1000,n2=1000,p=100,op=0.1,det1=0.5,det2=0.5,Sig1,OC='gaussian',OCmu='I'){
  Sig2=Sig1
  ##
  mu1=mu2=rep(0,p)
  if(OCmu=='I'){
    ##change only in mean
    Ind = sort(sample(1:p,ceiling(p*op))) #randomly choose mm occured change streams
    #random magnitudes with random signs#
    mu2[Ind]=mu1[Ind]+runif(length(Ind),det1,det2)*(2*rbinom(length(Ind),1,1/2)-1)
    THETA<-rep(0,p)
    THETA[Ind]=1
  }
  if(OCmu=='II'){
    sig12=0.05;
    THETA=rbinom(p,1,op)
    for(i in 1:p){
      if(THETA[i]==1){
        ind=rbinom(1,1,0.5)#rbinom(n,size,prob)
        mu2[i]=ind*rnorm(1,det1,sig12)+(1-ind)*rnorm(1,det2,sig12)
      }    
    }
    Ind<-which(THETA==1)
  }
  X1<-rmvnorm(n1,mu1,Sig1)
  if(OC=='gaussian'){
    X2<-rmvnorm(n2,mu2,Sig2)  
  }else if(OC=='t'){##df>1, the expectation is 0.## 
    X2<-rmvt(n2,Sig2,df=3)+rep(mu2,each=n2)
  }else if(OC=='gamma'){#expectation=shape*scale=2, so demean data by subtracting 2.
    X2 <- matrix(rgamma(n2*p,shape=2,scale=1),n2,p)-2
    eig <- eigen(Sig2)
    Sigsqrt <- eig$vectors%*%diag((eig$values)^0.5)%*%t(eig$vectors)
    X2 <- t(Sigsqrt%*%t(X2))
    X2<-X2+rep(mu2,each=n2)
  }
  return(list(X1=X1,X2=X2,Sig1=Sig1,Sig2=Sig2,mu1=mu1,mu2=mu2,Ind=Ind,THETA=THETA))
}
###
Gen_linearprofile<-function(X,n1=100,n2=50,m=500,p=100,rho=0.5,op=30,det1=0.5,det2=1,model='lm',cauchy_ratio=0.2){
    beta1 = rep(0,p)
    beta1[round(seq(1,p,length=10))] = 1##beta is sparse when p is large
    beta2=beta1
    Ind <- sort(sample(1:p,op))
    beta2[Ind]<-beta1[Ind]+runif(length(Ind),det1,det2)
  ####observations before and after change points#####
  if(model=='lm'){
      Y1<-t(do.call(rbind,lapply(1:n1,function(x){
        return(as.numeric(X %*% beta1 + rnorm(m)))  
      })))
      #Out of control phase#
      Y2<-t(do.call(rbind,lapply(1:n2,function(x){
        return(as.numeric(X %*% beta2 + rnorm(m)))
      })))
   }
  if(model%in%c('glm')){
    inv.logit<-function(x){return(exp(x)/(1+exp(x)))}
    
    Y1<-t(do.call(rbind,lapply(1:n1,function(x){
      
      pi <- as.vector(inv.logit(X%*%beta1))
      Y<-rbinom(m,1,pi)
      return(Y)
    })))
    Y2<-t(do.call(rbind,lapply(1:n2,function(x){
      pi <- as.vector(inv.logit(X%*%beta2))
      Y<-rbinom(m,1,pi)
      return(Y)
    })))
  }  
  if(model%in%c('rlm2')){
    Y1<-t(do.call(rbind,lapply(1:n1,function(x){
      cauchy_size<-ceiling(0.2*m)
      e1<-rcauchy(cauchy_size)#The error from cauchy error
      e2<-rnorm(m-cauchy_size)
      e<-c(e1,e2)
      Y<-as.numeric(X%*%beta1 + e )
      return(Y)
    })))
    Y2<-t(do.call(rbind,lapply(1:n2,function(x){
      cauchy_size<-ceiling(0.2*m)
      e1<-rcauchy(cauchy_size)#The error from cauchy error
      e2<-rnorm(m-cauchy_size)
      e<-c(e1,e2)
      Y<-as.numeric(X%*%beta2 + e )
      return(Y)
    })))
  }
  return(list(Y1=Y1,Y2=Y2,Ind=Ind))
}
########simplified multistage generation#########
Gen_multistage<-function(n1=1000,n2=1000,p=100,op=10,det1=0.5,det2=0.5,Ak=1,Ck=1,a0=0,sig_epsilon=1,sig_vk=1,sig_wk=1){
  delta=rep(0,p)##normal state
  Ind <- sort(sample(1:p,ceiling(op)))
  delta[Ind]=runif(length(Ind),det1,det2)*(2*rbinom(length(Ind),1,1/2)-1)
  V=matrix(rnorm((n1+n2)*p,0,sqrt(sig_vk)),nrow=p)#error term of Y
  W=matrix(rnorm((n1+n2)*p,0,sqrt(sig_wk)),nrow=p)#error term of X
  X0=rnorm((n1+n2),a0,sqrt(sig_epsilon))
  ###
  X=matrix(0,p,(n1+n2))##The first row means the initial state  
  Y=X
  ##
  X[1,(1:n1)]=Ak*X0[(1:n1)]+W[1,(1:n1)]##in-control sample
  X[1,(n1+1):(n1+n2)]=Ak*X0[(n1+1):(n1+n2)]+W[1,(n1+1):(n1+n2)]+delta[1]*ifelse(1%in%Ind,1,0)###杩欒竟娌¤繖涓猨鍝?
  Y[1,]=Ck*X[1,]+V[1,]
  for(k in 2:p){
    X[k,(1:n1)]=Ak*X[(k-1),(1:n1)]+W[k,(1:n1)]##in-control sample
    X[k,(n1+1):(n1+n2)]=Ak*X[(k-1),(n1+1):(n1+n2)]+W[k,(n1+1):(n1+n2)]+delta[k]*ifelse(k%in%Ind,1,0)###杩欒竟娌¤繖涓猨鍝?
  }
  Y[-1,]=Ck*X[-1,]+V[-1,]  
  ##
  W =V=G= rep(0,p)
  #initial of W
  W[1]=sig_wk+Ak^2*sig_epsilon
  for(k in 2:p){
    W[k]=Ak^2*W[k-1]-Ak*Ck*W[k-1]*G[k-1]+sig_wk  
  }
  ##
  V=Ck^2*W+sig_vk
  G=Ak*Ck*W/V
  e=matrix(0,p,(n1+n2))#e_nj denotes the j-th product of the n-th stage 
  v=u=e
  v[1,]=Y[1,]-Ck*Ak*a0
  e[1,]=v[1,]/V[1]
  for(k in 2:p){
    u[k,]=Ak*u[(k-1),]+G[(k-1)]*v[(k-1),]
    v[k,]=Y[k,]-Ck*u[k,]
    e[k,]=v[k,]*(V[k])^{-1/2}#return e
  }
  X1 <- t(Y[,1:n1])#in-control sample of multistage process
  X2 <- t(Y[,(n1+1):(n1+n2)])#out-of-control sample of multistage process
  e1 <- t(e[,1:n1])#one-step standardized residual of in-control sample
  e2 <- t(e[,(n1+1):(n1+n2)])#residual of out-of-control sample
  ###delta is the true beta difference###
  return(list(X1=X1,e1=e1,X2=X2,e2=e2,Ind=Ind,delta=delta))
}  
###
Gen_multicov<-function(p=100,Ak=1,Ck=1,a0=0,sig_epsilon=1,sig_vk=1,sig_wk=1){
  Gam <- matrix(0,p,p)
  for(j in 1:(p-1)){
    for(i in (j+1):p){
      Gam[i,j]=Ck*(Ak)^{i-j}#prod(A[(j+1):i])
    }
  }
  diag(Gam)=rep(Ck,p)
  invGam <- solve(Gam)
  VarY<-rep(0,p)
  VarY[1]=Ck^2*(Ak^2*sig_epsilon+sig_wk)+sig_vk
  for(k in 2:p){
    #k=p
    sum_prodA<-0
    for(l in 1:(k-1)){
      #l=1
      tmp<-sig_wk*(Ak^{k-l})^2        #(prod(A[(l+1):(k)])^2)
      sum_prodA<-sum_prodA+tmp  
    }
    VarY[k]=Ck^2*(Ak^{2*k}*sig_epsilon+sum_prodA)+sig_vk+Ck^2*sig_wk
  }
  covY<-matrix(0,p,p)
  diag(covY)=VarY
  ###triple loop,very slow###
  for(l in 2:(p-1)){
    for(k in (l+1):p){#l<k
      ##
      sum_prodA<-0
      for(m in 1:(l-1)){
        tmp<-sum(Ck*Ak^{k-m}*Ck*Ak^{l-m}*sig_wk)
        sum_prodA<-sum_prodA+tmp
      }
      covY[l,k]=Ck*Ak^{k}*Ck*Ak^{l}*sig_epsilon+sum_prodA+Ck*Ak^{k-l}*Ck*sig_wk
    }
  }
  ###
  for(k in 2:p){#l<k
    covY[1,k]=Ck*Ak^{k}*Ck*Ak^{1}*sig_epsilon+
      Ck*Ak^{k-1}*Ck*sig_wk
  }
  ###make covY symmetric### 
  covY=covY+t(covY)-diag(diag(covY))####涔嬪墠鏄眰浜哻ovY鐨勪竴鍗婏紝鍐嶅姞涓婂彟涓€鍗婏紝澶氬姞浜嗗瑙掔嚎锛屽啀鍓帀瀵硅绾?###
  #covY=(covY+t(covY))/2
  covYY<-invGam%*%covY%*%t(invGam)
  covYY=(covYY+t(covYY))/2#make it symmetric
  if(!(matrixcalc::is.positive.definite(covYY))){
    covYY<-as.matrix(Matrix::nearPD(covYY, conv.tol = 1e-15, conv.norm = "F")$mat)
  }
  return(list(Gam=Gam,invGam=invGam,covY=covY,covYY=covYY))
}  
#=====#
OGA <- function(X, y, Kn = NULL, c1 = 5){
  if (!is.vector(y)) stop("y should be a vector")
  if (!is.matrix(X)) stop("X should be a matrix")
  
  n = nrow(X)
  p = ncol(X)
  if (n != length(y)) stop("the number of observations in y is not equal to the number of rows of X")
  if (n == 1) stop("the sample size should be greater than 1")
  if (is.null(Kn)) K = max(1, min(floor(c1 * sqrt(n / log(p))), p))
  else{
    if ((Kn < 1) | (Kn > p)) stop(paste("Kn should between 1 and ", p, sep = ""))
    if ((Kn - floor(Kn)) != 0) stop("Kn should be a positive integer")
    K = Kn
  }
  
  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))
  
  Jhat = sigma2hat = rep(0, K)
  XJhat = matrix(0, n, K)
  u = as.matrix(dy)
  xnorms = sqrt(colSums((dX) ^ 2))
  
  aSSE = (abs(t(u) %*% dX) / xnorms)
  Jhat[1] = which.max(aSSE)
  XJhat[, 1] = (dX[, Jhat[1]] / sqrt(sum((dX[, Jhat[1]]) ^ 2)))
  u = u - XJhat[, 1] %*% t(XJhat[, 1]) %*% u
  sigma2hat[1] = mean(u ^ 2)
  
  if (K > 1){
    for (k in 2:K) {
      aSSE = (abs(t(u) %*% dX) / xnorms)
      aSSE[Jhat[1:(k-1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dX[, Jhat[k]] - XJhat[, 1:(k-1)] %*% t(XJhat[, 1:(k-1)]) %*% dX[, Jhat[k]]
      XJhat[, k] = (rq / sqrt(sum((rq) ^ 2)))
      u = u - XJhat[, k] %*% t(XJhat[, k]) %*% u
      sigma2hat[k] = mean(u ^ 2)
    }
  }
  
  return(list("n" = n, "p" = p, "Kn" = K, "J_OGA" = Jhat))
}
#===========#
Ohit <- function(X, y, Kn = NULL, c1 = 5, HDIC_Type = "HDBIC", c2 = 2, c3 = 2.01, intercept = TRUE){
  if (!is.vector(y)) stop("y should be a vector")
  if (!is.matrix(X)) stop("X should be a matrix")
  
  n = nrow(X)
  p = ncol(X)
  if (n != length(y)) stop("the number of observations in y is not equal to the number of rows of X")
  if (n == 1) stop("the sample size should be greater than 1")
  if (is.null(Kn)) K = max(1, min(floor(c1 * sqrt(n / log(p))), p))
  else{
    if ((Kn < 1) | (Kn > p)) stop(paste("Kn should between 1 and ", p, sep = ""))
    if ((Kn - floor(Kn)) != 0) stop("Kn should be a positive integer")
    K = Kn
  }
  
  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))
  
  Jhat = sigma2hat = rep(0, K)
  XJhat = matrix(0, n, K)
  u = as.matrix(dy)
  xnorms = sqrt(colSums((dX) ^ 2))
  
  aSSE = (abs(t(u) %*% dX) / xnorms)
  Jhat[1] = which.max(aSSE)
  XJhat[, 1] = (dX[, Jhat[1]] / sqrt(sum((dX[, Jhat[1]]) ^ 2)))
  u = u - XJhat[, 1] %*% t(XJhat[, 1]) %*% u
  sigma2hat[1] = mean(u ^ 2)
  
  if (K > 1){
    for (k in 2:K) {
      aSSE = (abs(t(u) %*% dX) / xnorms)
      aSSE[Jhat[1:(k-1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dX[, Jhat[k]] - XJhat[, 1:(k-1)] %*% t(XJhat[, 1:(k-1)]) %*% dX[, Jhat[k]]
      XJhat[, k] = (rq / sqrt(sum((rq) ^ 2)))
      u = u - XJhat[, k] %*% t(XJhat[, k]) %*% u
      sigma2hat[k] = mean(u ^ 2)
    }
  }
  
  if ((HDIC_Type != "HDAIC") & (HDIC_Type != "HDBIC") & (HDIC_Type != "HDHQ")) stop("HDIC_Type should be \"HDAIC\", \"HDBIC\" or \"HDHQ\"")
  if (HDIC_Type == "HDAIC") omega_n = c2
  if (HDIC_Type == "HDBIC") omega_n = log(n)
  if (HDIC_Type == "HDHQ") omega_n = c3 * log(log(n))
  
  hdic = (n * log(sigma2hat))+((1:K) * omega_n * (log(p)))
  kn_hat = which.min(hdic)
  benchmark = hdic[kn_hat]
  J_HDIC = sort(Jhat[1:kn_hat])
  
  J_Trim = Jhat[1:kn_hat]
  trim_pos = rep(0, kn_hat)
  if (kn_hat > 1){
    for (l in 1:(kn_hat-1)){
      JDrop1 = J_Trim[-l]
      fit = lm(dy~.-1, data = data.frame(dX[, JDrop1]))
      uDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(uDrop1 ^ 2)))+((kn_hat - 1) * omega_n * (log(p)))
      if (HDICDrop1 > benchmark) trim_pos[l] = 1
    }
    trim_pos[kn_hat] = 1
    J_Trim = J_Trim[which(trim_pos==1)]
  }
  J_Trim = sort(J_Trim)
  
  X_HDIC = as.data.frame(as.matrix(X[, J_HDIC]))
  X_Trim = as.data.frame(as.matrix(X[, J_Trim]))
  X = data.frame(X)
  colnames(X_HDIC) = names(X)[J_HDIC]
  colnames(X_Trim) = names(X)[J_Trim]
  
  if (intercept == TRUE){
    fit_HDIC = lm(y~., data = X_HDIC)
    fit_Trim = lm(y~., data = X_Trim)
  }else{
    fit_HDIC = lm(y~.-1, data = X_HDIC)
    fit_Trim = lm(y~.-1, data = X_Trim)
  }
  
  betahat_HDIC = summary(fit_HDIC)
  betahat_Trim = summary(fit_Trim)
  
  return(list("n" = n, "p" = p, "Kn" = K, "J_OGA" = Jhat, "HDIC" = hdic, "J_HDIC" = J_HDIC, "J_Trim" = J_Trim, "betahat_HDIC" = betahat_HDIC, "betahat_Trim" = betahat_Trim))
}
#================#
predict_Ohit <- function(object, newX){
  if (!is.matrix(newX)) stop("newX should be a matrix")
  if (ncol(newX) != object$p) stop(paste("the number of columns of newX is not equal to ", object$p, sep = ""))
  
  if (length(object$J_HDIC) == nrow(object$betahat_HDIC$coefficients)){
    pred_HDIC = newX[, object$J_HDIC] %*% as.matrix(object$betahat_HDIC$coefficients[, 1])
  }else if (nrow(newX) == 1){
    pred_HDIC = sum(c(1, newX[, object$J_HDIC]) * object$betahat_HDIC$coefficients[, 1])
  }else{
    pred_HDIC = cbind(rep(1, nrow(newX)), newX[, object$J_HDIC]) %*% as.matrix(object$betahat_HDIC$coefficients[, 1])
  }
  
  if (length(object$J_Trim) == nrow(object$betahat_Trim$coefficients)){
    pred_Trim = newX[, object$J_Trim] %*% as.matrix(object$betahat_Trim$coefficients[, 1])
  }else if (nrow(newX) == 1){
    pred_Trim = sum(c(1, newX[, object$J_Trim]) * object$betahat_Trim$coefficients[, 1])
  }else{
    pred_Trim = cbind(rep(1, nrow(newX)), newX[, object$J_Trim]) %*% as.matrix(object$betahat_Trim$coefficients[, 1])
  }
  
  return(list("pred_HDIC" = as.vector(pred_HDIC), "pred_Trim" = as.vector(pred_Trim)))
}
##=========##
sigmarcv <-
  function(y,X,cv=FALSE,fit=NA,intercept=TRUE){
    
    #Taken from:
    
    #Variance estimation using refitted cross-validation in
    #ultrahigh dimensional regression
    
    #Jianqing Fan,
    #Princeton University, USA
    #Shaojun Guo
    #Chinese Academy of Sciences, Beijing, PRC
    #and Ning Hao
    #University of Arizona, Tucson, USA
    
    #Journal of Royal Statistics Society
    #2011
    
    
    #kindly proportioned by the paper authors
    
    
    n=nrow(X)
    p=ncol(X)
    
    # we assume the data have been resampled by random permutation. If not, resample the label here.
    
    k          <- floor(n/2)
    x1         <- X[1:k,]
    y1         <- y[1:k]
    x2         <- X[(k+1):n,]
    y2         <- y[(k+1):n]
    
    out=NULL
    
    if(cv==TRUE){
      #CV
      flasso=fit
      sigma.est1=NA
      maxiter=0
      while(is.na(sigma.est1)|sigma.est1==0){
        if(is.na(fit)) flasso     <- lars(X,y,type="lasso",eps=1e-10,use.Gram=FALSE,intercept=intercept);
        object     <- cv.lars(X,y,type="lasso",use.Gram=FALSE,plot.it=FALSE,se=FALSE)
        beta.est   <- predict.lars(flasso,s=object$index[which.min(object$cv)],
                                   type="coefficients",mode="fraction")$coef
        
        df         <- length(which(abs(beta.est)>0))
        if(intercept) beta0=mean(y)-mean(X%*%beta.est) else beta0=0
        
        sigma.est1 <- sum((y-X%*%beta.est-beta0)^2)/(n-df)
        out$sigmahat=sqrt(sigma.est1)
        out$ind=which(beta.est!=0)
        maxiter=maxiter+1
        if(maxiter==10) stop("The maximum number of iterations where reached and all cross validation models were full rank. No estimation of sigma was able to be done by sigma='cv'.  Try sigma='qut'.")
      }
    }
    else if(cv==FALSE){
      #RCV
      maxiter=0
      sigma.est3=NA
      sigma.est4=NA
      while(is.na(sigma.est3)|is.na(sigma.est4)|(sigma.est3==0&sigma.est4==0)){
        flasso     <- lars(x1,y1,type="lasso",eps=1e-10,use.Gram=FALSE,intercept=intercept);
        object     <- cv.lars(x1,y1,type="lasso",use.Gram=FALSE,plot.it=FALSE,se=FALSE)
        beta.est1  <- predict.lars(flasso,s=object$index[which.min(object$cv)],
                                   type="coefficients",mode="fraction")$coef
        ind1       <- which(abs(beta.est1)>0)
        if ( length(ind1) == 0 ) {
          sigma.est3 <- var(y2)
        } else {
          if(intercept) object     <- lm(y2 ~ x2[,ind1] + 1) else object     <- lm(y2 ~ x2[,ind1] - 1)
          sigma.est3 <- sum((object$resid)^2)/(n-k-length(ind1))
        }
        
        flasso     <- lars(x2,y2,type="lasso",eps=1e-10,use.Gram=FALSE,intercept=intercept);
        object     <- cv.lars(x2,y2,type="lasso",use.Gram=FALSE,plot.it=FALSE,se=FALSE)
        beta.est2  <- predict.lars(flasso,s=object$index[which.min(object$cv)],
                                   type="coefficients",mode="fraction")$coef
        ind2       <- which(abs(beta.est2)>0)
        if ( length(ind2) ==0 ) {
          sigma.est4 <- var(y1)
        } else {
          if(intercept) object     <- lm(y1 ~ x1[,ind2] + 1) else object     <- lm(y1 ~ x1[,ind2] - 1) 
          sigma.est4 <- sum((object$resid)^2)/(k-length(ind2))
        }
        maxiter=maxiter+1
        if(maxiter==10) stop("The maximum number of iterations where reached and all RCV models were full rank. No estimation of sigma was able to be done by sigma='rcv'.  Try sigma='qut'.")
      }
      out$sigmahat=sqrt((sigma.est3+sigma.est4)/2)
      out$ind=which(abs(beta.est1*beta.est2)>0)
    }
    
    return(out)
  }
