###type 1 error and power : same banded covariance matrix between compared groups
type1error_equal_banded<-function(p,n1,n2,alpha,N)
  
{
  p_value_TMECAF=rep(0,N)
  p_value_TMEC=rep(0,N)
  p_value_oracle=rep(0,N)
  p_value_log=rep(0,N)
  p_value_raw=rep(0,N)
  for (j in 1:N) {
    mu1=runif(n=p,min=1,max=10)
    un=runif(n=p,min=1,max=3)
    D=diag(x=un,nrow = p,ncol=p)
    library(pracma)
    D_0.5=sqrtm(D)$B
    A=diag(1,p,p)
    diag(A[-nrow(A),-1])<-(-0.5)
    A=A+t(A)-diag(diag(A))
    omega1=D_0.5%*%A%*%D_0.5
    library(mvtnorm)
    Z_1=rmvnorm(n=n1,mean = mu1,sigma=omega1)
    W_1=exp(Z_1)
    X_1=matrix(0,n1,p)
    for (i in 1:n1) { 
      X_1[i,]=W_1[i,]/sum(W_1[i,])
    }
    Y_1=matrix(0,n1,p)
    for (i in 1:n1) { 
      Y_1[i,]=(diag(1,p,p)-(1/p)*matrix(1,nrow = p,ncol = 1)%*%matrix(1,ncol = p,nrow  = 1))%*%log(X_1[i,])
    }
    
    omega2<-omega1
    mu2=mu1
    Z_2=rmvnorm(n=n2,mean = mu2,sigma=omega2)
    W_2=exp(Z_2)
    X_2=matrix(0,n2,p)
    for (i in 1:n2) { 
      X_2[i,]=W_2[i,]/sum(W_2[i,])
    }
    Y_2=matrix(0,n2,p)
    for (i in 1:n2) { 
      Y_2[i,]=(diag(1,p,p)-(1/p)*matrix(1,nrow = p,ncol = 1)%*%matrix(1,ncol = p,nrow  = 1))%*%log(X_2[i,])
    }
    
    
    Y_1_bar=colSums(Y_1)/n1
    Y_2_bar=colSums(Y_2)/n2
    
    L_1_bar=colSums(log(X_1))/n1
    L_2_bar=colSums(log(X_2))/n2
    
    X_1_bar=colSums(X_1)/n1
    X_2_bar=colSums(X_2)/n2
    
    Z_1_bar=colSums(Z_1)/n1
    Z_2_bar=colSums(Z_2)/n2
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Z_c[,mm]))/(n1)+(var(Z_a[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)
    
    Tssmm=(max((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma)))^2+(min((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma)))^2
    tssmm=(Tssmm-2*c)/(2-((log(p-1))^(-1)))
    f <- function(t) exp(exp(-tssmm)/log(t))
    p_value_MECAF[j]=1- integrate( f,0,1)$value
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Y_1[,mm]))/(n1)+(var(Y_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_TMEC[j]=1-exp(-exp(-tssmm))
    
    #oracle,log,raw
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Z_1[,mm]))/(n1)+(var(Z_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((Z_1_bar[1:(p-1)]-Z_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_oracle[j]=1-exp(-exp(-tssmm))
    
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(L_1[,mm]))/(n1)+(var(L_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((L_1_bar[1:(p-1)]-L_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_log[j]=1-exp(-exp(-tssmm))
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(X_1[,mm]))/(n1)+(var(X_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((X_1_bar[1:(p-1)]-X_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_raw[j]=1-exp(-exp(-tssmm))
  }
  
  errorTMECAF=length(which(p_value_TMECAF<alpha))/N
  errorTMEC=length(which(p_value_TMEC<alpha))/N
  oracle=length(which(p_value_oracle<alpha))/N
  log=length(which(p_value_log<alpha))/N
  raw=length(which(p_value_raw<alpha))/N
  
  return(c(errorTMECAF,errorTMEC,oracle,log,raw))
}





power_equal_banded<-function(p,n1,n2,alpha,N,ss)
  
{
  p_value_TMECAF=rep(0,N)
  p_value_TMEC=rep(0,N)
  p_value_oracle=rep(0,N)
  p_value_log=rep(0,N)
  p_value_raw=rep(0,N)
  for (j in 1:N) {
    mu1=runif(n=p,min=1,max=10)
    un=runif(n=p,min=1,max=3)
    D=diag(x=un,nrow = p,ncol=p)
    library(pracma)
    D_0.5=sqrtm(D)$B
    A=diag(1,p,p)
    diag(A[-nrow(A),-1])<-(-0.5)
    A=A+t(A)-diag(diag(A))
    omega1=D_0.5%*%A%*%D_0.5
    library(mvtnorm)
    Z_1=rmvnorm(n=n1,mean = mu1,sigma=omega1)
    W_1=exp(Z_1)
    X_1=matrix(0,n1,p)
    for (i in 1:n1) { 
      X_1[i,]=W_1[i,]/sum(W_1[i,])
    }
    Y_1=matrix(0,n1,p)
    for (i in 1:n1) { 
      Y_1[i,]=(diag(1,p,p)-(1/p)*matrix(1,nrow = p,ncol = 1)%*%matrix(1,ncol = p,nrow  = 1))%*%log(X_1[i,])
    }
    
    omega2<-omega1
    mu2=mu1
    Z_2=rmvnorm(n=n2,mean = mu2,sigma=omega2)
    W_2=exp(Z_2)
    X_2=matrix(0,n2,p)
    for (i in 1:n2) { 
      X_2[i,]=W_2[i,]/sum(W_2[i,])
    }
    Y_2=matrix(0,n2,p)
    for (i in 1:n2) { 
      Y_2[i,]=(diag(1,p,p)-(1/p)*matrix(1,nrow = p,ncol = 1)%*%matrix(1,ncol = p,nrow  = 1))%*%log(X_2[i,])
    }
    
    
    Y_1_bar=colSums(Y_1)/n1
    Y_2_bar=colSums(Y_2)/n2
    
    L_1_bar=colSums(log(X_1))/n1
    L_2_bar=colSums(log(X_2))/n2
    
    X_1_bar=colSums(X_1)/n1
    X_2_bar=colSums(X_2)/n2
    
    Z_1_bar=colSums(Z_1)/n1
    Z_2_bar=colSums(Z_2)/n2
    
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Z_c[,mm]))/(n1)+(var(Z_a[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)
    
    Tssmm=(max((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma)))^2+(min((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma)))^2
    tssmm=(Tssmm-2*c)/(2-((log(p-1))^(-1)))
    f <- function(t) exp(exp(-tssmm)/log(t))
    p_value_MECAF[j]=1- integrate( f,0,1)$value
    
    
    
    
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Y_1[,mm]))/(n1)+(var(Y_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_TMEC[j]=1-exp(-exp(-tssmm))
    
    #oracle,log,raw
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Z_1[,mm]))/(n1)+(var(Z_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((Z_1_bar[1:(p-1)]-Z_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_oracle[j]=1-exp(-exp(-tssmm))
    
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(L_1[,mm]))/(n1)+(var(L_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((L_1_bar[1:(p-1)]-L_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_log[j]=1-exp(-exp(-tssmm))
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(X_1[,mm]))/(n1)+(var(X_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((X_1_bar[1:(p-1)]-X_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value_raw[j]=1-exp(-exp(-tssmm))
  }
  
  powerTMECAF=length(which(p_value_TMECAF<alpha))/N
  powerTMEC=length(which(p_value_TMEC<alpha))/N
  oracle=length(which(p_value_oracle<alpha))/N
  log=length(which(p_value_log<alpha))/N
  raw=length(which(p_value_raw<alpha))/N
  
  return(c(powerTMECAF,powerTMEC,oracle,log,raw))
}
