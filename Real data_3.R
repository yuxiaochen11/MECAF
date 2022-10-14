#####DA analysis in female and male mice of the murine microbiome study
library(phyloseq)
load("~/Desktop/repeated.measure.gut.RData")#phy: a phyloseq object
N=1000
n_1=20
n_2=10
alpha=0.05
tax_glom(phy,taxrank= rank_names(phy)[2])  -> phy.type
map = sample_data(phy.type) ## the meta data
otu.tab = otu_table(phy.type)
ind_1=map$Week== 13 & map$Treatment =='STAT' & map$Sex == 'F'
ind_2=map$Week== 13 & map$Treatment =='Control'  & map$Sex =='F'

sample_name_1=rownames(map)[ind_1]
sample_name_2=rownames(map)[ind_2]

phy_1 = subset_samples(phy.type,ind_1)
phy_2 = subset_samples(phy.type,ind_2)

p_value_TMECAF=rep(0,N)
p_value_TMEC=rep(0,N)
p_value_oracle=rep(0,N)
p_value_log=rep(0,N)
p_value_raw=rep(0,N)

for (j in 1:N) {
  
  Z_1=otu_table(phy_1)[,sample_name_1]
  Z_2=otu_table(phy_2)[,sample_name_2]
  
  
  
  N_1=length(rownames(map)[ind_1])
  N_2=length(rownames(map)[ind_2])
  p=length(rownames(Z_2))
  U=diag(1,p,p)-(1/p)*matrix(1,nrow = p,ncol = 1)%*%matrix(1,ncol = p,nrow  = 1)
  
  
  c1=sample(1:N_1,n_1,replace = F)
  c2=sample(1:N_2,n_2,replace = F)
  
  Z_1=Z_1[,c1]
  Z_2=Z_2[,c2]
  
  Z_1[Z_1==0]<-0.5
  Z_2[Z_2==0]<-0.5
  
  
  X_1=matrix(0,p,n_1)
  for (i in 1:n_1) { 
    X_1[,i]=Z_1[,i]/sum(Z_1[,i])
  }
  
  Y_1=matrix(0,p,n_1)
  for (i in 1:n_1) { 
    Y_1[,i]=U%*%log(X_1[,i])
  }
  
  
  
  X_2=matrix(0,p,n_2)
  for (i in 1:n_2) { 
    X_2[,i]=Z_2[,i]/sum(Z_2[,i])
  }
  Y_2=matrix(0,p,n_2)
  for (i in 1:n_2) { 
    Y_2[,i]=U%*%log(X_2[,i])
  }
  
  Y_1_bar=rowSums(Y_1)/n_1
  Y_2_bar=rowSums(Y_2)/n_2
  
  L_1_bar=rowSums(log(X_1))/n1
  L_2_bar=rowSums(log(X_2))/n2
  
  X_1_bar=rowSums(X_1)/n1
  X_2_bar=rowSums(X_2)/n2
  
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
log=length(which(p_value_log<alpha))/N
raw=length(which(p_value_raw<alpha))/N

