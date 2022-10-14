###### DA analysis of the CDI microbiome study 
load("~/Desktop/VincentC.phyloseq.obj.RData")
NN=1000
bugs<-c(bugs.kingdom,bugs.phylum,bugs.class,bugs.order,bugs.family,bugs.genus,bugs.species.and.strains)
p_value=array(1:4*7*NN,c(4,7,NN))
for (rank in 1:7) {
  map = bugs[rank][[1]]@sam_data ## the meta data
  otu.tab = bugs[rank][[1]]@otu_table
  ind_control=map$study_condition =='control' & 0==map$days_from_first_collection 
  # ind_control= map$study_condition =='control' & 0<map$days_from_first_collection &map$days_from_first_collection<8
  # ind_control= map$study_condition =='control' & map$days_from_first_collection >7
  ind_uniq_control=duplicated(map$subjectID[ind_control])
  sample_name_control=rownames(map)[ind_control][!ind_uniq_control]
  for (N in 1:NN) {
    sample_ind=sample(seq(1,length(sample_name_control),1),floor(length(sample_name_control)/2),replace=FALSE)
    sample_name_1=sample_name_control[sample_ind]
    sample_name_2=sample_name_control[-sample_ind]
    
    phy_control = subset_samples(bugs[rank][[1]],ind_control)
    
    
    Z_1=otu_table(phy_control)[,sample_name_1]
    Z_2=otu_table(phy_control)[,sample_name_2]
    
    n_1=ncol(Z_1)
    n_2=ncol(Z_2)
    p=nrow(Z_1)
    U=diag(1,p,p)-(1/p)*matrix(1,nrow = p,ncol = 1)%*%matrix(1,ncol = p,nrow  = 1)
    
    
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
    
    L_1_bar=rowSums(log(X_1))/n_1
    L_2_bar=rowSums(log(X_2))/n_2
    
    X_1_bar=rowSums(X_1)/n_1
    X_2_bar=rowSums(X_2)/n_2
    
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Z_c[,mm]))/(n1)+(var(Z_a[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)
    
    Tssmm=(max((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma)))^2+(min((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma)))^2
    tssmm=(Tssmm-2*c)/(2-((log(p-1))^(-1)))
    f <- function(t) exp(exp(-tssmm)/log(t))
    p_value[1,rank,N]=1- integrate( f,0,1)$value
    
    
    
    
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(Y_1[,mm]))/(n1)+(var(Y_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((Y_1_bar[1:(p-1)]-Y_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value[2,rank,N]=1-exp(-exp(-tssmm))
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(L_1[,mm]))/(n1)+(var(L_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((L_1_bar[1:(p-1)]-L_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value[3,rank,N]=1-exp(-exp(-tssmm))
    
    gamma=rep(0,(p-1))
    for (mm in 1:(p-1)) {
      gamma[mm]=(var(X_1[,mm]))/(n1)+(var(X_2[,mm]))/(n2)
    }
    c=2*log(p-1)-(log(log(p-1))+log(4*pi))+(log(log(p-1))+log(4*pi))*(2*log(p-1))^(-1)+log(4)-log(4)/(2*log(p-1))
    Tssmm=max(((X_1_bar[1:(p-1)]-X_2_bar[1:(p-1)])/sqrt(gamma))^2)
    tssmm=(Tssmm-c)/(2-((log(p-1))^(-1)))
    p_value[4,rank,N]=1-exp(-exp(-tssmm))
  }
  
  
}

error=matrix(0,nrow=4,ncol=7)

for (M1 in 1:4) {
  for (M2 in 1:7) {
    error[M1,M2] =length(which(p_value[M1,M2,]<0.05))/NN
  }
}
