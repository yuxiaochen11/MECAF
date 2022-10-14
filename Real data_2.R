## procanova
load("~/Desktop/VincentC.phyloseq.obj.RData")
NN=1000
bugs<-c(bugs.kingdom,bugs.phylum,bugs.class,bugs.order,bugs.family,bugs.genus,bugs.species.and.strains)
p_v=array(1:1*7*NN,c(1,7,NN))
for (rank in 1:7) {
  map = bugs[rank][[1]]@sam_data ## the meta data
  otu.tab = bugs[rank][[1]]@otu_table
  
  ind_control=map$study_condition =='control' & 0==map$days_from_first_collection 
  
  ind_uniq_control=duplicated(map$subjectID[ind_control])
  sample_name_control=rownames(map)[ind_control][!ind_uniq_control]
  
  for (N in 1:NN) {
    sample_ind=sample(seq(1,length(sample_name_control),1),floor(length(sample_name_control)/2),replace=FALSE)
    sample_name_1=sample_name_control[sample_ind]
    sample_name_2=sample_name_control[-sample_ind]
    
    phy_control = subset_samples(bugs[rank][[1]],ind_control)
    
    
    Z_1=otu_table(phy_control)[,sample_name_1]
    Z_2=otu_table(phy_control)[,sample_name_2]
    
    n1=ncol(Z_1)
    n2=ncol(Z_2)
    p=nrow(Z_1)
    
    
    
    Z_1[Z_1==0]<-0.5
    Z_2[Z_2==0]<-0.5
    library(vegan)
    otu=abs(data.frame(t(cbind(Z_1,Z_2))))
    
    a=rep(0,n1)
    for (i in 1:n1) {
      a[i]=paste('A',i,seq='')
    }
    
    rownames(otu)[1:n1]<-a
    
    b=rep(0,n2)
    for (i in 1:n2) {
      b[i]=paste('B',i,seq='')
    }
    
    rownames(otu)[(n1+1):(n1+n2)]<-b
    group=data.frame(rownames(otu))
    group$site=c(rep('A',n1),rep('B',n2))
    
    
    r<-adonis2(otu~site,group,distance='bray',permutations=999)  
    p_v[1,rank,N]=r$`Pr(>F)`[1]
    
  }
  
  
}

error=matrix(0,nrow=1,ncol=7)

for (M1 in 1:1) {
  for (M2 in 1:7) {
    error[M1,M2]=length(which(p_v[M1,M2,]<0.05))/NN
  }
}
