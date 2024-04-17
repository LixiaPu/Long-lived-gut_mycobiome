################################################################################
#main codes for enterotype analysis
################################################################################
library(cluster)
library(clusterSim)
library(vegan)
library(ade4)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsci)
library(scales)

genus.its<- read.csv("Long-lived-gut_mycobiome/level-genus.df.csv",sep=",",header=T,strip.white = T, stringsAsFactors = F,check.names=F)
metadata <- read.csv("Long-lived-gut_mycobiome_project/metadata.csv", header=TRUE, sep = ",",row.names = "SampleID",stringsAsFactors = F,fileEncoding = "GBK")
df<-merge(metadata,genus.its,by = "SampleID")
row.names(df)<-df$SampleID
data.df<-data.df%>%dplyr::select(!(1:92))
data.df<-t(data.df)
data.genus<-apply(data.df,2,as.numeric)
row.names(data.genus)<-row.names(data.genus)
data.genus<-apply(data.genus,2,function(x) x/sum(x)*100)
data.denoized<-data.genus[apply(data.genus,1,function(x) sum(x>0)>=0.2*length(x)),]
#PAM cluster +bary
data.dist.bc <- vegdist(t(data.denoized),method = "bray")
write.table(as.matrix(data.dist.bc),"bray_curtis-genus.txt",sep = '\t',quote = FALSE,col.names = NA)
#Calinski-Harabasz 
nclusters.bc = index.G1(t(data.denoized), data.cluster, d = data.dist.bc, enterotypes = "medoids")
nclusters.bc=NULL
for (k in 1:10) { 
  if (k==1) {
    nclusters.bc[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist.bc, k)
    nclusters.bc[k]=index.G1(t(data.denoized),data.cluster_temp,  d = data.dist.bc,
                          centrotypes = "medoids")
  }
}
write.table(nclusters.bc,'CH_index-pam-bray-2.tsv',sep = '\t',quote=F)
pdf("CH_index-pam-bray-2.pdf",height=4,width=5)
plot(nclusters.bc, type="b", xlab="k clusters", ylab="Calinski-Harabasz index",main="Optimal number of clusters",lwd = 2)
dev.off()
data.cluster.bc=pam.clustering(data.dist.bc, k=2) 
write.table(data.cluster.16s.bc,"cluster_pam-2.tsv",sep = '\t',quote=F)
#silhouett
obs.silhouette=mean(silhouette(data.cluster.bc, data.dist.bc)[,3])
cat(obs.silhouette) #0.1351978
##pca
obs.pca.bc=dudi.pca(data.frame(t(data.denoized)), scannf=F, nf = 2)
s.class(obs.pca.bc$li, fac=as.factor(data.cluster.bc), grid=F,sub="Principal component analysis", col=c(3,2,4))
##bca
obs.bet.bc=bca(obs.pca.bc, fac=as.factor(data.cluster.bc), scannf=F, nf=2) 
#s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis", col=c(3,2))
#pcoa
obs.pcoa.bc=dudi.pco(data.dist.bc,scannf=F, nf=2)
s.class(obs.pcoa.bc$li, fac=as.factor(data.cluster.bc), grid=F,sub="Principal coordiante analysis", col=c(3,2,4))
