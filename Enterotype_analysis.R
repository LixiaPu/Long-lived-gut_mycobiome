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
df<-merge(metadata,genus.its,by = "Sample_ID")
row.names(df)<-df$SampleID
id<-metadata.noerr$sampleID
data.16s.df <- subset(df, SampleID %in% id)
data.16s.df<-data.16s.df%>%dplyr::select(!(1:92))
df.16s.genus<-t(data.16s.df)
data.16s.genus<-apply(df.16s.genus,2,as.numeric)
row.names(data.16s.genus)<-row.names(df.16s.genus)
data.16s.genus<-apply(data.16s.genus,2,function(x) x/sum(x)*100)
data.16s.denoized<-data.16s.genus[apply(data.16s.genus,1,function(x) sum(x>0)>=0.2*length(x)),]
#data.16s.denoized=noise.removal(data.16s.genus,percent=0.0001)
data.16s.denoized<-data.16s.denoized[!grepl("g__metagenome.1",row.names(data.16s.denoized)),]
#PAM cluster +bary
data.dist.bc <- vegdist(t(data.16s.denoized),method = "bray")
write.table(as.matrix(data.dist.bc),"bray_curtis-genus_16s.txt",sep = '\t',quote = FALSE,col.names = NA)
#Calinski-Harabasz 
nclusters.bc = index.G1(t(data.16s.denoized), data.cluster, d = data.dist.bc, centrotypes = "medoids")
nclusters.bc=NULL
for (k in 1:10) { 
  if (k==1) {
    nclusters.bc[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist.bc, k)
    nclusters.bc[k]=index.G1(t(data.16s.denoized),data.cluster_temp,  d = data.dist.bc,
                          centrotypes = "medoids")
  }
}
write.table(nclusters.bc,'CH_index_16s-pam-bray-3.tsv',sep = '\t',quote=F)
pdf("CH_index_16s-pam-bray-3.pdf",height=4,width=5)
plot(nclusters.bc, type="b", xlab="k clusters", ylab="Calinski-Harabasz index",main="Optimal number of clusters",lwd = 2)
dev.off()
data.cluster.16s.bc=pam.clustering(data.dist.bc, k=2) 
write.table(data.cluster.16s.bc,"cluster_16s_pam-3.tsv",sep = '\t',quote=F)
#silhouett
obs.silhouette=mean(silhouette(data.cluster.16s.bc, data.dist.bc)[,3])
cat(obs.silhouette) #0.1351978
##pca
obs.pca.16s.bc=dudi.pca(data.frame(t(data.16s.denoized)), scannf=F, nf = 2)
s.class(obs.pca.16s.bc$li, fac=as.factor(data.cluster.16s.bc), grid=F,sub="Principal component analysis", col=c(3,2,4))
##bca
obs.bet.16s.bc=bca(obs.pca.16s.bc, fac=as.factor(data.cluster.16s.bc), scannf=F, nf=2) 
#s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis", col=c(3,2))
#pcoa
obs.pcoa.16s.bc=dudi.pco(data.dist.bc,scannf=F, nf=2)
s.class(obs.pcoa.16s.bc$li, fac=as.factor(data.cluster.16s.bc), grid=F,sub="Principal coordiante analysis", col=c(3,2,4))
