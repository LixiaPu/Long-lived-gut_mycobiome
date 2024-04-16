################################################################################
#main codes for alpha、belta diversity
################################################################################
rm(list = ls())
##r packages
library("tidyverse")
library("ggplot2")
library("phyloseq")
library("ape")

#Importing phyloseq Data
#OTU # corrected
OTU_table<- read.csv("Long-lived-gut_mycobiome_project/otu.xls",sep="\t",header=T, row.names = 1,strip.white = T, stringsAsFactors = F)
otu<-data.frame(apply(OTU_table, 2, function(x) as.numeric(as.character(x))))
rownames(otu) <- row.names(OTU_table)
OTU = otu_table(otu, taxa_are_rows = TRUE)
#Taxa table
tax_table<- read.csv("Long-lived-gut_mycobiome_project/taxonomy.csv",sep=",",header=T, row.names = 1,stringsAsFactors=FALSE)
TAX <- tax_table(as.matrix(tax_table))
#mapping file
metadata <- read.csv("Long-lived-gut_mycobiome_project/metadata.csv", header=TRUE, sep = ",",row.names = "SampleID",stringsAsFactors = F,fileEncoding = "GBK")
sampledata <- sample_data(metadata)
#tree file
random_tree <- read_tree ("Long-lived-gut_mycobiome_project/tree.nwk")
rooted_tree<-random_tree
physeq = phyloseq(OTU, TAX)
physeq1 = merge_phyloseq(physeq, sampledata, rooted_tree)                      
physeq1_rarefy<-rarefy_even_depth(phyloseqin_sub,rngseed = 711,sample.size = min(sample_sums(physeq1)))
###
richness<-estimate_richness(physeq1_rarefy, measures=c("Observed", "Shannon"))
evenness<-evenness(physeq1_rarefy, index = "pielou", zeroes = TRUE, detection = 0)
alpha.its<-cbind(richness,evenness)
mypal = pal_npg("nrc",alpha = 1)(3)
alpha.its$Type <- metadata$Type
alpha.its$Type<-factor(alpha.its$Type,levels = c("Young","Old","Long-lived"))
#plot
df1 <- data.frame(a=c(1,1,2,2), b=c(1255,1270,1270,1255))
df2 <- data.frame(a=c(2,2,3,3), b=c(1455,1470,1470,1455))
df3 <- data.frame(a=c(1,1,3,3), b=c(1655,1670,1670,1655))
figure2c.Chao1<-ggplot(alpha.its,aes(Type, as.numeric(Chao1), color = Type))+
   geom_boxplot(aes(color=Type),outlier.shape = NA)+
  geom_jitter(shape=16,position=position_jitter(0.1))+
  ylab("Alpha Diversity Measure")+xlab("")+
  scale_fill_manual(values =mypal)+
  scale_color_manual(values= mypal)+
  theme_bw()+ggtitle("Chao1")+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
       axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text())+
  geom_line(data=df1,aes(a,b),cex=0.5,color="black")+
  geom_line(data=df2,aes(a,b),cex=0.5,color="black")+
  geom_line(data=df3,aes(a,b),cex=0.5,color="black")+
  annotate('text',x=1.5,y=1300,label='p-value:0.137',size=4,color="black")+
  annotate('text',x=2.5,y=1500,label='p-value:1.44e-06',size=4,color="black")+
  annotate('text',x=2.0,y=1700,label='p-value:3.52e-08',size=4,color="black")

##belta diversity
library(GUniFrac)
library(multcomp)
                      
phyloseqin.pcoa <- transform_sample_counts(physeq1_rarefy, function(x) x/sum(x) * 100)
phyloseqin_sub.pcoa<-filter_taxa(phyloseqin.pcoa, function(x) sum(x>0)>=0.2*length(x),prune = TRUE)
#协变量
data <- data.frame(Sex = metadata$Sex,Age = metadata$Age,BMI=metadata$BMI, Smoking =metadata$smoking,drinking_alcohol=metadata$drinking_alcohol,drinking_tea=metadata$drinking_tea,hypertension=metadata$hypertension)
rownames(data)<-rownames(metadata)
data.nona<-na.omit(data)
rownames(data)<-rownames(metadata)
id<-rownames(data.nona)
phyloseqin_sub.pcoa.df<-subset_samples(phyloseqin_sub.pcoa,sampleID %in% id)  
  
                      
                    
                    
               
