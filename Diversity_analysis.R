################################################################################
#main codes for alpha、belta diversity
################################################################################
rm(list = ls())
##r packages
library("tidyverse")
library("ggplot2")
library("phyloseq")
library("ape")
library("ggthemes")
#Importing phyloseq Data
#OTU # corrected
OTU_table<- read.csv("Long-lived-gut_mycobiome_project/otu.csv",sep=",",header=T, row.names = 1,stringsAsFactors = FALSE)
OTU = otu_table(as.matrix(OTU_table), taxa_are_rows = TRUE)
#Taxa table
tax_table<- read.csv("Long-lived-gut_mycobiome_project/taxonomy.csv",sep=",",header=T, row.names = 1,stringsAsFactors=FALSE)
TAX <- tax_table(as.matrix(tax_table))
#mapping file
metadata <- read.csv("Long-lived-gut_mycobiome_project/metadata.csv", header=TRUE, sep = ",",row.names = "SampleID",stringsAsFactors = F,fileEncoding = "GBK")
sampledata <- sample_data(metadata)
#tree file
random_tree <- read_tree ("Long-lived-gut_mycobiome_project/tree.nwk")
physeq = phyloseq(OTU, TAX)
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)                      
physeq1_rarefy<-rarefy_even_depth(physeq1,rngseed = 711,sample.size = min(sample_sums(physeq1)))
#----------------------------------------------------------------------alpha diversity----------------------------------------------------------------#
richness<-estimate_richness(physeq1_rarefy, measures=c("Observed", "Shannon"))
evenness<-evenness(physeq1_rarefy, index = "pielou", zeroes = TRUE, detection = 0)
#alpha.its<-cbind(richness,evenness)
cbbPalette<-c("#E64B35FF","#4DBBD5FF","#3C5488FF")
richness$Type <- metadata$Type
richness$Type<-factor(richness$Type,levels = c("Young","Old","Long-lived"))
evenness$Type <- evenness$Type
evenness$Type<-factor(evenness$Type,levels = c("Young","Old","Long-lived"))
##
a_my_comparisons <- list(c("Young","Old"),c("Young","Long-lived"),c("Old","Long-lived"))
#plot
figure.Observed<-ggplot(richness,aes(Type, as.numeric(Observed), color = Type))+
  geom_boxplot(aes(color=Type),outlier.shape = NA)+
  geom_jitter(shape=16,position=position_jitter(0.1))+
  ylab("Alpha Diversity Measure")+xlab("")+
  scale_fill_manual(values = cbbPalette)+
  scale_color_manual(values=  cbbPalette)+
  theme_bw()+ggtitle("Observed")+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
       axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text())+
  stat_compare_means(method="wilcox.test",comparisons =a_my_comparisons,label = "p.label")

figure.Shannon<-ggplot(richness,aes(Type, as.numeric(Shannon), color = Type))+
  geom_boxplot(aes(color=Type),outlier.shape = NA)+
  geom_jitter(shape=16,position=position_jitter(0.1))+
  ylab("")+xlab("")+
  scale_fill_manual(values =cbbPalette)+
  scale_color_manual(values=  cbbPalette)+
  theme_bw()+ggtitle("Shannon")+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
       axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text())+
  stat_compare_means(method="wilcox.test",comparisons =a_my_comparisons,label = "p.label")
figure.Pielou<-ggplot(evenness,aes(Type, as.numeric(pielou), color = Type))+
  geom_boxplot(aes(color=Type),outlier.shape = NA)+
  geom_jitter(shape=16,position=position_jitter(0.1))+
  ylab("")+xlab("")+
  scale_fill_manual(values =cbbPalette)+
  scale_color_manual(values=  cbbPalette)+
  theme_bw()+ggtitle("Pielou")+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
       axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text())+
  stat_compare_means(method="wilcox.test",comparisons =a_my_comparisons,label = "p.label")
figure.alpha_diversity<-ggarrange(figure.Observed,figure.Shannon,figure.Pielou,ncol = 3)

#-----------------------------------------------------------------belta diversity----------------------------------------------------------------------------------------------------------------------------#
library(GUniFrac)
library(multcomp)                  
phyloseqin.pcoa <- transform_sample_counts(physeq1_rarefy, function(x) x/sum(x) * 100)
phyloseqin_sub.pcoa<-filter_taxa(phyloseqin.pcoa, function(x) sum(x>0)>=0.2*length(x),prune = TRUE)
#cofunder
data <- data.frame(Sex = metadata$Sex,Age = metadata$Age,BMI=metadata$BMI, Smoking =metadata$smoking,drinking_alcohol=metadata$drinking_alcohol,drinking_tea=metadata$drinking_tea,hypertension=metadata$hypertension)
rownames(data)<-rownames(metadata)
data.nona<-na.omit(data)
rownames(data)<-rownames(metadata)
id<-rownames(data.nona)
phyloseqin_sub.pcoa.df<-subset_samples(phyloseqin_sub.pcoa,sampleID %in% id)
phyloseqin_sub.pcoa.df_tree <- phy_tree(phyloseqin_sub.pcoa.df)
phyloseqin_sub.pcoa.df_otu<-t(otu_table(phyloseqin_sub.pcoa.df))
phyloseqin_sub.pcoa.df_unifracs <- GUniFrac(phyloseqin_sub.pcoa.df_otu, phyloseqin_sub.pcoa.df_tree, alpha=c(0, 0.5, 1))
#unwighted
pcoa.unwighted.unifrac<-phyloseqin_sub.pcoa.df_unifracs$unifracs[,,"d_UW"]
#wighted
#pcoa.wighted.unifrac<-phyloseqin_sub.pcoa.df_unifracs$unifracs[,,"d_1"]
rownames(data.nona) <- rownames(as.matrix(pcoa.unwighted.unifrac))
#apcoa ajust
pcoa.unwighted.unifrac.result <- aPCoA(pcoa.unwighted.unifrac~Sex+BMI+Smoking+drinking_alcohol+drinking_tea+hypertension, data.nona, maincov = Type)
pcoa.unwighted.unifrac.aPCoA.score <- data.frame(pcoa.unwighted.unifrac.result$plotMatrix)                                 
PERMANOVA_metadata <- data.frame(sample_data(phyloseqin_sub.pcoa))
PERMANOVA_dis<-phyloseq::distance(phyloseqin_sub.pcoa, method="unifrac", weighted=FALSE)
adonis_result_dis<-adonis(unname(pcoa.unwighted.unifrac)~ Type, data = data.nona, permutations = 999)
##plot
plot.data<-cbind(pcoa.unwighted.unifrac.aPCoA.score,data.nona)
PC1 = plot.data[,1]
PC2 = plot.data[,2]
plotdata.its <- data.frame(rownames(plot.data),plot.data[,1:2],plot.data$Type)
colnames(plotdata.its) <-c("sample","PC1","PC2","Group")
plotdata.its$Group <- factor(plotdata.its$Group,levels = c("Young","Old","Longevity"),labels = c("Young","Old","Long-lived"))
yf <- plotdata.its
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~Group,data = plotdata.its)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~Group,data = plotdata.its)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
#plot
cbbPalette<-c("#E64B35FF","#4DBBD5FF","#3C5488FF")
plotdata.its$Group=factor(plotdata.its$Group,levels = c("Young","Old","Long-lived"),labels = c("Young","Old","Long-lived"))
p1 <- ggplot(plotdata.its,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  ylim(-0.4,0.5)+
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(),
        axis.line = element_line(), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=20,face = "bold"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        legend.position = "none")

p3 <- ggplot(plotdata.its,aes(Group,PC2,fill=Group)) +
  geom_boxplot() +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 7,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none")

p2<-ggplot(plotdata.its, aes(PC1, PC2,))+
  geom_point(aes(fill=Group),size=6,pch = 21)+
  stat_ellipse(aes(color=Group),level = 0.9)+
  scale_fill_manual(values=cbbPalette)+ 
  scale_color_manual(values = cbbPalette)+
  xlab(paste("PCoA1 ( ",9.43,"%"," )",sep="")) + 
  ylab(paste("PCoA2 ( ",7.20,"%"," )",sep=""))+
  xlim(-0.4,0.5) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=34),
        axis.ticks.length = unit(0.5,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=34,vjust = 7),
        axis.title.y=element_text(colour='black', size=34,vjust = -2),
        axis.text=element_text(colour='black',size=26),legend.position="none") 

 p4 <- ggplot(plotdata.its, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",adonis_result_dis$aov.tab$Df[1],  "\nR2 = ",round(adonis_result_dis$aov.tab$R2[1],4),  "\np-value = ",adonis_result_dis$aov.tab$`Pr(>F)`[1],sep = "")),
            size = 7) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(5,1),ncol = 2,nrow = 2)
                              
                    
                    
               
