library("randomForest")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("caret")
library("pROC")

rm(list = ls())

feature <- read.table('Long-lived-gut_mycobiome_project/level-6-20-44_66-85_100-117.txt',header = T,sep = "\t",check.names = T)
set.seed(151)
train_use <- sample(nrow(otu_table_scaled_all_used), nrow(otu_table_scaled_all_used)*0.7)
otu_table_scaled_train <- otu_table_scaled_all_used[train_use,]
