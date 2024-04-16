library("randomForest")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("caret")
library("pROC")

rm(list = ls())

OTU_table<- read.csv("Long-lived-gut_mycobiome_project/otu.csv",sep=",",header=T, row.names = 1,stringsAsFactors = FALSE)
metadata <- read.csv("Long-lived-gut_mycobiome_project/metadata.csv", header=TRUE, sep = ",",row.names = "SampleID",stringsAsFactors = F,fileEncoding = "GBK")
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = FALSE, scale = TRUE)
richness<-estimate_richness(physeq1_rarefy, measures=c("Observed"))
otu_table_scaled<-cbind(t(otu_table_scaled),richness)
otu_table_scaled_all <- data.frame((otu_table_scaled))
otu_table_scaled_all$Type <- as.factor(metadata[rownames(otu_table_scaled_all), "Type"])
otu_table_scaled_all_used <- data.frame(otu_table_scaled_all)
set.seed(151)
train_use <- sample(nrow(otu_table_scaled_all_used), nrow(otu_table_scaled_all_used)*0.7)
otu_table_scaled_train <- otu_table_scaled_all_used[train_use,]
set.seed(151)
##The first step is defining the parameters we want to use for training: Grid search for mtry
##choose the best mtry number trControl 
train(x=otu_table_scaled_train[,1:(ncol(otu_table_scaled_train)-1)], y=factor(otu_table_scaled_train$Type),
      method = 'rf',# Use the 'random forest' algorithm
      trControl = trainControl(method='cv', 
                               number=10,
                               search='grid'))
###choose the best tree number
set.seed(151)
RF_state_classify <- randomForest(x=otu_table_scaled_train[,1:(ncol(otu_table_scaled_train)-1)] , y=factor(otu_table_scaled_train$Type),ntree=5000, importance=TRUE, proximity=TRUE,replace=TRUE,mtry=47,na.action=na.omit)
plot(RF_state_classify)
#Identify important features
set.seed(151)
RF_state_classify_s <- randomForest(x=otu_table_scaled_train[,1:(ncol(otu_table_scaled_train)-1)] , y=factor(otu_table_scaled_train$Type),ntree=1000, importance=TRUE, proximities=TRUE,replace=TRUE,mtry=47,na.action=na.omit,oob_score=TRUE)
imp <- importance(RF_state_classify_s)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseAccuracy))
importance_otu.select <- imp.sort[1:28,]##28 MeanDecreaseAccuracy>5
##prepare dataset based on the best model
otu_id.select <- importance_otu.select$predictors
otu.select <- subset(otu_table_scaled_all_used,select=c(as.character(otu_id.select),"Type"))
##l vs o select old and long-lived individuals
otu.select.lo<-otu.select[1:156,]
##training dataset(70%) and validation dataset (30%)
set.seed(151)
train <- sample(nrow(otu.select.lo), nrow(otu.select.lo)*0.7)
otu_train.select <- otu.select.lo[train, ]
otu_train.select.data<-apply(otu_train.select[,-ncol(otu_train.select)], c(1,2), function(x) as.numeric(as.character(x)))
otu_test.select <- otu.select.lo[-train, ]
otu_test.select.data<-apply(otu_test.select[,-ncol(otu_test.select)], c(1,2), function(x) as.numeric(as.character(x)))
###
##randomForest model
##find the best mtry
set.seed(151)
##choose the best mtry number ##15
train(x=otu_train.select.data , y=factor(otu_train.select$Type),
      method = 'rf',# Use the 'random forest' algorithm
      trControl = trainControl(method='cv', 
                               number=10, 
                               search='grid'))
                            
set.seed(151)
otu_train.select.forest <- randomForest(x = otu_train.select.data,y=factor(otu_train.select$Type), importance = TRUE, ntree=1000,mtry=15, proximity=TRUE)
##
train_predict <- predict(otu_train.select.forest, otu_train.select,type="response")
test_predict <- predict(otu_train.select.forest, otu_test.select,type="response")
##
confusionMatrix(train_predict, factor(otu_train.select$Type))
confusionMatrix(test_predict, factor(otu_test.select$Type))
##
library("pROC")
roc.info_test.OvsL<-roc(factor(otu_test.select$Type),
    test_predict[,1],
    plot=TRUE,
    ci.method="bootstrap",
    ci =TRUE,
    conf.level =0.95
    )
#
roc.info_train.OvsL<-roc(factor(otu_train.select$Type), 
    train_predict[,1], 
    plot=FALSE, 
    ci=TRUE,
    ci.method="bootstrap",
    conf.level =0.95
)


