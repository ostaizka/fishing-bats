install.packages("caret")
install.packages("rpart")
install.packages("xgboost")
install.packages("randomForest")
install.packages("kernlab")
install.packages("LiblineaR")
install.packages("pROC")
install.packages("dyverse")
install.packages("cowplot")
install.packages("reshape2")


library(caret)
library(rpart)
library(xgboost)
library(randomForest)
library(kernlab)
library(LiblineaR)
library(pROC)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(vegan)
library(gtools)
library(reshape2)
library(hilldiv)
library(caretEnsemble)

setwd("/Users/anttonalberdi/Downloads/")

counts <- read.table("ASVs_counts_bacteriaclass.tsv",sep="\t",header=TRUE,row.names=1)
hierarchy <- read.csv("samples.csv",sep=";")
hierarchy[,1] <- as.character(hierarchy[,1])
hierarchy[,3] <- as.character(hierarchy[,3])

#Remove singletons
species <- as.character(unique(hierarchy$Species))
for (sp in species){
samples <- as.character(hierarchy[hierarchy$Species == sp,"Samples"])
counts[rowSums(counts[,samples] != 0) == 1,samples] <- 0
}
counts <- counts[apply(counts, 1, function(x) !all(x==0)),]

#Transform to 0-1
counts_tss <- tss(counts)

#Transpose
counts_tss_t <- t(counts_tss)

#Create working data table (group classification (diet) in the end)
data <- merge(counts_tss_t,hierarchy[,c(1,3)],by.x="row.names",by.y="Samples")
rownames(data) <- data[,1]
data <- data[,-1]
# We want the diagnosis column to be a factor
data$diet <- factor(data$diet)


uncertain <- hierarchy[hierarchy$Species == "Myotis_capaccinii",1]
certain <- hierarchy[hierarchy$Species != "Myotis_capaccinii",1]
trainTransformed <- data[certain,]
testTransformed  <- data[uncertain,]

folds <- 7
cvIndex <- createMultiFolds(factor(trainTransformed$diet), folds, times=100)
cv <- trainControl(method="repeatedcv",
                     number=folds,
                     index = cvIndex,
                     returnResamp="final",
                     classProbs=TRUE,
                     summaryFunction=twoClassSummary,
                     indexFinal=NULL,
                     savePredictions = TRUE)

#TuneGrids
grid_rf <-  expand.grid(mtry = c(80,500,1000,1500))
grid_lr <- expand.grid(cost = c(0.0001, 0.001, 0.0025, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 10),
                     loss = "L2_primal",
                     epsilon = 0.01)
grid_xgb <-  expand.grid(nrounds=500,
              gamma=0,
              eta=c(0.001, 0.01, 0.1, 1),
              max_depth=8,
              colsample_bytree= 0.8,
              min_child_weight=1,
              subsample=c(0.4, 0.5, 0.6, 0.7))

tuneList <- list(
    rf=caretModelSpec(method="rf", tuneGrid=grid_rf),
    lr=caretModelSpec(method="regLogistic", tuneGrid=grid_lr),
    xgb=caretModelSpec(method="xgbTree", tuneGrid=grid_xgb))

model_list <- caretList(
    diet~., data=trainTransformed,
    trControl=cv,
    methodList=c("rf", "regLogistic","xgbTree"),
    tuneList=tuneList)

results <- resamples(model_list[c(1:3)])
summary(results)

# Ensemble the predictions of `models` to form a new combined prediction based on glm
glm_ensemble <- caretStack(
  model_list[c(1:3)],
  method="glm",
  metric="ROC",
  trControl=trainControl(
    method="boot",
    number=100,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=twoClassSummary
  )
)

# Variable importance
varimp_ensemble <- varImp(caretEnsemble(model_list[c(1:3)]))
varimp_ensemble <- varimp_ensemble[order(varimp_ensemble$overall,decreasing=TRUE),]
head(varimp_ensemble,20)

# Predict test data
predicted_testdata_ensemble <- predict(glm_ensemble,testTransformed)

#Individual predictions
sample_prediction_ensemble <- cbind(Sample=rownames(testTransformed),Prediction=as.character(predicted_testdata_ensemble))
sample_prediction_ensemble <- merge(sample_prediction_ensemble,hierarchy[,c("Samples","Species")],by.x="Sample",by.y="Samples")

# Compute the confusion matrix
confmat_ensemble_train <- confusionMatrix(reference = trainTransformed$diet, data = predicted_traindata_ensemble, mode='everything', positive='Insectivorous')
confmat_ensemble <- confusionMatrix(reference = testTransformed$diet, data = predicted_testdata_ensemble, mode='everything', positive='Insectivorous')

##########
# Sankey plots
##########

library(data.table)
library(cluster)
library(factoextra)
library(plotly)

setwd("/Users/anttonalberdi/Downloads/")
pairdis.q1 <- readRDS("pairdis.q1.RData")
pairdis.q1.UqN <- pairdis.q1$L1_UqN
hierarchy <- read.csv("samples.csv",sep=";")
hierarchy[,1] <- as.character(hierarchy[,1])
hierarchy[,3] <- as.character(hierarchy[,3])

#####
# Piscivorous
#####

#Subset distance table
pairdis.q1.UqN_fish <- as.dist(pairdis.q1.UqN[hierarchy[hierarchy$Group == "Fishing",1],hierarchy[hierarchy$Group == "Fishing",1]])

#Determine optimal number of clusters
clus_fish <- c()
for (k in c(2:10)){
pam_fish <- pam(pairdis.q1.UqN_fish, k, stand = FALSE)
pam_fish_sil <- silhouette(pam_fish, pairdis.q1.UqN_fish)
sil <- mean(aggregate(pam_fish_sil[,3],by=list(pam_fish_sil[,1]),FUN=mean)[,2])
clus_fish <- rbind(clus_fish,c(k=k,sil=sil))
}

# Create cluster table
k=4
clusters_fish <- pam(pairdis.q1.UqN_fish, k, stand = FALSE)$clustering
clusters_fish <- as.data.frame(cbind(as.character(names(clusters_fish)),as.factor(clusters_fish)))
clusters_fish2 <- merge(clusters_fish,hierarchy[,c(1,5)],by.x="V1",by.y="Samples")
clusters_fish2[,3] <- as.character(clusters_fish2[,3])
clusters_fish_table <- table(clusters_fish2[,c(2:3)])

# Create Sanky plot
clusnumber <- 4
clusters <- paste("Cluster",c(1:clusnumber),sep="")

source = sort(rep(c(0:(ncol(clusters_fish_table)-1)),clusnumber))
target = rep(c(ncol(clusters_fish_table):(ncol(clusters_fish_table)+(clusnumber-1))),ncol(clusters_fish_table))
value =  melt(clusters_fish_table)$value

sankey_fish <- plot_ly(
    type = "sankey",
    orientation = "h",

    node = list(
      label = c("Myotis_capaccinii", "Myotis_pilosus", "Myotis_vivesi", "Noctilio_leporinus","Cluster1","Cluster2","Cluster3","Cluster4"),
      color = c("#1f3d8a", "#0a9a9a", "#4dc2f1", "#5582c2", "red", "red", "red", "red"),
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),

    link = list(
      source = source,
      target = target,
      value =  value
    )
  )

sankey_fish

#####
# Arthropodivorous
#####

#Subset distance table
pairdis.q1.UqN_insect <- pairdis.q1.UqN[hierarchy[(hierarchy$Group == "Insectivorous") | (hierarchy$Group == "Uncertain"),1], hierarchy[(hierarchy$Group == "Insectivorous") | (hierarchy$Group == "Uncertain"),1]]
pairdis.q1.UqN_insect <- as.dist(pairdis.q1.UqN_insect[-which(rownames(pairdis.q1.UqN_insect) %in% c("PO254_S24","GM1.H17")),-which(colnames(pairdis.q1.UqN_insect) %in% c("PO254_S24","GM1.H17"))])

#Determine optimal number of clusters
clus_insect <- c()
for (k in c(2:20)){
pam_insect <- pam(pairdis.q1.UqN_insect, k, stand = FALSE)
pam_insect_sil <-silhouette(pam_insect, pairdis.q1.UqN_insect)
sil <- mean(aggregate(pam_insect_sil[,3],by=list(pam_insect_sil[,1]),FUN=mean)[,2])
clus_insect <- rbind(clus_insect,c(k=k,sil=sil))
}
clus_insect
