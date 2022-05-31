# xgboost ML for the data prediction on GDSC
library("cmapR")
library(ggplot2)
library(limma)
library(ggfortify)
library(ActivePathways)
library(ArrayExpress)
library("fgsea")
library(BiocParallel)
library(edgeR)
library(readr)
library("readxl")
library(stringr)
library(pROC)
library(tidyverse)
library(caret)
library(rpart)
library(xgboost)
library(DiagrammeR)
library(caret)


# Reading the data of Sanger dataset to predict drug candidates for the patients and do systematic validation
sysverif = read_delim("Validation dataset from Sanger/Cell_line_RMA_proc_basalExp.txt")
sysverifmat = as.matrix(sysverif, rownames=T)
rownames(sysverifmat) = sysverifmat[,1]
sysverifmat = sysverifmat[,-c(1,2)]
sysnumericmat = as.numeric(sysverifmat)
sysactualnumericmat <- matrix(data=sysnumericmat, ncol=ncol(sysverifmat), nrow=nrow(sysverifmat))
dimnames(sysactualnumericmat) <- list(rownames(sysverifmat), colnames(sysverifmat))



#Read the cell line annotations to match the response ids
sangercellineannotations = read_csv("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Validation\ dataset\ from\ Sanger/Cell_listTue\ Mar\ 22\ 20_34_26\ 2022.csv")
cellinenamesfortable = sangercellineannotations$Name[match(colnames(sysactualnumericmat), sangercellineannotations$COSMIC_ID...1)]

# substitute cell-line ids with cell-line names in the expression matrix
colnames(sysactualnumericmat) = cellinenamesfortable

# remove duplicate NAs
sysactualnumericmat = sysactualnumericmat[,which(!is.na(colnames(sysactualnumericmat)))]

#read the sensitivity matrix
sensitivitytablesanger = read_xlsx("Validation dataset from Sanger/TableS5C.xlsx", sheet=2)
sensitivitytablesangermat = as.matrix(sensitivitytablesanger, rownames=T)
rownames(sensitivitytablesangermat) = sensitivitytablesanger$`Screened Compounds:`
sensitivitytablesangermat = sensitivitytablesangermat[,-1]
factorsensitivitytablesangermat = as.factor(sensitivitytablesangermat) 
dim(factorsensitivitytablesangermat) <- c(nrow(sensitivitytablesangermat), ncol(sensitivitytablesangermat))
dimnames(factorsensitivitytablesangermat) <- list(rownames(sensitivitytablesangermat), colnames(sensitivitytablesangermat))
factorsensitivitytablesangermat = factorsensitivitytablesangermat[, !duplicated(colnames(factorsensitivitytablesangermat))]


drugnamealias = read_delim("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Validation\ dataset\ from\ Sanger/listofdrugs.txt")

# fixing the naming of drugs in the sensitivity matrix
for(drn in 1:ncol(factorsensitivitytablesangermat)){
  #print(colnames(factorsensitivitytablesangermat)[drn])
  
  print(strsplit(colnames(factorsensitivitytablesangermat)[drn], "[..]")[[1]][1])
  for(i in 1:length(drugnamealias$drug_id)){
    if(grepl(strsplit(colnames(factorsensitivitytablesangermat)[drn], "[..]")[[1]][1], drugnamealias$drug_name[i], fixed = TRUE) || grepl(strsplit(colnames(factorsensitivitytablesangermat)[drn], "[..]")[[1]][1], drugnamealias$synonyms[i], fixed = TRUE)){
      colnames(factorsensitivitytablesangermat)[drn] = drugnamealias$drug_name[i]
      break
    }
  }
  
  
}


# match the expression table cell lines with the response table
sysactualnumericmat = sysactualnumericmat[,colnames(sysactualnumericmat)[which(colnames(sysactualnumericmat)%in%rownames(factorsensitivitytablesangermat))]]
sysactualnumericmat_train = sysactualnumericmat[,1:622]
sysactualnumericmat_test = sysactualnumericmat[,623:922]

factorsensitivitytablesangermat = factorsensitivitytablesangermat[colnames(sysactualnumericmat)[which(colnames(sysactualnumericmat)%in%rownames(factorsensitivitytablesangermat))],]
factorsensitivitytablesangermat_train = factorsensitivitytablesangermat[1:622,]
factorsensitivitytablesangermat_test = factorsensitivitytablesangermat[623:922,]

write.table(sysactualnumericmat_train, file="sysactualnumericmat_train.txt", quote = F)
write.table(sysactualnumericmat_test, file="sysactualnumericmat_test.txt", quote = F)
write.table(factorsensitivitytablesangermat_train-1, file="factorsensitivitytablesangermat_train.txt", quote = F)
write.table(factorsensitivitytablesangermat_test, file="factorsensitivitytablesangermat_test.txt", quote = F)




for(i in colnames(factorsensitivitytablesangermat_train)[125:265]){
  tryCatch({
  labls_train = factorsensitivitytablesangermat_train[!is.na(factorsensitivitytablesangermat_train[,i]),i]
  labls_train = as.numeric(labls_train)-1
  matrix_train = t(sysactualnumericmat_train[,!is.na(factorsensitivitytablesangermat_train[,i])])
  
  labls_test = factorsensitivitytablesangermat_test[!is.na(factorsensitivitytablesangermat_test[,i]),i]
  labls_test = as.numeric(labls_test)-1
  matrix_test = t(sysactualnumericmat_test[,!is.na(factorsensitivitytablesangermat_test[,i])])
  
  write.table(labls_train, file=paste("xgboostinput/labls_train_", i, ".txt", sep=""), quote = F)
  write.table(matrix_train, file=paste("xgboostinput/matrix_train_", i, ".txt", sep=""), quote = F)
  write.table(labls_test, file=paste("xgboostinput/labls_test_", i, ".txt", sep=""), quote = F)
  write.table(matrix_test, file=paste("xgboostinput/matrix_test_", i, ".txt", sep=""), quote = F)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


write.table(colnames(factorsensitivitytablesangermat_train), file="xgboostinput/drugnames.txt", quote=F)


train_control = trainControl(method = "cv", number = 5)
# Customsing the tuning grid
gbmGrid <-  expand.grid(max_depth = c(3, 5, 7), 
                        nrounds = (1:10)*50,    # number of trees
                        # default values below
                        eta = 0.3,
                        gamma = 0,
                        subsample = 1,
                        min_child_weight = 1,
                        colsample_bytree = 0.6)

# training a XGboost Regression tree model while tuning parameters
model = train(Cost~., data = train, method = "xgbTree", trControl = train_control, tuneGrid = gbmGrid)







dtrain <- xgb.DMatrix(data = round(matrix_train), label=labls_train)
dtest <- xgb.DMatrix(data = round(matrix_test), label=labls_test)

watchlist <- list(train=dtrain, test=dtest)

bst <- xgb.train(data=dtrain, max.depth=15, eta=1, nthread = 2, nrounds=15, watchlist=watchlist, objective = "binary:logistic")

bst <- xgb.train(data=dtrain, max.depth=2, eta=1, nthread = 2, nrounds=2, watchlist=watchlist, eval.metric = "error", eval.metric = "logloss", objective = "binary:logistic")



model = xgb.train(data = dtrain, max.depth = 3, watchlist=watchlist, nrounds = 100)


model_xgboost = xgboost(data = dtrain, max.depth = 3, nrounds = 16, verbose = 0)
summary(model_xgboost)

importance_matrix = xgb.importance(colnames(dtrain), model = model)


dtrain <- xgb.DMatrix(data = matrix, label=labls)

bstSparse <- xgboost(data = matrix, label = labls, max.depth = 6, eta = 1, nthread = 2, nrounds = 3, objective = "binary:logistic", verbose = 2)
bst <- xgb.train(data = xgb.DMatrix(matrix), label = labls, max.depth = 6, eta = 1, nthread = 2, nrounds = 3, objective = "binary:logistic", verbose = 2)

bst <- xgb.train(data=dtrain, booster = "gbtree", nthread = 2, nrounds=2, objective = "binary:logistic")

importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)


data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')

bst1 <- xgboost(data = agaricus.train$data, label = agaricus.train$label, max_depth = 2,
               eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")
xgb.importance(model = bst1)


bst1 <- xgboost(data = matrix, label = labls, max_depth = 2,
                eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")

importance_matrix = xgb.importance(colnames(matrix), model = bst)




labels = c(0,
            1,
            1,
            0,
            0,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            1,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            1,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0)

predicted = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


roc(labels, predicted, plot=T)

