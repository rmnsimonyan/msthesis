# Training on GDSC validating on CTRP (same cell lines as in CCLE but proper RNA seq and suseptibility data)
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

# Run DE analysis to create signitures

# figure out matrix
cellline = colnames(sysactualnumericmat_train)

# fix the top genes from DE to be considered
topgenes = 10

for(i in colnames(factorsensitivitytablesangermat_train)){
  tryCatch({
  print(i)
  
  status_ = factorsensitivitytablesangermat_train[!is.na(factorsensitivitytablesangermat_train[,i]),i]
  
  mm <- model.matrix(~status_) #cellline + 
  
  fit <- lmFit(sysactualnumericmat_train[,!is.na(factorsensitivitytablesangermat_train[,i])], mm)
  
  head(coef(fit))
  fit <- eBayes(fit)
  
  top.table <- topTable(fit, coef = "status_S", n = Inf)
  
  #rownames(top.table) = genes[rownames(top.table),"gene_symbol"]
  write_delim(top.table, paste("Trained_and_validated_on_GDSC/", i, ".txt", sep = ""))
  #head(top.table, 20)
  down = as.data.frame(t(c(paste(i, "_down", sep = ""), "", top.table$ID[top.table$adj.P.Val<0.01 & top.table$logFC<0][1:topgenes])))
  names(down) = NULL
  write_delim(down, paste("Trained_and_validated_on_GDSC/",i, "_down.txt", sep = ""), col_names=F)
  
  up = as.data.frame(t(c(paste(i, "_up", sep = ""), "", top.table$ID[top.table$adj.P.Val<0.01 & top.table$logFC>0][1:topgenes])))
  names(up) = NULL
  write_delim(up, paste("Trained_and_validated_on_GDSC/",i, "_up.txt", sep = ""), col_names=F)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }



# filtering out 10 (or n) top genes here for each signiture
#topgenes = 10

#for(drugname in sangerdrugsforgmt){
#  tryCatch({
#    print(drugname)
#    
#    DEgenetable = read.table(paste("DE_outputs_sanger_whole_regression/",drugname, ".txt", sep = ""), header=TRUE, sep=" ", row.names = 1)
#    DEgenetable = DEgenetable[DEgenetable$P.Value<0.01,]
    #    DEgenetable = DEgenetable[DEgenetable$logFC>0,]
#    
#    up = as.data.frame(t(c(paste(drugname, "_up", sep = ""), "", rownames(DEgenetable[DEgenetable$logFC>0,])[1:topgenes][!is.na(rownames(DEgenetable[DEgenetable$logFC>0,])[1:topgenes])])))
#    names(up) = NULL
#    write_delim(up, paste("DE_outputs_sanger_whole_regression_top", topgenes, "/",drugname, "_up.txt", sep = ""), col_names=F)
#    
#    
#    down = as.data.frame(t(c(paste(drugname, "_down", sep = ""), "", rownames(DEgenetable[DEgenetable$logFC<0,])[1:topgenes][!is.na(rownames(DEgenetable[DEgenetable$logFC<0,])[1:topgenes])])))
#    names(down) = NULL
#    write_delim(down, paste("DE_outputs_sanger_whole_regression_top", topgenes, "/",drugname, "_down.txt", sep = ""), col_names=F)
#  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#  
#}


# Then run the following in termnial to produce the combined gmt file for further analysis in fGSEA
# cat $(ls | grep "up\|down")> combined.gmt.txt


# mean to compare to
sysDEmean = rowMeans(sysactualnumericmat_test)
sysFoldchagen = sysactualnumericmat_test/sysDEmean
syslogFoldchange = log(sysFoldchagen)
cellinenamesforfoldchangetable = colnames(sysactualnumericmat_test)


fGSEAhuman = function(results, filter="",samplename)  {
  pathways = pathways <- gmtPathways("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Trained_and_validated_on_GDSC/combined.gmt.txt")  #("msigdb.v7.0.symbols.gmt.txt")
  pathways = pathways[grepl(filter, names(pathways))]
  # str(head(pathways))
  fgseaRes <- suppressWarnings(fgsea(pathways, stats=results, nperm=100000, minSize = 8, maxSize = 2000, BPPARAM=MulticoreParam(parallel::detectCores())))
  #  fgseaRes = fgseaRes[padj < 0.05, ]
  topPathwaysUp <- fgseaRes[ES > 0, ][head(order(pval), n=15), pathway]
  topPathwaysDown <- fgseaRes[ES < 0, ][head(order(pval), n=15), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  if (length(topPathways) > 0) {
    #' `#+ fig.width=10, fig.height=10`
    pdf(file = paste("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Trained_and_validated_on_GDSC", "/fgsea/", samplename, ".pdf", sep = ""), width = 8, height = 11) # defaults to 7 x 7 inches
    plotGseaTable(pathways[topPathways], results, fgseaRes, gseaParam = 0.5, colwidths = c(7, 3, 0.8, 1.2, 1.2))
    knitr::kable(format(paste0("[", fgseaRes$pathway,"](", "http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",fgseaRes$pathway, ") ", fgseaRes$pval, " ", ifelse(fgseaRes$ES < 0, "downregulated", "upregulated")), caption = main))
    dev.off()
    write_delim(fgseaRes, paste("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Trained_and_validated_on_GDSC", "/fgsea/", samplename, ".txt", sep = ""))
  } else {print(paste("out of",length(pathways),"no significant", filter, "pathways"))}
}


for(i in 1:length(cellinenamesforfoldchangetable)){
  print(cellinenamesforfoldchangetable[i])
  fGSEAhuman(syslogFoldchange[,i], samplename = cellinenamesforfoldchangetable[i])
}


predictionmatrix = matrix(,nrow=length(cellinenamesforfoldchangetable),ncol=length(colnames(factorsensitivitytablesangermat_test)))
rownames(predictionmatrix) = cellinenamesforfoldchangetable
colnames(predictionmatrix) = colnames(factorsensitivitytablesangermat_test)

for(enrichedcellline in cellinenamesforfoldchangetable){
  tryCatch({
    print(enrichedcellline)
    enrichment = read.table(paste("Trained_and_validated_on_GDSC", "/fgsea/",enrichedcellline, ".txt", sep = ""), header=TRUE, sep=" ", row.names = 1)
    
    # the following two lines take only the signiture for each drug (up or down) that was enriched with a lower p value
    for(j in unique(unlist(strsplit(rownames(enrichment), "_"))[c(TRUE, FALSE)])){
      
      if(sum(unlist(strsplit(rownames(enrichment), "_"))[c(TRUE, FALSE)] %in% j)==2){
        coordis = which(unlist(strsplit(rownames(enrichment), "_"))[c(TRUE, FALSE)] %in% j)
        if(enrichment[coordis[1],]$padj>enrichment[coordis[2],]$padj){
          enrichment[coordis[1],]$padj=111
        }else{
          enrichment[coordis[2],]$padj=111
        }
      }
    }
    
    enrichment = enrichment[enrichment$padj!=111,]
    #enrichment = enrichment[enrichment$padj<0.05,]
    # cut from here if unneccesary
    
    enpathway = unlist(strsplit(rownames(enrichment), "_"))[c(TRUE, FALSE)]
    updown = as.numeric(as.factor(unlist(strsplit(rownames(enrichment), "_"))[c(FALSE, TRUE)]))*2-3
    #    updownregulated = as.numeric(enrichment$NES>0)
    #    suspttable = ifelse(updown != updownregulated, "S", "R")
    suspttable = (-updown*enrichment$ES+1)/2
    names(suspttable) = enpathway
    predictionmatrix[enrichedcellline,enpathway] = suspttable
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

#predictionmatrix[is.na(predictionmatrix)] = 0

# overlaptables
truematrix = factorsensitivitytablesangermat_test[which(rownames(factorsensitivitytablesangermat_test) %in% rownames(predictionmatrix)),which(colnames(factorsensitivitytablesangermat_test) %in% colnames(predictionmatrix))]
predictmatrix = predictionmatrix[which(rownames(predictionmatrix) %in% rownames(factorsensitivitytablesangermat_test)),which(colnames(predictionmatrix) %in% colnames(factorsensitivitytablesangermat_test))]
predictmatrix = predictmatrix[rownames(truematrix),colnames(truematrix)]
levels(truematrix) <- c(0,1)

# it gives worse results with binarizing the prediction table
#binpredictmatrix = predictmatrix
#binpredictmatrix[which(binpredictmatrix>0.5)]=1
#binpredictmatrix[which(binpredictmatrix<0.5)]=0

roc(truematrix, predictmatrix, plot=T)


for(roci in 1:dim(truematrix)[2]){
  tryCatch({
    pdf(file = paste("Trained_and_validated_on_GDSC", "/rocs/", colnames(truematrix)[roci], ".pdf", sep = ""), width = 8, height = 11) # defaults to 7 x 7 inches
    roccurve = roc(truematrix[,roci], predictmatrix[,roci], plot=T)
    title(main = paste(colnames(truematrix)[roci],roccurve$auc, sep=" "))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


for(roci in c("Dabrafenib","Palbociclib","Trametinib", "Talazoparib")){
  pdf(file = paste("Trained_and_validated_on_GDSC/rocs/", roci, ".pdf", sep = ""), width = 8, height = 11) # defaults to 7 x 7 inches
  roccurve = roc(truematrix[,roci], predictmatrix[,roci], plot=T)
  title(main = paste(roci,roccurve$auc, sep=" "))
  dev.off()
}

roci ="Dabrafenib"
roc(truematrix[,roci], predictmatrix[,roci], plot=T)






