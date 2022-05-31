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


# Reading the data of Sanger dataset to predict drug candidates for the patients and do systematic validation
sysverif = read_delim("Validation dataset from Sanger/Cell_line_RMA_proc_basalExp.txt")
sysverifmat = as.matrix(sysverif, rownames=T)
rownames(sysverifmat) = sysverifmat[,1]
sysverifmat = sysverifmat[,-c(1,2)]
sysnumericmat = as.numeric(sysverifmat)
sysactualnumericmat <- matrix(data=sysnumericmat, ncol=ncol(sysverifmat), nrow=nrow(sysverifmat))
dimnames(sysactualnumericmat) <- list(rownames(sysverifmat), colnames(sysverifmat))

# remove the zero sum rows
# sysactualnumericmat = sysactualnumericmat[rowSums(sysactualnumericmat[])>0,]
# actualnumericmat = actualnumericmat + 2

# mean to compare to
sysDEmean = rowMeans(sysactualnumericmat)
sysFoldchagen = sysactualnumericmat/sysDEmean
syslogFoldchange = log(sysFoldchagen)


#Read the cell line annotations to match the response ids
sangercellineannotations = read_csv("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Validation\ dataset\ from\ Sanger/Cell_listTue\ Mar\ 22\ 20_34_26\ 2022.csv")
cellinenamesforfoldchangetable = sangercellineannotations$Name[match(colnames(sysFoldchagen), sangercellineannotations$COSMIC_ID...1)]

# Run the enrichment and save results for each cell line

topgenes = "_hand_curated"

fGSEAhuman = function(results, filter="",samplename)  {
  pathways = pathways <- gmtPathways("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Curated\ gene\ signitures.gmt.txt")  #("msigdb.v7.0.symbols.gmt.txt")
  pathways = pathways[grepl(filter, names(pathways))]
  # str(head(pathways))
  fgseaRes <- suppressWarnings(fgsea(pathways, results, nperm=100000, minSize = 8, maxSize = 2000, BPPARAM=MulticoreParam(parallel::detectCores())))
  #  fgseaRes = fgseaRes[padj < 0.05, ]
  topPathwaysUp <- fgseaRes[ES > 0, ][head(order(pval), n=15), pathway]
  topPathwaysDown <- fgseaRes[ES < 0, ][head(order(pval), n=15), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  if (length(topPathways) > 0) {
    #' `#+ fig.width=10, fig.height=10`
    pdf(file = paste("DE_outputs_sanger_whole_regression_top",topgenes, "/fgsea/", samplename, ".pdf", sep = ""), width = 8, height = 11) # defaults to 7 x 7 inches
    plotGseaTable(pathways[topPathways], results, fgseaRes, gseaParam = 0.5, colwidths = c(7, 3, 0.8, 1.2, 1.2))
    knitr::kable(format(paste0("[", fgseaRes$pathway,"](", "http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",fgseaRes$pathway, ") ", fgseaRes$pval, " ", ifelse(fgseaRes$ES < 0, "downregulated", "upregulated")), caption = main))
    dev.off()
    write_delim(fgseaRes, paste("DE_outputs_sanger_whole_regression_top", topgenes, "/fgsea/", samplename, ".txt", sep = ""))
  } else {print(paste("out of",length(pathways),"no significant", filter, "pathways"))}
}


for(i in 1:length(cellinenamesforfoldchangetable)){
  print(cellinenamesforfoldchangetable[i])
  fGSEAhuman(syslogFoldchange[,i], samplename = cellinenamesforfoldchangetable[i])
}


# Read the data produced by enrichment analysis to produce roc curve

#read the sensitivity matrix
sensitivitytablesanger = read_xlsx("Validation dataset from Sanger/Drug_resistance_for_hand_cur.xlsx", sheet=2)
sensitivitytablesangermat = as.matrix(sensitivitytablesanger, rownames=T)
rownames(sensitivitytablesangermat) = sensitivitytablesanger$`Screened Compounds:`
sensitivitytablesangermat = sensitivitytablesangermat[,-1]
factorsensitivitytablesangermat = as.factor(sensitivitytablesangermat) 
dim(factorsensitivitytablesangermat) <- c(nrow(sensitivitytablesangermat), ncol(sensitivitytablesangermat))
dimnames(factorsensitivitytablesangermat) <- list(rownames(sensitivitytablesangermat), colnames(sensitivitytablesangermat))
factorsensitivitytablesangermat = factorsensitivitytablesangermat[, !duplicated(colnames(factorsensitivitytablesangermat))]


#drugnamealias = read_delim("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/Validation\ dataset\ from\ Sanger/listofdrugs.txt")

#read drug info from Sanger dataset
#sangerdrugs = read.delim(file = "Validation dataset from Sanger/listofdrugs.txt")
#sangerdrugs[sangerdrugs[,"drug_name"] %in% compounds$cmap_name,"drug_name"]

#lowersangerdrugs=NULL
#for(i in 1:length(sangerdrugs[,"drug_name"])){
#  lowersangerdrugs[i] = tolower(sangerdrugs[i,"drug_name"])
#}

#sangerdrugsforgmt = c(lowersangerdrugs, sangerdrugs[,"drug_name"])




# Renaming alias table with the annotations used before (some drugs strat in lowercase, im replacing the id with other name, bcs im too lazy to make a new column)
#for(drn in 1:length(drugnamealias$drug_name)){
#  tryCatch({
#    print(drugnamealias$drug_name[drn])
#    
#    drugnamealias$drug_id[drn] = sangerdrugsforgmt[str_detect(sangerdrugsforgmt, fixed(drugnamealias$drug_name[drn], ignore_case=TRUE))]
#  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#  
#}


for(drn in 1:ncol(factorsensitivitytablesangermat)){
  #print(colnames(factorsensitivitytablesangermat)[drn])
  
  print(strsplit(colnames(factorsensitivitytablesangermat)[drn], "[..]")[[1]][1])
  for(i in 1:length(drugnamealias$drug_id)){
    if(grepl(strsplit(colnames(factorsensitivitytablesangermat)[drn], "[..]")[[1]][1], drugnamealias$drug_name[i], fixed = TRUE) || grepl(strsplit(colnames(factorsensitivitytablesangermat)[drn], "[..]")[[1]][1], drugnamealias$synonyms[i], fixed = TRUE)){
      colnames(factorsensitivitytablesangermat)[drn] = drugnamealias$drug_id[i]
      break
    }
  }
  
  
}


#read the predictions
predictionmatrix = matrix(,nrow=length(cellinenamesforfoldchangetable),ncol=length(colnames(factorsensitivitytablesangermat)))
rownames(predictionmatrix) = cellinenamesforfoldchangetable
colnames(predictionmatrix) = colnames(factorsensitivitytablesangermat)

for(enrichedcellline in cellinenamesforfoldchangetable){
  tryCatch({
    print(enrichedcellline)
    enrichment = read.table(paste("DE_outputs_sanger_whole_regression_top", topgenes, "/fgsea/",enrichedcellline, ".txt", sep = ""), header=TRUE, sep=" ", row.names = 1)
    
    # the following two lines take only the signiture for each drug (up or down) that was enriched with a lower p value
      if(enrichment["erlotinib_down",]$padj>enrichment["erlotinib_up",]$padj){
        enrichment["erlotinib_down",]$padj=111
      }else{
        enrichment["erlotinib_up",]$padj=111
      }
    
    if(enrichment["lapatinib_down",]$padj>enrichment["lapatinib_up",]$padj){
      enrichment["lapatinib_down",]$padj=111
    }else{
      enrichment["lapatinib_up",]$padj=111
    }
    
    if(enrichment["talazoparib_down",]$padj>enrichment["talazoparib_up",]$padj){
      enrichment["talazoparib_down",]$padj=111
    }else{
      enrichment["talazoparib_up",]$padj=111
    }
    
    enrichment = enrichment[enrichment$padj!=111,]
    #enrichment = enrichment[enrichment$ES,]
    # cut from here if unneccesary
   
    #taking only significant ones
    #enrichment = enrichment[enrichment$padj<0.05,]
    
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
truematrix = factorsensitivitytablesangermat[which(rownames(factorsensitivitytablesangermat) %in% rownames(predictionmatrix)),which(colnames(factorsensitivitytablesangermat) %in% colnames(predictionmatrix))]
predictmatrix = predictionmatrix[which(rownames(predictionmatrix) %in% rownames(factorsensitivitytablesangermat)),which(colnames(predictionmatrix) %in% colnames(factorsensitivitytablesangermat))]
predictmatrix = predictmatrix[rownames(truematrix),colnames(truematrix)]
levels(truematrix) <- c(0,1)

# it gives worse results with binarizing the prediction table
#binpredictmatrix = predictmatrix
#binpredictmatrix[which(binpredictmatrix>0.5)]=1
#binpredictmatrix[which(binpredictmatrix<0.5)]=0

roc(as.numeric(truematrix)-1, predictmatrix, plot=T)


for(roci in 1:dim(truematrix)[2]){
  tryCatch({
    pdf(file = paste("DE_outputs_sanger_whole_regression_top",topgenes, "/rocs_test_rem<0/", colnames(truematrix)[roci], ".pdf", sep = ""), width = 8, height = 11) # defaults to 7 x 7 inches
    roccurve = roc(truematrix[,roci], predictmatrix[,roci], plot=T)
    title(main = paste(colnames(truematrix)[roci],roccurve$auc, sep=" "))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

