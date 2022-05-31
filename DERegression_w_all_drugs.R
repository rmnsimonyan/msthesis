# training on CCLE validating on GDSC (Sanger dataset)

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

# Running comprehensive Regression on the whole thing (all drug-treated cell-lines)

#reading metadata of LINCS dataset - out training dataset
metag = read.delim(file = "Level3L1000/instinfo_beta.txt", row.names = "sample_id")
compounds = read.delim(file = "Level3L1000/compoundinfo_beta.txt")
genes = read.delim(file = "Level3L1000/geneinfo_beta.txt", row.names = 1)

#reading the ids of control plates to later filter only the ones that we need (according to the compound) and read only them
normnames = read_gctx_meta("Level3L1000/level3_beta_ctl_n188708x12328.gctx", dim = "column")
normcids = unlist(strsplit(unlist(normnames), ":"))[ c(TRUE,FALSE) ]


#read drug info from Sanger dataset
sangerdrugs = read.delim(file = "Validation dataset from Sanger/listofdrugs.txt")
sangerdrugs[sangerdrugs[,"drug_name"] %in% compounds$cmap_name,"drug_name"]

lowersangerdrugs=NULL
for(i in 1:length(sangerdrugs[,"drug_name"])){
  lowersangerdrugs[i] = tolower(sangerdrugs[i,"drug_name"])
}

sangerdrugsforgmt = c(lowersangerdrugs[which(lowersangerdrugs %in% compounds$cmap_name)], sangerdrugs[sangerdrugs[,"drug_name"] %in% compounds$cmap_name,"drug_name"])


#gcpnames = read_gctx_meta("Level3L1000/level3_beta_trt_cp_n1805898x12328.gctx", dim = "column")
#gcpcids = unlist(strsplit(unlist(gcpnames), ":"))[ c(TRUE,FALSE) ]


gcp <- parse_gctx("Level3L1000/level3_beta_trt_cp_n1805898x12328.gctx", cid= rownames(metag[which(metag$cmap_name %in% sangerdrugsforgmt),]))
cmpdids = unlist(strsplit(gcp@cid, ":"))[ c(TRUE,FALSE) ]
gcpcelllines = metag$cell_iname[match(cmpdids, metag$det_plate)]
sangdrugnames = metag$cmap_name[match(gcp@cid, rownames(metag))]

norm = parse_gctx("Level3L1000/level3_beta_ctl_n188708x12328.gctx", cid = which(normcids %in% unique(metag$det_plate[which(metag$cell_iname %in% metag$cell_iname[which(metag$det_plate %in% normcids[which(normcids %in% cmpdids)])])])))
nids = unlist(strsplit(unlist(norm@cid), ":"))[ c(TRUE,FALSE) ]
normcelllines = metag$cell_iname[match(nids, metag$det_plate)]

mat = cbind(gcp@mat, norm@mat)#[,randomids])
mat = mat[rownames(genes[which(genes$feature_space=="landmark"),]),]

group = c(sangdrugnames, rep("1norm", length(nids)))
cellline = c(gcpcelllines, normcelllines)

mm <- model.matrix(~cellline + group)

fit <- lmFit(mat, mm)

head(coef(fit))
fit <- eBayes(fit)

for(i in colnames(fit$coefficients)[grepl("group", colnames(fit$coefficients))]){
  
  coefname = unlist(strsplit(i, "group"))[2]
  print(coefname)
  top.table <- topTable(fit, coef = i,  n = Inf)
  
  rownames(top.table) = genes[rownames(top.table),"gene_symbol"]
  write_delim(cbind(rownames(top.table),top.table), paste("DE_outputs_sanger_whole_regression/", coefname, ".txt", sep = ""))
  #head(top.table, 20)
  down = as.data.frame(t(c(paste(coefname, "_down", sep = ""), "", rownames(top.table)[top.table$adj.P.Val<0.01 & top.table$logFC<0])))
  names(down) = NULL
  write_delim(down, paste("DE_outputs_sanger_whole_regression/",coefname, "_down.txt", sep = ""), col_names=F)
  
  up = as.data.frame(t(c(paste(coefname, "_up", sep = ""), "", rownames(top.table)[top.table$adj.P.Val<0.01 & top.table$logFC>0])))
  names(up) = NULL
  write_delim(up, paste("DE_outputs_sanger_whole_regression/",coefname, "_up.txt", sep = ""), col_names=F)
}


# filtering out 50 (or n) top genes here for each signiture
topgenes = 5

for(drugname in sangerdrugsforgmt){
  tryCatch({
    print(drugname)
    
    DEgenetable = read.table(paste("DE_outputs_sanger_whole_regression/",drugname, ".txt", sep = ""), header=TRUE, sep=" ", row.names = 1)
    DEgenetable = DEgenetable[DEgenetable$P.Value<0.01,]
    #    DEgenetable = DEgenetable[DEgenetable$logFC>0,]
    
    up = as.data.frame(t(c(paste(drugname, "_up", sep = ""), "", rownames(DEgenetable[DEgenetable$logFC>0,])[1:topgenes][!is.na(rownames(DEgenetable[DEgenetable$logFC>0,])[1:topgenes])])))
    names(up) = NULL
    write_delim(up, paste("DE_outputs_sanger_whole_regression_top", topgenes, "/",drugname, "_up.txt", sep = ""), col_names=F)
    
    
    down = as.data.frame(t(c(paste(drugname, "_down", sep = ""), "", rownames(DEgenetable[DEgenetable$logFC<0,])[1:topgenes][!is.na(rownames(DEgenetable[DEgenetable$logFC<0,])[1:topgenes])])))
    names(down) = NULL
    write_delim(down, paste("DE_outputs_sanger_whole_regression_top", topgenes, "/",drugname, "_down.txt", sep = ""), col_names=F)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


# Then run the following in termnial to produce the combined gmt file for further analysis in fGSEA
# cat $(ls | grep "up\|down")> combined.gmt.txt

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

listof1000genes = read.delim("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/DE_outputs/abiraterone.txt", sep = " ")

# Run the enrichment and save results for each cell line

fGSEAhuman = function(results, filter="",samplename)  {
  pathways = pathways <- gmtPathways(paste("/Users/armansimonyan/Documents/Copenhagen/Rigshospitalet/TaskX.\ Thesis/DESEQ2/DE_outputs_sanger_whole_regression_top", topgenes, "/combined.gmt.txt", sep=""))  #("msigdb.v7.0.symbols.gmt.txt")
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
  fGSEAhuman(syslogFoldchange[rownames(syslogFoldchange) %in% listof1000genes[,1],i], samplename = cellinenamesforfoldchangetable[i])
}


# Read the data produced by enrichment analysis to produce roc curve

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

# Renaming alias table with the annotations used before (some drugs strat in lowercase, im replacing the id with other name, bcs im too lazy to make a new column)
for(drn in 1:length(drugnamealias$drug_name)){
  tryCatch({
    print(drugnamealias$drug_name[drn])
    
    drugnamealias$drug_id[drn] = sangerdrugsforgmt[str_detect(sangerdrugsforgmt, fixed(drugnamealias$drug_name[drn], ignore_case=TRUE))]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


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
predictionmatrix = matrix(,nrow=length(cellinenamesforfoldchangetable),ncol=length(unique(sangerdrugsforgmt)))
rownames(predictionmatrix) = cellinenamesforfoldchangetable
colnames(predictionmatrix) = unique(sangerdrugsforgmt)

for(enrichedcellline in cellinenamesforfoldchangetable){
  tryCatch({
    print(enrichedcellline)
    enrichment = read.table(paste("DE_outputs_sanger_whole_regression_top", topgenes, "/fgsea/",enrichedcellline, ".txt", sep = ""), header=TRUE, sep=" ", row.names = 1)
    
    # the following two lines take only the signiture for each drug (up or down) that was enriched with a lower p value
    for(j in 1:(dim(enrichment)[1]/2)){
      if(enrichment[j*2-1,]$padj>enrichment[j*2,]$padj){
        enrichment[j*2-1,]$padj=111
      }else{
        enrichment[j*2,]$padj=111
      }
    }
    enrichment = enrichment[enrichment$padj!=111,]
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
truematrix = factorsensitivitytablesangermat[which(rownames(factorsensitivitytablesangermat) %in% rownames(predictionmatrix)),which(colnames(factorsensitivitytablesangermat) %in% colnames(predictionmatrix))]
predictmatrix = predictionmatrix[which(rownames(predictionmatrix) %in% rownames(factorsensitivitytablesangermat)),which(colnames(predictionmatrix) %in% colnames(factorsensitivitytablesangermat))]
predictmatrix = predictmatrix[rownames(truematrix),colnames(truematrix)]
levels(truematrix) <- c(0,1)

# it gives worse results with binarizing the prediction table
#binpredictmatrix = predictmatrix
#binpredictmatrix[which(binpredictmatrix>0.5)]=1
#binpredictmatrix[which(binpredictmatrix<0.5)]=0

roc(truematrix, predictmatrix, plot=T)


for(roci in 1:dim(truematrix)[2]){
  tryCatch({
  pdf(file = paste("DE_outputs_sanger_whole_regression_top",topgenes, "/rocs/", colnames(truematrix)[roci], ".pdf", sep = ""), width = 8, height = 11) # defaults to 7 x 7 inches
  roccurve = roc(truematrix[,roci], predictmatrix[,roci], plot=T)
  title(main = paste(colnames(truematrix)[roci],roccurve$auc, sep=" "))
  dev.off()
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



