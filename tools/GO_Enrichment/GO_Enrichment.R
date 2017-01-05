#############################################################################
##
## What: GO_Enrichment.R
## Who: 
## When: 
##
## Script to perform Gene Ontology enrichment analysis
##
## Args: 
## -----
## rData=DE_DESeq2.RData    # differentially expressed gene result from DESeq2
## log2Fold=2               # log2 foldchange parameter to filter out gene list
## pvalue=0.01              # pvalue parameter to filter out gene list
## orgDb=org.Hs.eg.db       # genome wide annotation database for human
## univ=all                 # universe to enrich with DE genes
## type=gene_name           # type of key to be converted into Entrez ID
## plotCategory=5           # number of category to be shown in the barplot
## out=GO_Analysis          # output directory
##
##
############################################################################
options(stringsAsFactors=FALSE)
library(clusterProfiler)
library(DESeq2)

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
  
}

args <- commandArgs(T)
rData <- parseArgs(args,"rData=") # result DESeq2.RData
log2Fold <- parseArgs(args, "log2Fold=", 2, "as.numeric") # log2 foldchange parameter for filtering out gene list
pvalue <- parseArgs(args, "pvalue=", 0.01, "as.numeric") # pvalue parameter for filtering gene list
orgDb <- parseArgs(args, "orgDb=", "org.Hs.eg.db") # genome wide annotation for human 
univ <- parseArgs(args, "univ=", "all") # universe to compare with DE genes
type <- parseArgs(args, "type=", "gene_name") # type of key to be converted into entrez id
plotCategory <- parseArgs(args, "plotCategory=", 5, "as.numeric") # number of category to be shown in the barplot
out <- parseArgs(args,"out=", "GO_Analysis") # output directory

runstr <- "Call with: Rscript goEnrichment.R [rData=DESeq2.RData] [log2Fold=2] [pvalue=0.01] [orgDb=org.Hs.eg.db] 
[univ=all|express] [type=gene_name|ensembl_id] [plotCategory=5] [out=GO_Analysis] "
if (!is.numeric(log2Fold)) stop("log2 foldchange not numeric. Run with:\n",runstr)
if (!is.numeric(pvalue)) stop("pvalue not numeric. Run with:\n", runstr)
if (!is.numeric(plotCategory)) stop("plotCategory not numeric. Run with:\n",runstr)
if (!require(orgDb, character.only=TRUE)) stop ("Annotation DB", orgDb, " not installed\n")
if(!file.exists(out)) dir.create(out) 


runstr <- paste("Call with: Rscript goEnrichment.R [rData=",rData,"] [log2Fold=",log2Fold, "] [pvalue=",pvalue,"] 
                [orgDb=",orgDb,"] [univ=",univ,"] [plotCategory=",plotCategory,"] [type=",type,"] [out=",out,"]")
cat(runstr)

##
# gene ontology enrichment and plots
##
load(rData) # load the RData
for (i in 1:length(res)) { # looping through every element in the list
  resultData <- as.data.frame(res[[i]], row.names = NULL)
  collect.padj <- grep("padj", colnames(resultData)) # grep column number for padj
  resultData<- resultData[-which(is.na(resultData[ ,collect.padj])), ] # remove all rows with NA in padj column 

  keyType <- grep(type, colnames(resultData)) # grep the column number for gene name or ensembl id
  collect.fc <- grep("log2FoldChange", colnames(resultData)) # grep the column number for log2 foldchange
  collect.pvalue <- grep("pvalue", colnames(resultData)) # grep the column number for pvalues
  DE <- resultData[ ,collect.pvalue] < pvalue & abs(resultData[ ,collect.fc]) > log2Fold
  geneFiltered <- (resultData[ ,keyType][DE]) # list of all filtered genes
  entrezID <- bitr(geneFiltered, fromType ="SYMBOL", toType ="ENTREZID", OrgDb=orgDb) # convert gene name into entrez id
  collect.expression <- grep("baseMean", colnames(resultData)) # grep the column for expression value
  express <- resultData[ ,collect.expression] > 0 
  geneExpress <- (resultData[ ,keyType][express]) # list genes with expression value > 0

  if(type=="gene_name") {
      fType <- "SYMBOL"
  } else {
      fType <- "ENSEMBL"
  }
  # convert gene name / ensembl id (with expression value >0 for the universe) into entrez id
  entrezUnivID <- bitr(geneExpress, fromType =fType, toType ="ENTREZID", OrgDb=orgDb) 
 

  if (univ=="all") {
    uni <- orgDb 
  } else {
    uni <- entrezUnivID[[2]]
  }
  enrich <- enrichGO(entrezID[[2]], universe = uni, OrgDb = orgDb, keytype = "ENTREZID",ont="BP", readable = TRUE)
  
  if(!(nrow(enrich)==0)){

  # write GO enrichment table into output file 
  write.csv(enrich, file=paste0(out, "/", gsub("=.*", "", conts[i,1]), "_GO_Enrichment.csv"), row.names=T)
  
  # create barplot showing GO category
  png(file=paste0(out, "/", gsub("=.*", "", conts[i,1]), "_GO_Barplot.png"), width = 700, height = 500)
  plot(barplot(enrich, showCategory = 10))
  dev.off()

  # create network plot for the results
  png(file=paste0(out, "/", gsub("=.*", "", conts[i,1]), "_GO_network.png"), width = 700, height = 500)
  plot(enrichMap(enrich))
  dev.off()
  } else {
      errorstr <- "No GO enrichment result available"
      cat(errorstr)
  }
}

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GO_Enrichment_session_info.txt", sep=""))
save(entrezUnivID, enrich, file=paste0(out, "/GO_Enrichment.RData"))

