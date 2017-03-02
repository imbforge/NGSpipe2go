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
## plotCategory=10          # number of category to be shown in the barplot
## out=GO_Analysis          # output directory
## cores=1                  # number of cores to use
##
############################################################################
options(stringsAsFactors=FALSE)
library(DESeq2)
library(parallel)

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  if(length(i <- grep(string,args,fixed=T))== 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args     <- commandArgs(T)
rData    <- parseArgs(args,"rData=") # result DESeq2.RData
log2Fold <- parseArgs(args, "log2Fold=", 2, "as.numeric") # log2 foldchange parameter for filtering out gene list
pvalue   <- parseArgs(args, "pvalue=", 0.01, "as.numeric") # pvalue parameter for filtering gene list
orgDb    <- parseArgs(args, "orgDb=", "org.Hs.eg.db") # genome wide annotation for human 
univ     <- parseArgs(args, "univ=", "all") # universe to compare with DE genes
type     <- parseArgs(args, "type=", "gene_name") # type of key to be converted into entrez id
plotCategory <- parseArgs(args, "plotCategory=", 10, "as.numeric") # number of category to be shown in the barplot
out      <- parseArgs(args,"out=", "GO_Analysis") # output directory
cores    <- parseArgs(args, "cores=", 1, "as.numeric") # number of category to be shown in the barplot

runstr <- paste("Call with: Rscript goEnrichment.R rData=<DESeq2.RData> [log2Fold=2]",
                "[pvalue=0.01] [orgDb=org.Hs.eg.db] [univ=all|express] [type=gene_name|ensembl_id]",
                "[plotCategory=10] [out=GO_Analysis] [cores=1]")
if (!is.numeric(log2Fold)) stop("log2 foldchange not numeric. Run with:\n",runstr)
if (!is.numeric(pvalue)) stop("pvalue not numeric. Run with:\n", runstr)
if (!is.numeric(plotCategory)) stop("plotCategory not numeric. Run with:\n",runstr)
if (!require(orgDb, character.only=TRUE)) stop ("Annotation DB", orgDb, " not installed\n")
if(!file.exists(out)) dir.create(out) 

runstr <- paste0("Called with: Rscript goEnrichment.R rData=", rData, " log2Fold=", log2Fold, " pvalue=", pvalue,
                 " orgDb=", orgDb," univ=", univ," plotCategory=", plotCategory," type=", type," out=", out, " cores=", cores)
cat(runstr)

##
# gene ontology enrichment and plots
##
load(rData) # load the RData
res <- lapply(res, as.data.frame, row.names=NULL)

calculateGoEnrichment <-  function(x) {

    library(clusterProfiler)

    # unzip res/contrast
    resultData <- x[[1]]
    contrast   <- x[[2]]

    # grep where the interesting columns are...
    colpadj   <- grep("padj", colnames(resultData))
    genes     <- grep(type, colnames(resultData))
    colfc     <- grep("log2FoldChange", colnames(resultData))
    colpvalue <- grep("pvalue", colnames(resultData))
    colexpression <- grep("baseMean", colnames(resultData))

    # remove all genes that were not tested (pval=NA), and keep only those that
    # where expressed and pass the FDR and FC thresholds
    resultData <- resultData[-which(is.na(resultData[, colpadj])), ]
    resultData <- resultData[resultData[, colexpression] > 0, ]
    DE <- resultData[, colpvalue] < pvalue & abs(resultData[, colfc]) > log2Fold

    if(sum(DE) == 0) {
        writeLines("No DE genes found!", paste0(out, "/", contrast, "_GO_Enrichment.csv"))
        return()
    }

    # convert to entrezID DE/univers genes
    getEntrezId <- function(genes) {
        bitr(genes, fromType=if(type == "gene_name") "SYMBOL" else "ENSEMBL", toType="ENTREZID", OrgDb=orgDb)
    }
    entrezDeId   <- getEntrezId(resultData[DE, genes])
    entrezUnivId <- getEntrezId(resultData[  , genes])
    
    enriched <- enrichGO(entrezDeId[[2]], OrgDb=orgDb, keytype="ENTREZID", ont="BP", readable=TRUE,
                         universe=if(univ == "all") orgDb else entrezUnivId[[2]])
    
    # write GO enrichment table into output file 
    write.csv(enriched, file=paste0(out, "/", contrast, "_GO_Enrichment.csv"), row.names=T)
  
    if(nrow(enriched) > 0) {
        # create barplot showing GO category
        png(file=paste0(out, "/", contrast, "_GO_Barplot.png"), width=700, height=500)
        plot(barplot(enriched, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_GO_network.png"), width=700, height=500)
        plot(enrichMap(enriched))
        dev.off()
    }
}

if(cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, c("log2Fold", "pvalue", "orgDb", "univ", "type", "plotCategory", "out"))
    parLapply(cl, zipup(res, names(res)), calculateGoEnrichment)
    stopCluster(cl)
} else {
    lapply(zipup(res, names(res)), calculateGoEnrichment)
}

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GO_Enrichment_session_info.txt", sep=""))

