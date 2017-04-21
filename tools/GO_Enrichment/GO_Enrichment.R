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
## padj=0.01                # padj parameter to filter out gene list
## orgDb=org.Hs.eg.db       # genome wide annotation database for human
## organism="human"         # organism  	 	
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
log2Fold <- parseArgs(args, "log2Fold=", 0, "as.numeric") # log2 foldchange parameter for filtering out gene list
padj     <- parseArgs(args, "padj=", 0.01, "as.numeric") # padj parameter for filtering gene list
orgDb    <- parseArgs(args, "orgDb=", "org.Hs.eg.db") # genome wide annotation for human
org      <- parseArgs(args, "organism=", "human") # organism  	 	
univ     <- parseArgs(args, "univ=", "all") # universe to compare with DE genes
type     <- parseArgs(args, "type=", "gene_name") # type of key to be converted into entrez id
plotCategory <- parseArgs(args, "plotCategory=", 10, "as.numeric") # number of category to be shown in the barplot
out      <- parseArgs(args,"out=", "GO_Analysis") # output directory
cores    <- parseArgs(args, "cores=", 1, "as.numeric") # number of category to be shown in the barplot

runstr <- paste("Call with: Rscript goEnrichment.R rData=<DESeq2.RData> [log2Fold=2]",
                "[padj=0.01] [orgDb=org.Hs.eg.db] [org=human] [univ=all|express] [type=gene_name|ensembl_id]",
                "[plotCategory=10] [out=GO_Analysis] [cores=1]")
if (!is.numeric(log2Fold)) stop("log2 foldchange not numeric. Run with:\n",runstr)
if (!is.numeric(padj)) stop("padj not numeric. Run with:\n", runstr)
if (!is.numeric(plotCategory)) stop("plotCategory not numeric. Run with:\n",runstr)
if (!require(orgDb, character.only=TRUE)) stop ("Annotation DB", orgDb, " not installed\n")
if (!file.exists(out)) dir.create(out) 

runstr <- paste0("Called with: Rscript goEnrichment.R rData=", rData, " log2Fold=", log2Fold, " padj=", padj,
                 " orgDb=", orgDb," org=", org," univ=", univ," plotCategory=", plotCategory," type=", type," out=", out, " cores=", cores)
cat(runstr)


##
# gene ontology enrichment and plots
##
load(rData) # load the RData
res <- lapply(res, as.data.frame, row.names=NULL)

calculateGoEnrichment <-  function(x) {

    library(clusterProfiler)
    library(ReactomePA)

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
    
    # Separate enrichment analysis of GO terms and pathways for up- and downregulated genes
    up <- resultData[, colpadj] < padj & (resultData[, colfc]) > log2Fold
    down <- resultData[, colpadj] < padj & (resultData[, colfc]) < log2Fold


    # GO and Pathway enrichment analysis for the upregulated genes
    if(sum(up) == 0) {
        writeLines("No upregulated genes found!", paste0(out, "/", contrast, "_GO_Enrichment_upregulated_genes.csv"))
        return()
    }

    # convert to entrezID DE upregulated/univers genes
    getEntrezId <- function(genes) {
        bitr(genes, fromType=if(type == "gene_name") "SYMBOL" else "ENSEMBL", toType="ENTREZID", OrgDb=orgDb)
    }
    entrezDeId   <- getEntrezId(resultData[up, genes])
    entrezUnivId <- getEntrezId(resultData[  , genes])
    
    enriched <- enrichGO(entrezDeId[[2]], OrgDb=orgDb, keytype="ENTREZID", ont="BP", readable=TRUE,
                         universe=if(univ == "all") orgDb else entrezUnivId[[2]])
    
    enrichedDAVID <- enrichDAVID(entrezDeId[[2]], idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_ALL", 
                                 david.user = "g.petrosino@imb-mainz.de")

    enrichedKEGG <- enrichKEGG(entrezDeId[[2]], org, universe = entrezUnivId[[2]])

    enrichedReactome <- enrichPathway(entrezDeId[[2]], org, readable=TRUE, 
									 universe = entrezUnivId[[2]])

    # write GO and Pathway enrichment tables into output file 
    write.csv(enriched, file=paste0(out, "/", contrast, "_GO_Enrichment_upregulated_genes.csv"), row.names=T)
    write.csv(enrichedDAVID, file=paste0(out, "/", contrast, "_DAVID_GO_Enrichment_upregulated_genes.csv"), row.names=T)
    write.csv(enrichedKEGG, file=paste0(out, "/", contrast, "_KEGG_Pathway_Enrichment_upregulated_genes.csv"), row.names=T)
    write.csv(enrichedReactome, file=paste0(out, "/", contrast, "_Reactome_Pathway_Enrichment_upregulated_genes.csv"), row.names=T)

    if(nrow(enriched) > 0) {
        # create barplot showing GO category
        png(file=paste0(out, "/", contrast, "_GO_Barplot_upregulated_genes.png"), width=700, height=500)
        plot(barplot(enriched, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_GO_network_upregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enriched))
        dev.off()
    }

    if(nrow(enrichedDAVID) > 0) {
        # create barplot showing GO category
        png(file=paste0(out, "/", contrast, "_DAVID_GO_Barplot_upregulated_genes.png"), width=700, height=500)
        plot(barplot(enrichedDAVID, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_DAVID_GO_network_upregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enrichedDAVID))
        dev.off()
    }

    if(nrow(enrichedKEGG) > 0) {
        # create barplot showing Pathway terms
        png(file=paste0(out, "/", contrast, "_KEGG_Barplot_upregulated_genes.png"), width=700, height=500)
        plot(barplot(enrichedKEGG, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_KEGG_network_upregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enrichedKEGG))
        dev.off()
    }
    
    if(nrow(enrichedReactome) > 0) {
        # create barplot showing Pathway terms
        png(file=paste0(out, "/", contrast, "_Reactome_Barplot_upregulated_genes.png"), width=700, height=500)
        plot(barplot(enrichedReactome, showCategory=plotCategory))
        dev.off()
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_Reactome_network_upregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enrichedReactome))
        dev.off()
    }



    # GO and Pathway enrichment analysis for the downregulated genes
    if(sum(down) == 0) {
        writeLines("No downregulated genes found!", paste0(out, "/", contrast, "_GO_Enrichment_downregulated_genes.csv"))
        return()
    }
    # convert to entrezID downregulated/univers genes
    getEntrezId <- function(genes) {
        bitr(genes, fromType=if(type == "gene_name") "SYMBOL" else "ENSEMBL", toType="ENTREZID", OrgDb=orgDb)
    }
    entrezDeId   <- getEntrezId(resultData[down, genes])
    entrezUnivId <- getEntrezId(resultData[  , genes])
    
    enriched <- enrichGO(entrezDeId[[2]], OrgDb=orgDb, keytype="ENTREZID", ont="BP", readable=TRUE,
                         universe=if(univ == "all") orgDb else entrezUnivId[[2]])

    enrichedDAVID <- enrichDAVID(entrezDeId[[2]], idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_ALL", 
                                 david.user = "g.petrosino@imb-mainz.de")
    
    enrichedKEGG <- enrichKEGG(entrezDeId[[2]], org, universe = entrezUnivId[[2]])

    enrichedReactome <- enrichPathway(entrezDeId[[2]], org, readable=TRUE, 
									 universe = entrezUnivId[[2]])

    # write GO and Pathway enrichment tables into output file 
    write.csv(enriched, file=paste0(out, "/", contrast, "_GO_Enrichment_downregulated_genes.csv"), row.names=T)
    write.csv(enrichedDAVID, file=paste0(out, "/", contrast, "_DAVID_GO_Enrichment_downregulated_genes.csv"), row.names=T)
    write.csv(enrichedKEGG, file=paste0(out, "/", contrast, "_KEGG_Pathway_Enrichment_downregulated_genes.csv"), row.names=T)
    write.csv(enrichedReactome, file=paste0(out, "/", contrast, "_Reactome_Pathway_Enrichment_downregulated_genes.csv"), row.names=T)
  
    if(nrow(enriched) > 0) {
        # create barplot showing GO category
        png(file=paste0(out, "/", contrast, "_GO_Barplot_downregulated_genes.png"), width=700, height=500)
        plot(barplot(enriched, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_GO_network_downregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enriched))
        dev.off()
    }

    if(nrow(enrichedDAVID) > 0) {
        # create barplot showing GO category
        png(file=paste0(out, "/", contrast, "_DAVID_GO_Barplot_downregulated_genes.png"), width=700, height=500)
        plot(barplot(enrichedDAVID, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_DAVID_GO_network_downregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enrichedDAVID))
        dev.off()
    }

    if(nrow(enrichedKEGG) > 0) {
        # create barplot showing Pathway terms
        png(file=paste0(out, "/", contrast, "_KEGG_Barplot_downregulated_genes.png"), width=700, height=500)
        plot(barplot(enrichedKEGG, showCategory=plotCategory))
        dev.off()
  
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_KEGG_network_downregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enrichedKEGG))
        dev.off()
    }
    
    if(nrow(enrichedReactome) > 0) {
        # create barplot showing Pathway terms
        png(file=paste0(out, "/", contrast, "_Reactome_Barplot_downregulated_genes.png"), width=700, height=500)
        plot(barplot(enrichedReactome, showCategory=plotCategory))
        dev.off()
        # create network plot for the results
        png(file=paste0(out, "/", contrast, "_Reactome_network_downregulated_genes.png"), width=700, height=500)
        plot(enrichMap(enrichedReactome))
        dev.off()
    }

}


if(cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, c("log2Fold", "padj", "orgDb", "org", "univ", "type", "plotCategory", "out"))
    parLapply(cl, zipup(res, names(res)), calculateGoEnrichment)
    stopCluster(cl)
} else {
    lapply(zipup(res, names(res)), calculateGoEnrichment)
}

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GO_Enrichment_session_info.txt", sep=""))
