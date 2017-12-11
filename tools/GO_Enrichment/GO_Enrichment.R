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
## organism="human"         # organism. Currently, one of human, mouse, fly, worm, yeast
## univ=all                 # universe to enrich with DE genes
## type=gene_name           # type of key to be converted into Entrez ID
## plotCategory=10          # number of categories to be shown in the barplot
## out=GO_Analysis          # output directory
## cores=1                  # number of cores to use
##
############################################################################
options(stringsAsFactors=FALSE)
library(DESeq2)
library(parallel)

# suported organisms
orgDb <- c(human="org.Hs.eg.db",
           mouse="org.Mm.eg.db",
           fly  ="org.Dm.eg.db",
           worm ="org.Ce.eg.db",
           yeast="org.Sc.sgd.db")

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  if(length(i <- grep(string,args,fixed=T))== 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args     <- commandArgs(T)
rData    <- parseArgs(args, "rData=") # result DESeq2.RData
log2Fold <- parseArgs(args, "log2Fold=", 0, "as.numeric") # log2 foldchange parameter for filtering out gene list
padj     <- parseArgs(args, "padj=", 0.01, "as.numeric") # padj parameter for filtering gene list
org      <- parseArgs(args, "organism=", "human") # organism  	 	
univ     <- parseArgs(args, "univ=", "all") # universe to compare with DE genes
type     <- parseArgs(args, "type=", "gene_name") # type of key to be converted into entrez id
plotCategory <- parseArgs(args, "plotCategory=", 10, "as.numeric") # number of category to be shown in the barplot
out      <- parseArgs(args,"out=", "GO_Analysis") # output directory
cores    <- parseArgs(args, "cores=", 1, "as.numeric") # number of category to be shown in the barplot

runstr <- paste("Call with: Rscript goEnrichment.R rData=<DESeq2.RData> [log2Fold=2]",
                "[padj=0.01] [org=human] [univ=all|expressed] [type=gene_name|ensembl_id]",
                "[plotCategory=10] [out=GO_Analysis] [cores=1]")
if (!is.numeric(log2Fold)) stop("log2 foldchange not numeric. Run with:\n",runstr)
if (!is.numeric(padj)) stop("padj not numeric. Run with:\n", runstr)
if (!is.numeric(plotCategory)) stop("plotCategory not numeric. Run with:\n",runstr)
if (!(univ %in% c("all", "expressed"))) stop("univ must be one of 'all|expressed'")
if (!(org %in% names(orgDb))) stop("Organism must be one of", paste(names(orgDb), collapse=", "))
if (!require(orgDb[org], character.only=TRUE)) stop ("Annotation DB", orgDb[org], " not installed\n")
if (!file.exists(out)) dir.create(out) 

runstr <- paste0("Called with: Rscript goEnrichment.R rData=", rData, " log2Fold=", log2Fold, " padj=", padj,
                 " org=", org," univ=", univ," plotCategory=", plotCategory," type=", type," out=", out, " cores=", cores)
cat(runstr)

##
# gene ontology enrichment and plots
##
load(rData) # load the RData
res <- lapply(res, as.data.frame, row.names=NULL)

processContrast <-  function(x) {

    library(clusterProfiler)
    library(ReactomePA)
    library(Cairo)

    calculateGoEnrichment <-  function(de.genes, univ.genes, suffix) {
        # convert to entrezID downregulated/univers genes
        getEntrezId <- function(genes) {
            bitr(genes, fromType=if(type == "gene_name") "SYMBOL" else "ENSEMBL", toType="ENTREZID", OrgDb=orgDb[org])
        }
        entrezDeId   <- getEntrezId(de.genes)
        entrezUnivId <- getEntrezId(univ.genes)
        
        enriched         <- enrichGO(entrezDeId$ENTREZID, OrgDb=orgDb[org], keyType="ENTREZID", ont="BP", readable=TRUE,
                                     universe=if(univ == "all") orgDb[org] else entrezUnivId$ENTREZID)
        enrichedKEGG     <- enrichKEGG(entrezDeId$ENTREZID, org, universe=entrezUnivId$ENTREZID)
        enrichedReactome <- enrichPathway(entrezDeId$ENTREZID, org, readable=TRUE, universe=entrezUnivId$ENTREZID)

        # write GO and Pathway enrichment tables into output file 
        write.csv(as.data.frame(enriched),
                  file=paste0(out, "/", contrast, "_GO_Enrichment_", suffix, "_genes.csv"))
        write.csv(as.data.frame(enrichedKEGG),
                  file=paste0(out, "/", contrast, "_KEGG_Pathway_Enrichment_", suffix, "_genes.csv"))
        write.csv(as.data.frame(enrichedReactome),
                  file=paste0(out, "/", contrast, "_Reactome_Pathway_Enrichment_", suffix, "_genes.csv"))
      
        if(!is.null(enriched) && nrow(enriched) > 0) {
            # create barplot showing GO category
            CairoPNG(file=paste0(out, "/", contrast, "_GO_Barplot_", suffix, "_genes.png"), width=700, height=500)
            print(barplot(enriched, showCategory=plotCategory))
            dev.off()
      
            # create network plot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_GO_network_", suffix, "_genes.png"), width=700, height=500)
            print(enrichMap(enriched))
            dev.off()
        }

        if(!is.null(enrichedKEGG) && nrow(enrichedKEGG)> 0) {
            # create barplot showing Pathway terms
            CairoPNG(file=paste0(out, "/", contrast, "_KEGG_Barplot_", suffix, "_genes.png"), width=700, height=500)
            print(barplot(enrichedKEGG, showCategory=plotCategory))
            dev.off()
      
            # create network plot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_KEGG_network_", suffix, "_genes.png"), width=700, height=500)
            print(enrichMap(enrichedKEGG))
            dev.off()
        }
        
        if(!is.null(enrichedReactome) && nrow(enrichedReactome) > 0) {
            # create barplot showing Pathway terms
            CairoPNG(file=paste0(out, "/", contrast, "_Reactome_Barplot_", suffix, "_genes.png"), width=700, height=500)
            print(barplot(enrichedReactome, showCategory=plotCategory))
            dev.off()
            # create network plot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_Reactome_network_", suffix, "_genes.png"), width=700, height=500)
            print(enrichMap(enrichedReactome))
            dev.off()
        }
    }

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
    up   <- resultData[, colpadj] < padj & (resultData[, colfc]) > log2Fold
    down <- resultData[, colpadj] < padj & (resultData[, colfc]) < log2Fold

    # GO and Pathway enrichment analysis for the upregulated genes
    if(sum(up) == 0)
        writeLines("No upregulated genes found!", paste0(out, "/", contrast, "_GO_Enrichment_upregulated_genes.csv"))
    else
        calculateGoEnrichment(resultData[up, genes], resultData[, genes], "up")

    # GO and Pathway enrichment analysis for the downregulated genes
    if(sum(down) == 0)
        writeLines("No downregulated genes found!", paste0(out, "/", contrast, "_GO_Enrichment_downregulated_genes.csv"))
    else
        calculateGoEnrichment(resultData[down, genes], resultData[, genes], "down")

    invisible(0)
}

if(cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, c("log2Fold", "padj", "orgDb", "org", "univ", "type", "plotCategory", "out"))
    parLapply(cl, zipup(res, names(res)), processContrast)
    stopCluster(cl)
} else {
    lapply(zipup(res, names(res)), processContrast)
}

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GO_Enrichment_session_info.txt", sep=""))
