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
library(clusterProfiler)
library(ReactomePA)
library(Cairo)
library(ggplot2)

# supported organisms
orgDb <- c(human="org.Hs.eg.db",
           mouse="org.Mm.eg.db",
           fly  ="org.Dm.eg.db",
           worm ="org.Ce.eg.db",
           yeast="org.Sc.sgd.db",
           zebrafish="org.Dr.eg.db")

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

runstr <- paste("Call with: Rscript goEnrichment.R rData=<DESeq2.RData> [log2Fold=0]",
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

    calculateGoEnrichment <-  function(de.genes, univ.genes, de.genes.lfc, suffix) {
        # convert to entrezID downregulated/univers genes
        guessKeyType <- function(genes) {
            orgdb <- eval(parse(text=orgDb[org]))
            keytypes(orgdb)[which.max(sapply(keytypes(orgdb), function(x) sum(genes %in% keys(orgdb, keytype=x))))]
        }
        getEntrezId <- function(genes, keytype) {
            bitr( genes, fromType=keytype, toType="ENTREZID", OrgDb=orgDb[org])
        }
        keytype      <- guessKeyType(de.genes)
        entrezDeId   <- getEntrezId(de.genes, keytype)
        entrezUnivId <- getEntrezId(univ.genes, keytype)

        # merge entrezID and log2 foldchange
        entrezDeIdLfc <- merge(x = entrezDeId, y = de.genes.lfc[ , c(type, "log2FoldChange")], by.x=keytype, by.y= type)
        entrezDeIdLfc <- entrezDeIdLfc[!duplicated(entrezDeIdLfc[,1]),]
        rownames(entrezDeIdLfc) <- make.names(entrezDeIdLfc$ENTREZID, unique=TRUE)
        entrezDeIdLfc[,1] <- NULL
        entrezDeIdLfc$ENTREZID <- NULL
        entrezDeIdLfc <- apply(t(entrezDeIdLfc), 2,as.numeric)

        
        enriched         <- enrichGO(entrezDeId$ENTREZID, OrgDb=orgDb[org], keyType="ENTREZID", ont="BP",
                                     readable="SYMBOL" %in% keytypes(eval(parse(text=orgDb[org]))),
                                     universe=if(univ == "all") orgDb[org] else entrezUnivId$ENTREZID)
        enrichedKEGG     <- enrichKEGG(entrezDeId$ENTREZID, keyType="ncbi-geneid", organism=org, universe=entrezUnivId$ENTREZID)
        enrichedReactome <- enrichPathway(entrezDeId$ENTREZID, organism=org,
                                          readable="SYMBOL" %in% keytypes(eval(parse(text=orgDb[org]))),
                                          universe=entrezUnivId$ENTREZID)
             
        # filter enriched results by gene count
        enriched         <- gsfilter(enriched, by = "Count", min = 2)
        enrichedKEGG     <- gsfilter(enrichedKEGG, by = "Count", min = 2)
        enrichedReactome <- gsfilter(enrichedReactome, by = "Count", min = 2)
                                   
        # write GO and Pathway enrichment tables into output file 
        write.csv(as.data.frame(enriched), row.names=FALSE,
                  file=paste0(out, "/", contrast, "_GO_Enrichment_", suffix, "_genes.csv"))
        write.csv(as.data.frame(enrichedKEGG), row.names=FALSE,
                  file=paste0(out, "/", contrast, "_KEGG_Pathway_Enrichment_", suffix, "_genes.csv"))
        write.csv(as.data.frame(enrichedReactome), row.names=FALSE,
                  file=paste0(out, "/", contrast, "_Reactome_Pathway_Enrichment_", suffix, "_genes.csv"))


        if(!is.null(enriched) && nrow(enriched) > 0) {
            # calculate reduced terms and generate reduced term plots and table
            # reduction is only possible, if more than one term enriched
            if (nrow(enriched) > 1) {
                go_analysis  <- as.data.frame(enriched)
                scores       <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
                simMatrix    <- calculateSimMatrix(go_analysis$ID, orgdb=orgDb[org], ont="BP", method="Rel")
                reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb=orgDb[org])
                # this scatter plot can only be created when nrow(enriched)>2, otherwise the scatterPlot function crashes
                if (nrow(enriched) > 2) {
                    CairoPNG(file=paste0(out, "/", contrast, "_GO_scatterplot_", suffix, "_genes.png"), width=1200, height=1200)
                    print(scatterPlot(simMatrix, reducedTerms))
                    dev.off()
                }
                CairoPNG(file=paste0(out, "/", contrast, "_GO_treemap_", suffix, "_genes.png"), width=1200, height=1200)
                treemapPlot(reducedTerms)
                dev.off()

                # overwrite GO enrichment table with parent term information
                x <- merge(go_analysis, reducedTerms, by=1, all.x=TRUE)  # first column is term ID
                write.csv(x, row.names=FALSE, file=paste0(out, "/", contrast, "_GO_Enrichment_", suffix, "_genes.csv"))

                # remove redundant terms from enrichGO enrichment results, and do barplot and cnetplot of top terms
                enriched_reduced <- enriched
                enriched_reduced@result <- subset(enriched_reduced@result, ID %in% unique(reducedTerms$parent))
            } else {
                enriched_reduced <- enriched
            }

            # create barplot showing GO category
            CairoPNG(file=paste0(out, "/", contrast, "_GO_Barplot_", suffix, "_genes.png"), width=1000)
            print(barplot(enriched_reduced, showCategory=plotCategory) + ylab("Number of genes"))
            dev.off()
      
            # create network plot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_GO_network_", suffix, "_genes.png"), width=1000, height=1000)
            print(emapplot(enriched))
            dev.off()

            # create cnetplot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_GO_cnetplot_", suffix, "_genes.png"), width=1200, height=1200)
            print(cnetplot(enriched_reduced, showCategory=plotCategory, categorySize="pvalue", foldChange=entrezDeIdLfc))
            dev.off()
        }

        if(!is.null(enrichedKEGG) && nrow(enrichedKEGG)> 0) {
            # create barplot showing Pathway terms
            CairoPNG(file=paste0(out, "/", contrast, "_KEGG_Barplot_", suffix, "_genes.png"), width=1000)
            print(barplot(enrichedKEGG, showCategory=plotCategory) + ylab("Number of genes"))
            dev.off()
      
            # create network plot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_KEGG_network_", suffix, "_genes.png"), width=1000, height=1000)
            print(emapplot(enrichedKEGG))
            dev.off()
        }
        
        if(!is.null(enrichedReactome) && nrow(enrichedReactome) > 0) {
            # create barplot showing Pathway terms
            CairoPNG(file=paste0(out, "/", contrast, "_Reactome_Barplot_", suffix, "_genes.png"), width=1000)
            print(barplot(enrichedReactome, showCategory=plotCategory) + ylab("Number of genes"))
            dev.off()
            # create network plot for the results
            CairoPNG(file=paste0(out, "/", contrast, "_Reactome_network_", suffix, "_genes.png"), width=1000, height=1000)
            print(emapplot(enrichedReactome))
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
    up   <- resultData[, colpadj] < padj & (resultData[, colfc]) > abs(log2Fold)
    down <- resultData[, colpadj] < padj & (resultData[, colfc]) < -abs(log2Fold)

    # GO and Pathway enrichment analysis for the upregulated genes
    if(sum(up) == 0)
        writeLines("No upregulated genes found!", paste0(out, "/", contrast, "_GO_Enrichment_upregulated_genes.csv"))
    else
        calculateGoEnrichment(resultData[up, genes], resultData[, genes], resultData[up, cbind(genes,colfc)], "up")

    # GO and Pathway enrichment analysis for the downregulated genes
    if(sum(down) == 0)
        writeLines("No downregulated genes found!", paste0(out, "/", contrast, "_GO_Enrichment_downregulated_genes.csv"))
    else
        calculateGoEnrichment(resultData[down, genes], resultData[, genes], resultData[down, cbind(genes,colfc)], "down")

    invisible(0)
}

if(cores > 1 && length(res) > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, c("log2Fold", "padj", "orgDb", "org", "univ", "type", "plotCategory", "out"))
    clusterEvalQ(cl, {
        library(clusterProfiler)
        library(DOSE)
        library(ReactomePA)
        library(Cairo)
        library(ggplot2)
        library(rrvgo)
        library(orgDb[org], character.only=TRUE)
    })
    parLapply(cl, zipup(res, names(res)), processContrast)
    stopCluster(cl)
} else {
    library(clusterProfiler)
    library(DOSE)
    library(ReactomePA)
    library(Cairo)
    library(ggplot2)
    library(rrvgo)
    library(orgDb[org], character.only=TRUE)
    lapply(zipup(res, names(res)), processContrast)
}

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GO_Enrichment_session_info.txt", sep=""))
