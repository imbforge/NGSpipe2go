##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("knitr")        # for markdown output
library("edgeR")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("reshape2")
library("ggrepel")
library("VennDiagram")
library("grid")
library("plotly")
library("GeneOverlap")
library("pheatmap")
library("viridis")
library("gridExtra")
library("dplyr")
library("tidyr")
library("forcats")
library("ngsReports")
library("ggbeeswarm")
library("maser")
library("GenomeInfoDb")
library("kableExtra")

#' loadGlobalVars: read configuration from bpipe vars
#'
#' @param f - a file defining multiple variables for reporting to run. 
#'
#' @return A set of variables that are mentioned in input file, e.g.
#'         SHINYREPS_PROJECT <- "projects/example_project"
#'         SHINYREPS_ORG <- "human"
#'         SHINYREPS_DB <- "hg38"
#'              
#' @description File content should be:
#'              SHINYREPS_PROJECT <- "projects/example_project"
#'              SHINYREPS_ORG=human
#'              SHINYREPS_DB=hg38
#'              
loadGlobalVars <- function(f="shinyReports.txt") {

    # read in the conf file
    conf <- readLines(f)
    conf <- conf[grep("^SHINYREPS_", conf)]
    
    # create the vars
    sapply(conf, function(x) {
        x <- unlist(strsplit(x, "=", fixed=T))
        assign(x[1], x[2], envir=.GlobalEnv)
    })
    
    invisible(0)
}

##
## Some generic functions
##

#'
#' shorten: if a text string is longer than certain length, shorten by showing the first and last characters
#'
#' @param x - a character vector
#' @param max.len - a numeric vector of length 1 indicating the maximal length of a string to be tolerated
#' @param ini - a numeric vector of length 1 indicating how many characters at the beginning of the string should be conserved 
#' @param end - a numeric vector of length 1 indicating how many characters at the end of the string should be conserved 
#'
#' @return a character vector
#'
#' @examples long_string <- "This_is_a_very_long_string_and_should_be_shortened"
#'           short_string <- shorten(long_string, max.len = 20, ini = 5, end = 9)
#'           print(short_string)
#'           ## [1] "This_..._shortened"
#'           
shorten <- function(x, max.len=40, ini=20, end=15) {
    l <- nchar(x)
    if(l > max.len) paste(substr(x, 1, ini), substr(x, (l-end), l), sep="...") else x
}

#'
#' Define a group color palette.
#'
#' @param num - number of groups
#'
#' @return a color vector
#'
#' @examples group.colors <- define.group.palette(length(levels(colData(dds)[,"group"])))
#'           
define.group.palette <- function(num) {

    if (num <=9) {
        group.palette <- brewer.pal(9,"Set1")[1:num]
    } else {
        group.palette <- colorRampPalette(brewer.pal(9,"Set1"))(num)
    }
    return(group.palette)
}

#'
#' Define a replicate color palette.
#'
#' @param num - number of replicates
#'
#' @return a color vector
#'
#' @examples replicate.colors <- define.replicate.palette(length(levels(colData(dds)[,"replicate"])))
#'           
define.replicate.palette <- function(num) {

    if (num <=9) {
        replicate.palette <- brewer.pal(8,"Dark2")[1:num]
    } else {
        replicate.palette <- colorRampPalette(brewer.pal(8,"Dark2"))(num)
    }
    return(replicate.palette)
}

##
## DESeq2 DE analysis
##

#' Create MDS plot from DESeq2 object. 
#'
#' @param dds - a DESeq2 analysis S4 object
#' @param rld - a rlog transformed DESeq2 object
#'
#' @return MDS plot
#'
#' @description This function will use a DESeq2 differential analysis object to create an MDS plot,
#'              which is labelled non-overlapping gene names.
#' 
#' @return A printed plot (ggplot2) object.
#' 
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           rld <- rlog(dds)
#'           # rld and dds are taken from environment
#'           DEhelper.DESeq2.MDS()
#' 
DEhelper.DESeq2.MDS <- function() {

	pca.data <- plotPCA(rld, intgroup=colnames(colData(dds))[1], returnData=TRUE)
	percentVar <- round(100 * attr(pca.data, "percentVar"))
	ggplot(pca.data, aes(PC1, PC2, color=group)) +
 		geom_point(size=2,alpha=0.7) +
 		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		coord_fixed() + 
  		scale_color_manual(values=define.group.palette(length(levels(colData(dds)[,"group"])))) +
  		geom_text_repel(aes(label=colData(dds)[,"replicate"]), show.legend=FALSE) + 
  		theme_bw()

}


#' Create (pairwise) MDS plot from DESeq2 object.
#'
#' @param pairwise.dds - a list of DESeq2 analysis S4 objects
#' @param i - index of the list element to be plotted
#'
#' @return MDS plot
#'
#' @description This function will use a DESeq2 differential analysis object to create an MDS plot,
#'              which is labelled non-overlapping gene names.
#'
#' @return A printed plot (ggplot2) object.
#'
#' @examples Taken from DE_DESeq2.R
#'           pairwise.dds <- lapply(conts[,1],function(cont) {
#'                       ...
#'           })
#'           # pairwise.dds are taken from environment
#'           DEhelper.DESeq2.pairwisePCA(i)
#'
DEhelper.DESeq2.pairwisePCA <- function(i=1) {

	pca.data <- plotPCA(rlog(pairwise.dds[[i]]), intgroup=colnames(colData(pairwise.dds[[i]]))[1], returnData=TRUE)
     	percentVar <- round(100 * attr(pca.data, "percentVar"))
     	print(ggplot(pca.data, aes(PC1, PC2, color=group)) +
		   geom_point(size=2,alpha=0.7) +
		   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
		   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		   coord_fixed() + 
	 	   scale_color_manual(values=brewer.pal(9,"Set1")[1:2]) + 
		   geom_text_repel(aes(label=colData(pairwise.dds[[i]])[,"replicate"]), show.legend=FALSE) +
	 	   theme_bw())
}


## DEhelper.DESeq2.heatmap 
#' Heatmap of sample to sample distances, of top variant 'n' genes of the counts-per-million table or 
#' of top 'n' highest mean genes of the regularized log transformed counts-per-million table 
#' with or without rlog normalization.
#' This function replaces the helper functions DEhelper.DESeq2.corr, DEhelper.DESeq2.corr.pairwise,
#' DEhelper.DESeq2.cluster.sd, DEhelper.DESeq2.cluster.sd.pairwise, 
#' DEhelper.DESeq2.cluster.mean and DEhelper.DESeq2.cluster.mean.pairwise.
#'
#' @param i numeric index, only needed if dds is a list of objects to select a list element.
#' @param dds - rlog transformed DESeq2 object.
#' @param logTransform logical, if TRUE, dds will be transformed using DESeq2::rlog().
#' @param anno_factors character vector with factors given in dds to be used for heatmap column annotation.
#' @param type character with type of heatmap to be drawn. One of "distance", "cluster_sd" or "cluster_mean".
#' @param n amount of transcripts/rows to be plotted (default = 40). Applies only for type "cluster_sd" and "cluster_mean".
#' 
#' @return A plot (base plotting) object.
#'
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           rld <- rlog(dds)
#'           # rld is taken from environment
#'           DEhelper.DESeq2.heatmap()
#'           
DEhelper.DESeq2.heatmap <- function(i=NULL, dds=rld, logTransform=F, anno_factors = c("group","replicate"), type="distance", n=40) {
  
  # for pairwise object select ith element
  if(!is.null(i)) {dds <- dds[[i]]}
  
  # extract assay, optionally logtransform (as for the pairwise.dds objects) and over-write gene_id in rownames with gene_name
  if(logTransform) {
    if(class(dds) != "DESeqTransform") {
      dds <- DESeq2::rlog(dds)
    } else {
      warning("The dds object is of class DESeqTransform, which means it is already transformed. Therefore, the log transformation selected here is omitted.")
    }
  } 
  assay.rld <- assay(dds)
  
  gene_ids <- gtf$gene_name[match(rownames(assay(dds)), gtf$gene_id)]
  rownames(assay.rld) <- ifelse(is.na(gene_ids), rownames(assay(dds)), gene_ids)
  
  
  ## set parameter specific for heatmap type
  switch(type,
         distance={
           # extract sample to sample distances
           sampleDists <- dist(t(assay.rld))
           sampleDistMatrix <- as.matrix(sampleDists)
           # sample to sample distance heatmap w/ 0's on diagonal set to NA
           # since they otherwise shift the scale too much
           sampleDistMatrix.diagNA <- sampleDistMatrix
           sampleDistMatrix.diagNA[sampleDistMatrix.diagNA==0] <- NA
           data2plot <- sampleDistMatrix.diagNA
           # set color scheme
           col <- plasma(255)
           hmmain = NA
           hm_fontsize_row = 8
           hm_fontsize_col = 8
           clustering_distance_rows = sampleDists
           clustering_distance_cols = sampleDists
           cluster_rows=TRUE
         },
         cluster_sd={
           # pick top n most variable genes
           select.highestSD <- order(apply(assay.rld,1,sd),decreasing=TRUE)[1:n]
           data2plot <- assay.rld[select.highestSD,]
           # set color scheme
           col <- colorRampPalette(brewer.pal(9,"GnBu"))(255)
           hmmain=paste("Normalized expression values of",n,"most variable genes") # heatmap title
           hm_fontsize_row = 6
           hm_fontsize_col = 8
           clustering_distance_rows = "euclidean"
           clustering_distance_cols = "euclidean"
           cluster_rows=TRUE
         },
         cluster_mean={
           # pick top n most variable genes
           select.highestMean <- order(rowMeans(assay.rld),decreasing=TRUE)[1:n]
           data2plot <- assay.rld[select.highestMean,]
           # set color scheme
           col <- rev(heat.colors(255))
           hmmain=paste("Normalized expression values of",n,"genes with highest mean") # heatmap title
           hm_fontsize_row = 6
           hm_fontsize_col = 8
           clustering_distance_rows = "euclidean"
           clustering_distance_cols = "euclidean"
           cluster_rows=FALSE
         }
  )
  
  ## extract information for legend if available
  anno_factors <- anno_factors[anno_factors %in% colnames(colData(dds))]
  if(length(anno_factors)>0) {
    
    legend.df <- data.frame(colData(dds)[,c(anno_factors),drop=F],row.names=rownames(colData(dds)))
    
    # fix group colors for legend additional factors if available
    legend_colors <- list()
    
    for (f in anno_factors) {
      if (f == "group") {
        mypalette <- define.group.palette(length(levels(colData(dds)[,f])))
        legend_colors[[f]] <- setNames(mypalette, levels(colData(dds)[,f]))
      } else {
        mypalette <- define.replicate.palette(length(unique(colData(dds)[,f])))
        legend_colors[[f]] <- setNames(mypalette, sort(unique(colData(dds)[,f])))
      }
    }
  } else {
    legend.df <- NA
    legend_colors <- NA
  }  
  
  
  # plot heatmap
  pheatmap(data2plot,
           cluster_rows=cluster_rows,
           cluster_cols=TRUE,
           clustering_distance_rows=clustering_distance_rows,
           clustering_distance_cols=clustering_distance_cols,
           annotation_col=legend.df,
           annotation_colors=legend_colors,
           col=col,
           border_color=NA,
           main=hmmain, 
           fontsize_row=hm_fontsize_row,
           fontsize_col=hm_fontsize_col,
           treeheight_row=20,
           treeheight_col=20,
           annotation_names_col=T)
}





## DEhelper.MAplot: MA plots
#' MA plots of DESeq2 results
#'
#' @param i - iterator to select contrast from conts [default = 1]
#' @param fdr - FDR cutoff [default = 0.01]
#' @param conts - list of contrasts
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#' @return MA plot with genes padj < fdr hichlighted
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           conts <- c("shMed12vsshNMC=(shMed12-shNMC)")
#'           # res & conts are taken from environment
#'           DEhelper.DESeq2.MAplot()
#'           
DEhelper.DESeq2.MAplot <- function(i=1, fdr=.01) {
	cont.name <- gsub("(.+)=(.+)","\\1",conts[i,1])
	plotMA(res[[i]], main=cont.name, alpha=fdr)
}

## DEhelper.DEgenes: show the DE results
#' Show the DE results ordered by adjusted pValue.
#'
#' @param i - iterator to select precalculated DESeq2 result of contrast from result list
#'            [default = 1]
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#' @return data frame ordered by adjusted pvalue
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           # res is taken from environment
#'           DEhelper.DESeq2.DEgenes(i=2)
#'           
DEhelper.DESeq2.DEgenes <- function(i=1) {
    ord  <- order(-log(res[[i]]$padj), 
                   abs(res[[i]]$log2FoldChange), 
                  decreasing=TRUE)
    res[[i]][ord, ]
}

## DEhelper.DESeq2.VolcanoPlot: Volcano plots from DEseq2 results
#' Produce volcano plots from DEseq2 results.
#'
#' @param i - integer, iterator to determine which contrast's DESeq2 result to plot
#' @param fdr - float, FDR cut off of genes to be highlighted
#' @param top - integer, count of genes to be labelled in Volcano plot
#' @param web - boolean
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#'
#' @return ggplotly object, if web = TRUE
#'         printed ggplot2 object, if web = FALSE
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           # res is taken from environment
#'           DEhelper.DESeq2.VolcanoPlot(i = 1, fdr = 0.01, top = 20, web = FALSE)
#'           
DEhelper.DESeq2.VolcanoPlot <- function(i=1, fdr=.01, top=25, web=TRUE) {
    # gather results
    d <- as.data.frame(DEhelper.DESeq2.DEgenes(i))
    x.limit <- max(abs(d$log2FoldChange), na.rm=T) # find the maximum spread of the x-axis
    
    # plotting
    p <- ggplot(d) +
            geom_point(mapping=aes(log2FoldChange, -log10(padj), size=log10(baseMean+1), color=padj < fdr), alpha=0.1) +
            theme_bw() +
            xlim(-x.limit, x.limit) +
            ylab("-log10 adj. p-value") +
            xlab("log2 fold change") + 
            scale_color_manual(values=c("black", "red"), guide=FALSE) +
            scale_size_continuous("mean norm. counts (log10)") +
	    guides(size = guide_legend(nrow=1)) +
       	    theme(legend.position = "top")

    # add name of top genes
    if(top > 0) p <- p + geom_text_repel(data=d[1:min(top, nrow(d)),],
                                         mapping=aes(log2FoldChange, -log10(padj), label=gene_name),
                                         color="black")

    # return plot
    if(web)
        ggplotly(p)
    else
        print(p)
}

## DEhelper.DESeq2.ChrOverrepresentation: test if there are more genes in a chromosome than expected by chance.
#' Test if there are more genes in a chromosome than expected by chance.
#'
#' @param i - integer, iterator to determine which contrast's DESeq2 result to select
#' @param fdr_de_gene - FDR cut-off to determine differentially expressed genes
#' @param fdr_fisher_test - FDR cut-off to determine significantly overrepresented genes on chromosome
#' @param filter - boolean, apply fdr_fisher_test TRUE/FALSE
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#' @return data frame - with biased chromosome, gene count, pvalue and Jaccard index
#'                    - or empty if no bias is found
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           # res is taken from environment
#'           table_fisher <- DEhelper.DESeq2.ChrOverrepresentation(i = 1, fdr_de_gene = 0.1, fdr_fisher_test = 0.1, filter = TRUE)
#'           
DEhelper.DESeq2.ChrOverrepresentation <- function(i=1, fdr_de_gene=0.1, fdr_fisher_test=0.1, filter=TRUE) {
    genes <- data.frame(res[[i]])
    chromosomes <- unique(genes$chr)
    universe <- length(unique(rownames(genes)))

    l_res <- list()
    for (chromosome in chromosomes){
        de_genes <- rownames(subset(genes, padj < fdr_de_gene))
        chr_genes <- rownames(subset(genes, chr == chromosome))
        overl <- newGeneOverlap(
           unique(de_genes),
           unique(chr_genes),
           genome.size=universe)

        overl <- testGeneOverlap(overl)
        l_res[[chromosome]] <- data.frame(
            Chromosome=chromosome,
            TotalGenes=length(chr_genes),
            DEGenes=length(overl@intersection),
            PValue=overl@pval,
            OddsRatio=overl@odds.ratio,
            JaccardIndex=overl@Jaccard
            )
    }

    table_fisher <- do.call("rbind", l_res)

    # correct pvalue for multiple testing
    table_fisher$FDR <- p.adjust(table_fisher$PValue,method="fdr")

    if(filter) table_fisher <- table_fisher[table_fisher$FDR < fdr_fisher_test, ]

    return(table_fisher[c("Chromosome","TotalGenes","DEGenes","PValue","FDR","OddsRatio","JaccardIndex")])
}


##
## DEhelper.STARparms: parse STAR log files and extract all mapping parameters
##
DEhelper.STARparms <- function() {
    
    # log file
    LOG <- SHINYREPS_STAR_LOG
    SUFFIX <- paste0(SHINYREPS_STARparms_SUFFIX, '$')
    if(!file.exists(LOG)) {
        return("STAR statistics not available")
    }
    
    # look for the lines containing the strings and get the values associated with this strings
    parseLog <- function(f) {
        # read in the lines
        f <- file(paste0(LOG, "/", f))
        l <- readLines(f)
        close(f)
        
        # get the version number from the first line (STAR svn revision compiled=STAR_2.3.1z13_r470)
        v <- unlist(strsplit(l[1], "="))[2]
        
        # get the redifined parameters and parse them in a key-value data.frame
        redefined <- l[grep("\\s+\\~RE-DEFINED$", l)]
        redefined <- sapply(redefined, function(x) {
            x <- unlist(strsplit(gsub("\\s+\\~RE-DEFINED$", "", x), "\\s+"))
            x[2] <- if(length(x) < 2) "" else paste(x[-1],collapse=" ")
            x[2] <- shorten(x[2])
            x[c(1, 2)]
        })

        # take the 1st time a parm appears re-defined (skip re-defined in index generation)
        redefined <- redefined[, !duplicated(redefined[1, ])] 

        # create a vector with the values
        x <- redefined[2, ]
        names(x) <- redefined[1, ]
        
        # put the STAR version number
        x <- c(version=v, x)
        return(x)
    }
    df <- sapply(list.files(LOG, pattern=SUFFIX), parseLog)
    
    # remove variable lines (lines depending on the fastq.gz file name)
    # and check if all the columns contain the same value. Display a warning otherwise
    df <- df[!grepl("(outFileNamePrefix|outTmpDir|readFilesIn)", rownames(df)), ]
    l <- apply(df, 1, function(x) length(unique(x)))    # rows differing (l > 1)
    df <- as.data.frame(df[, 1, drop=F])    # keep only the first column
    colnames(df) <- "parms"
    df$warning[l > 1] <- "Some files aligned with a different parm. Check logs"
    
    # set row and column names, and output the md table
    if(all(is.na(df$warning))) {
        kable(df[, 1, drop=F], align=c("r"), output=F, format="markdown")
    } else {
        kable(df, align=c("r", "r"), output=F, format="markdown")
    }
}

##
## DEhelper.STAR: parse STAR log files and create a md table
##
DEhelper.STAR <- function() {
    
    # log file
    LOG <- SHINYREPS_STAR_LOG
    SUFFIX <- paste0(SHINYREPS_STAR_SUFFIX, '$')
    if(!file.exists(LOG)) {
        return("STAR statistics not available")
    }
    
    # look for the lines containing the strings
    # and get the values associated with this strings
    x <- sapply(list.files(LOG, pattern = SUFFIX), function(f) {
        f <- file(paste0(LOG, "/", f))
        l <- readLines(f)
        close(f)
        
        sapply(c("Number of input reads",
                 "Uniquely mapped reads number",
                 "Uniquely mapped reads %",
                 "Number of reads mapped to multiple loci",
		 "% of reads mapped to multiple loci",
                 "Number of reads mapped to too many loci",
                 "% of reads mapped to too many loci",
                 "% of reads unmapped: too many mismatches",
                 "% of reads unmapped: too short",
                 "% of reads unmapped: other"), function(x) {
                     as.numeric(gsub("%", "", gsub(".+\\|\t(.+)", "\\1", l[grep(x, l)])))
                 })    
    })
    
    # set row and column names, and output the md table
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))

    ## row.names need to be set in case of only one sample (if more than one sample,
    ## rownames would be set automatically, but doesnt harm to set them)
    df_values <- as.data.frame(t(x[1:7,]),row.names=colnames(x))

    df_values["unmapped"] <- x[1, ] - x[2, ] - x[4, ] - x[6,]
    df_values["% unmapped"] <- x[8, ] + x[9, ] + x[10, ]
    df_values$sample <- rownames(df_values)
    # we clean up the colnames a little to make them shorter and nicer
    colnames(df_values) <- gsub("of reads mapped to ", "",
				gsub("Number of reads mapped to ", "",
				     gsub("Uniquely mapped reads %","% unique",
					  gsub("Uniquely mapped reads number","unique",
					       gsub("Number of input reads","input",colnames(df_values))))))

    # if we have a differential expression analysis
    # we refactor the samples depending on group/subject or alternatively on the
    # amount of unique_mapping reads
    if(file.exists(SHINYREPS_TARGET)){

        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "",targets$file )
        add_factors <- colnames(targets)[!colnames(targets) %in% c("group", "sample", "file")]

        # replace files names with nicer sample names given in targets file 
        # if sample is missing in targets file, use reduced file name
        df_values$sample <- sapply(df_values$sample, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,   
                                                                          targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                                          gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else{

      # remove sample prefix from sample names (suffix was already removed before)
      df_values$sample <- gsub(paste0("^",SHINYREPS_PREFIX), "", df_values$sample)
    }

    # sort sample alphabetically or based on number of (uniquely) mapped reads
    if (SHINYREPS_SORT_ALPHA) {
        df_values$sample <- fct_reorder(df_values$sample, df_values$sample, .desc=TRUE)
    } else {
        df_values$sample <- fct_reorder(df_values$sample, df_values$`% unique`)
    }

    df_melt <- reshape2::melt(df_values, value.name = "reads", variable.name = "mapping_stat")
    df_melt$value_info <- ifelse(grepl("%", df_melt$mapping_stat), "perc", "reads")

    ## melt does not work properly in case of only one sample (sample name gets lost)
    if (nrow(df_values) == 1) {
	    df_melt$sample <- rownames(df_values)
    }

    ## invert order for plotting
    df_melt$mapping_stat <- fct_rev(df_melt$mapping_stat)
    
    ## define color scheme
    my.palette <- rev(c("#e08a28","#dfc27d","#acd2f7","#5ea1e3"))
    my.color   <- "#cc0e25"

    # plot showing percent mapped reads
    p_perc <- ggplot(df_melt[df_melt$value_info == "perc",], 
		     aes(x = sample, y = reads, fill = mapping_stat )) +
	    geom_bar(stat     = "identity", 
		     position = "stack") +
	    labs(x = "", 
		 y = "% of sequenced reads", 
		 title = "Mapping summary (percentage)") + 
	    theme(axis.text.y = element_text(size = 8),
		  plot.title = element_text(hjust = 0.5,size=12),
		  legend.position = "top") +
	    guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) +
	    scale_fill_manual(values=my.palette) +
	    coord_flip()

    # max number of reads in any sample
    max.reads <- max(df_melt[df_melt$mapping_stat=="input","reads"])

    # plot a combination of percent mapped reads and total sequenced reads
    p_perc_count <- ggplot() + 
      geom_bar(data=df_melt[df_melt$value_info == "reads" & df_melt$mapping_stat!="input",], 
               mapping=aes(x = sample, y = reads, fill = mapping_stat), 
               stat = "identity", position = "fill", width = 0.8) + 
      geom_point(data=df_melt[df_melt$mapping_stat=="input",], 
                 mapping=aes(x = sample, y = reads/max.reads), 
                 size=2, fill=my.color, color=my.color, shape=18) +
      geom_line(data=df_melt[df_melt$mapping_stat=="input",],
                mapping=aes(x = sample, y = reads/max.reads, group=1),
                color=my.color, linetype="dashed", size=0.2) +
      scale_y_continuous(sec.axis = sec_axis(~ . *max.reads, name="# sequenced reads"),
                         labels = scales::percent_format()) + 
      labs(x = "",
           y = "% of sequenced reads") + 
      theme(axis.text.y = element_text(size = 8),
            axis.title.x.top=element_text(color=my.color), 
            axis.text.x.top=element_text(color=my.color), 
            axis.ticks.x.top=element_line(color=my.color),
            plot.title = element_text(hjust = 0.5,size=12),
            legend.position = "top") + 
      guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) + 
      scale_fill_manual(values=my.palette) + 
      coord_flip()

    # in case of non-alphabetically sorting, re-order based on total number of reads
    # before plotting reads counts
    if (SHINYREPS_SORT_ALPHA == FALSE) {
        df_melt$sample <- factor(df_melt$sample, levels=df_values$sample[order(df_values$input)])
    }

    p_count <- p_perc %+%
               df_melt[df_melt$value_info == "reads" & df_melt$mapping_stat != "input",] +
	       labs(x = "",
		    y = "# sequenced reads",
		    title = "Mapping summary (read counts)") 
    
    # prepare data for summary table  
    rownames(df_values) <- df_values$sample
    df_values <- df_values[, colnames(df_values) != "sample"]

    # reformat individual columns
    df_values[, grepl("%", colnames(df_values))] <- as.data.frame(
      lapply(
        df_values[, grepl("%", colnames(df_values))], function(x){
          paste0(format(x, nsmall=2), "%") 
        }))
    df_values[, !grepl("%", colnames(df_values))] <- as.data.frame(
      lapply(
        df_values[, !grepl("%", colnames(df_values))], function(x){
           format(x, big.mark=",") 
        }))

    # create data frame to print as table
    df_values_print <- data.frame(all.reads     = df_values$input,
				  unique        = paste0(df_values[,"unique"]," (",df_values[,"% unique"],")"),
				  multiple.loci = paste0(df_values[,"multiple loci"]," (",df_values[,"% multiple loci"],")"),
				  too.many.loci = paste0(df_values[,"too many loci"]," (",df_values[,"% too many loci"],")"),
				  unmapped      = paste0(df_values[,"unmapped"]," (",df_values[,"% unmapped"],")"), 
				  row.names     = rownames(df_values))

    # sort alphabetically based on sample names (=rownames) or based on % unique
    if (SHINYREPS_SORT_ALPHA == TRUE) {
        df_values_print <- df_values_print[order(rownames(df_values_print)),]
    } else {
        df_values_print <- df_values_print[rownames(df_values[order(df_values$`% unique`, decreasing=TRUE),]),]
    }
    colnames(df_values_print) <- gsub("\\."," ",colnames(df_values_print))
                     
    return( list(p_perc = p_perc,
		 p_perc_count = p_perc_count,
                 p_count = p_count,
		 stat = kable(df_values_print, align=c("r", "r", "r", "r", "r"), format="markdown", output=F))
    )
    
}

##
## DEhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
DEhelper.Fastqc <- function(web=FALSE, subdir="") {
    
  # logs folder
  if(!all(sapply(file.path(SHINYREPS_FASTQC_OUT, subdir), file.exists))) {
    return(paste("Fastqc statistics not available for", names(which(!sapply(file.path(SHINYREPS_FASTQC_OUT, subdir), file.exists)))))
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.dirs(QC, recursive=F, full.names = T)
  samples <- samples[sapply(samples, function(x) {file.exists(file.path(x, "fastqc_data.txt"))})] # exclude potential subdir which is also listed by list.dirs
  
  df <- sapply(samples, function(f) {
        c(paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"),
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_sequence_gc_content.png)"))
    })

    # set row and column names, and output the md table
    df <- as.data.frame(t(df))
    rownames(df) <-  basename(samples)
    colnames(df) <- c("Read qualities", "Sequence bias", "GC content")
    
    if(file.exists(SHINYREPS_TARGET)){
      
      # get target names
      targets <- read.delim(SHINYREPS_TARGET)
      targets$sample_ext <- gsub("\\..*$", "",targets$file )
      
      # replace files names with nicer sample names given in targets file
      # if sample is missing in targets file, use reduced file name
      rownames(df) <- sapply(rownames(df), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                                ifelse(sapply("R1", grepl, i), 
                                                                       paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R1"),
                                                                       ifelse(sapply("R2", grepl, i), 
                                                                              paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R2"),
                                                                              targets[sapply(targets$sample_ext, grepl, i),"sample"])),
                                                                gsub(paste0("^",SHINYREPS_PREFIX),"",i))})                                                    
      
    } else {
      if(!is.na(SHINYREPS_PREFIX)) {
        rownames(df) <- gsub(paste0("^",SHINYREPS_PREFIX), "", rownames(df))
      }
      rownames(df) <- gsub("_fastqc$", "", rownames(df))
      rownames(df) <- sapply(rownames(df), shorten)
    }
    
    # add a row with the sample name (as given in the rownames) before every row
    df.new <- do.call(rbind,lapply(1:nrow(df),function(i) {rbind(c("",rownames(df)[i],""),df[i,])}))
    rownames(df.new) <- NULL
    # kable(df.new, output=F, align="c", format="markdown") # print sample names in additional rows
    kable(df, output=F, align="c") # print sample names as rownames
}


##
## DEhelper.ngsReports.Fastqc: joint FastQC report of all samples in the experiment
##
DEhelper.ngsReports.Fastqc <- function(subdir="") {
	
	# output folder
	if(!file.exists(SHINYREPS_FASTQC_OUT)) {
		return("Fastqc statistics not available")
	}

    # Loading FastQC Data 
  f <- list.files(file.path(SHINYREPS_FASTQC_OUT, subdir), pattern="fastqc.zip$", full.names=TRUE)
  x <- ngsReports::FastqcDataList(f)
	lbls <- gsub(paste0("(^", SHINYREPS_PREFIX, "|.fastqc.zip$)"), "", names(x))
    names(lbls) <- gsub(".fastqc.zip", ".fastq.gz", names(x))

    print(ngsReports::plotBaseQuals(x, labels=lbls))
    print(ngsReports::plotSeqContent(x, labels=lbls) +
            theme(legend.position="right") +
            guides(fill=FALSE, color="legend") +
            geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                       data=data.frame(base=c("T", "A", "C", "G")),
                       inherit.aes=FALSE, show.legend=TRUE) +
            scale_color_manual("", values=c("red", "green", "blue", "black"))
    )
    print(ngsReports::plotGcContent(x, plotType="line", gcType="Genome", labels=lbls))  
}

##
## DEhelper.Fastqc.custom: prepare Fastqc summary plots
##
DEhelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, subdir="") {
  
    # logs folder
    if(!file.exists(SHINYREPS_FASTQC_OUT)) {
        return("Fastqc statistics not available")
    }
  
    # construct the folder name, which is different for web and noweb
    QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
    # read fastqc results in the appropriate format
    f <- list.files(QC, pattern="\\.zip$",full.names=T)
    fastqc.stats <- ngsReports::FastqcDataList(f)

    # create proper name vectoir as labels
    lbls <- gsub("_fastqc.zip$", "", names(fastqc.stats))
    names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", names(fastqc.stats))

    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "",targets$file )
        
        # replace files names with nicer sample names given in targets file 
        # if sample is missing in targets file, use reduced file name
        lbls <- sapply(lbls, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                  targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                                  gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
        
        if(SHINYREPS_PAIRED == "yes") {
            x <- names(lbls)
            lbls <- paste0(lbls, ifelse(grepl("R1", names(lbls)), ".R1", ".R2"))
            names(lbls) <- x
        }
    } else {
        
        if(!is.na(SHINYREPS_PREFIX)) {
            lbls <- gsub(paste0("^",SHINYREPS_PREFIX), "", lbls)
        }
    }

    # change names also in fastqc.stats (needed for seq. quality plot)
    names(fastqc.stats) <- lbls

    if (summarizedPlots == TRUE) {

        # prepare for plotting  
        df <- reshape2::melt(lapply(fastqc.stats , function(x) x@Per_base_sequence_quality[, c("Base","Mean")]))
        names(df)[names(df)=="L1"] <- "samplename"
  
        # color code the samples as done by fastqc:
        # A warning will be issued if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
        # A failure will be raised if the lower quartile for any base is less than 5 or if the median for any base is less than 20. 
        # (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)
        cols <- c(pass    = "#5cb85c",
                  warning = "#f0ad4e",
                  fail    = "#d9534f")
        colorcode <- do.call(rbind,lapply(names(fastqc.stats),
                                          function(i) {
                                              min.l.quart <- min(fastqc.stats[[i]]@Per_base_sequence_quality$Lower_Quartile)
                                              min.med <- min(fastqc.stats[[i]]@Per_base_sequence_quality$Median)
                                              col.sample <- ifelse(min.l.quart>=10 & min.med>=25, 
                                                                   cols["pass"],
                                                                   ifelse(min.l.quart>=5 & min.med>=20,
                                                                          cols["warning"],
                                                                          cols["fail"]))
                                              return(data.frame(sample=i,
                                                                min.lower.quart=min.l.quart,
                                                                min.median=min.med,
                                                                col=col.sample))
                                          }
        ))
  
        ## only label "warning"s and "fail"ures 
        to.be.labelled <- colorcode[colorcode$col == cols["warning"] | colorcode$col == cols["fail"],]

        ## in case all samples "pass"ed, label "all" in legend
        if (nrow(to.be.labelled) == 0) {
            to.be.labelled <- data.frame(sample="overlay of all samples", col=cols["pass"], row.names=NULL)
        } else {
            to.be.labelled <- rbind(to.be.labelled,c("all other samples","","",col=cols["pass"]))
        }

        # fix position on x-axis since there can be intervals of positions summarized in fastqc
        df$position <- factor(df$Base, levels=unique(df$Base))
        xlen <- length(unique(df$Base))

        p.qual <- ggplot(df, aes(x=as.numeric(position), y=value)) +
           labs(x = "position in read (bp)",
                y = "mean quality score (Phred)") +
           geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 20), fill = "#edc0c4", alpha = 0.3, color=NA) +
           geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 28), fill = "#f0e2cc", alpha = 0.3, color=NA) +
           geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 28, ymax = Inf), fill = "#ceebd1", alpha = 0.3, color=NA) +
           geom_line(aes(color=samplename)) +
           scale_color_manual(values = as.character(colorcode$col),
                              breaks = colorcode$sample,
                              labels = colorcode$sample) +
           geom_point(data=to.be.labelled,
                      mapping=aes(x=NaN, y=NaN, fill=sample),
                      inherit.aes=FALSE, show.legend=TRUE, size = 1.5, shape = 21, color = "white") +
           scale_fill_manual(values = as.character(to.be.labelled$col),
                             breaks = to.be.labelled$sample,
                             labels = to.be.labelled$sample) +
           geom_hline(yintercept = c(0,10,20,30,40),color="white",alpha=0.3) +
           geom_vline(xintercept = seq(0,xlen,10),color="white",alpha=0.3) +
           coord_cartesian(xlim = c(1,xlen), ylim = c(0,42)) +
           guides(color=FALSE,
                  fill=guide_legend(title="",ncol=3)) +
           theme(axis.text.x = element_text(size=6,angle=90,hjust=0.5,vjust=0.5),
                 axis.text.y = element_text(size=8),
                 axis.title  = element_text(size=10),
                 plot.title  = element_text(size=12),
                 legend.text = element_text(size=7),
                 legend.position = "top") +
           scale_x_continuous(breaks=unique(as.numeric(df$position)),
                              labels=unique(df$Base))

        ## use rev(lbls) to plot in reverse order
        p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=rev(lbls)) +
          labs(y = "") +
          theme(legend.position="right") +
          guides(fill=FALSE, color="legend") +
          geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                     data=data.frame(base=c("T", "A", "C", "G")),
                     inherit.aes=FALSE, show.legend=TRUE) +
          scale_color_manual("", values=c("red", "green", "blue", "black")) 

    } else {

        p.qual <- ngsReports::plotBaseQuals(fastqc.stats, labels=lbls, plotType="boxplot") +
          theme(axis.text.x = element_text(size=5))
        p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=lbls, plotType="line") +
          theme(axis.text.x = element_text(size=5), legend.position = "top")
    }

    # GC content line plot 
    # in case you want to add a theoretical distribution to the plot, use function plotGcContent with 
    # the following settings:
    # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=TRUE, species=SPECIES)
    # the default value for SPECIES is "Hsapiens", thus, if you don't specify it, human will be used as a default
    p.gc <- ngsReports::plotGcContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=FALSE) 
    if(!summarizedPlots) {
      p.gc <- p.gc + guides(color=guide_legend(title="",ncol=4)) + 
        theme(legend.position = "top", legend.text = element_text(size=8)) 
    }
    
    return(list(no.of.samples=length(f), p.qual=p.qual, p.content=p.content, p.gc=p.gc))
}


##
## DEhelper.fastqscreen: add FastqScreen data to and plot it as a barplot
##
#' DEhelper.fastqscreen: summarizes FastQScreen results, creates summarized barplots, only relevant contanimants shown
#'
#' @param perc.to.plot - a numeric vector of length 1 setting the percent cutoff of relevant contaminants, if any sample
#'                       shows more than perc.to.plot, contaminant will be shown in plot
#'
#' @return a list including a plot, the number of samples, and the number of plotted contaminants
#'
DEhelper.fastqscreen <- function(perc.to.plot = 1) {
  
  # logs folder
  if(!file.exists(SHINYREPS_FASTQSCREEN_OUT)) {
    return(list(errortext="FastQScreen statistics not available",
                no.of.genomes=1,
                no.of.samples=1))
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- SHINYREPS_FASTQSCREEN_OUT
   
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(SHINYREPS_FASTQSCREEN_OUT, pattern="_screen.txt$", recursive=T, full.names=T)
  df <- lapply(samples, function(f) {
    #we read in the file 
    screen_data <- read.delim(f, header=T, skip=1)
    #the last line is our %hit_no_genomes
    no_hit <- as.numeric(gsub("%Hit_no_genomes: ", "",screen_data[nrow(screen_data),1]))
    screen_data <- screen_data[-nrow(screen_data),]
    rownames(screen_data) <- screen_data$Genome
    screen_data <- screen_data[, !(colnames(screen_data)=="Genome")]
    #we get the number of unique/multiple hits per one genome and calculate sum (of the percentages)
    one_genome <- data.frame(perc=rowSums(screen_data[, 
                                                      grepl("one_genome.1",
                                                            colnames(screen_data))]))
    one_genome$genome <- rownames(one_genome)
    one_genome$category <- "one genome"
    multi_genome <- data.frame(perc=rowSums(screen_data[, colnames(screen_data) %in% 
                                                          c("X.One_hit_multiple_genomes.1",
                                                            "X.Multiple_hits_multiple_genomes")]))
    multi_genome$genome <- rownames(multi_genome)
    multi_genome$category <- "multiple genomes"
    
    mapping_info <- rbind(one_genome, multi_genome)
    mapping_info <- rbind(mapping_info, c(perc=no_hit, genome="no hit", category="no hit"))
    mapping_info$perc <- as.numeric(mapping_info$perc)
    mapping_info$category <- factor(mapping_info$category, levels=c("one genome","multiple genomes","no hit"))
    return(mapping_info)
  })
  
  # sample names
  samples <- gsub("_screen.txt", "", basename(samples))
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET)
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    samples <- sapply(samples, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                    ifelse(sapply("R1", grepl, i), 
                                                           paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R1"),
                                                           ifelse(sapply("R2", grepl, i), 
                                                                  paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R2"),
                                                                  targets[sapply(targets$sample_ext, grepl, i),"sample"])),
                                                    gsub(paste0("^",SHINYREPS_PREFIX),"",i))})                                                    
  } else {
    if(!is.na(SHINYREPS_PREFIX)) {
      samples <- gsub(paste0("^",SHINYREPS_PREFIX), "", samples)
    }
  }
  
  names(df) <- samples
  #createing a df out of the lists
  df <- reshape2::melt(df, value.name="perc")
  colnames(df) <- gsub("L1", "sample", colnames(df))

  # filter for relevant contaminants (e.g. showing >=1% (perc.to.plot) in any sample)
  max.per.genome <- aggregate(df$perc,list(df$genome),max)
  relevant.genomes <- max.per.genome[max.per.genome[,2]>=perc.to.plot,1]
  df <- df[df$genome %in% relevant.genomes,]

  # sort alphabetically
  df$sample <- factor(df$sample, levels=unique(df$sample)[order(unique(df$sample),decreasing=TRUE)])
  
  # replace "Mycoplasma" by "Mycoplasma species", FastQScreen itself cannot deal with space characters
  df$genome <- gsub("Mycoplasma","Mycoplasma species",df$genome)

  # split/wrap per genome
  df$genome <- factor(df$genome, levels=unique(df$genome))

  p.category.wrap <- ggplot(df, aes(x=sample, y=perc, fill=category)) +
          geom_col(position=position_stack(reverse=T),width=0.8) +
          scale_fill_manual(values=c(alpha("#4281a4",0.8),alpha("#ffa62b",0.8),"gray60")) +   
          scale_y_continuous(breaks=seq(0,100,by=10)) +
          theme_bw(base_size=10) +
          labs(x = "",
               y = "% mapped") +
          theme(axis.text.x = element_text(vjust=0.5, angle=90),
                legend.position = "top") +
          guides(fill=guide_legend(title="mapped to", ncol=3)) +
          facet_wrap(~genome,ncol=2) +
          coord_flip()
  
  return(list(p.category.wrap=p.category.wrap,
              no.of.genomes=length(unique(df$genome)),
              no.of.samples=length(unique(df$sample))))      
}



##
## DEhelper.dupRadar: go through dupRadar output dir and create a md table with
##     the duplication plots
##
DEhelper.dupRadar <- function(web=FALSE) {
    
    # logs folder
    if(!file.exists(SHINYREPS_DUPRADAR_LOG)) {
        return("DupRadar statistics not available")
    }

    SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){3})
    if(SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 3L    # default to 4 columns
    }

    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/dupRadar" else SHINYREPS_DUPRADAR_LOG
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_DUPRADAR_LOG, pattern=".png$")
    df <- sapply(samples, function(f) {
        paste0("![dupRadar img](", QC, "/", basename(f), ")")
    })
    
    # sample names 
    samples <- sapply(df, function(x) {
        gsub("_dupRadar.png)$", "", basename(x))
    })

    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "",targets$file )
        
        # replace files names with nicer sample names given in targets file
        # if sample is missing in targets file, use reduced file name
        samples <- sapply(samples, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,   
                                                        targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                        gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {
        if(!is.na(SHINYREPS_PREFIX)) {
            samples <- gsub(paste0("^",SHINYREPS_PREFIX), "", samples)
        }
    }

    # sort alphabetically
    samples <- samples[order(samples)]
    df <- df[names(samples)]

    # fill up additional columns if number of samples is not a multiple of SHINYREPS_PLOTS_COLUMN
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    while(length(samples) %% SHINYREPS_PLOTS_COLUMN != 0) samples <- c(samples, "")

    df      <- matrix(df     , ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    samples <- matrix(samples, ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    colnames(df.names) <- rep(" ", SHINYREPS_PLOTS_COLUMN)
    
    kable(as.data.frame(df.names), align="c", output=F, format="markdown")
}

##
## DEhelper.RNAtypes: parse Subread count results for RNAtypes usage
##
DEhelper.RNAtypes <- function() {
    
    FOLDER <- SHINYREPS_RNATYPES
    SUFFIX <- paste0(SHINYREPS_RNATYPES_SUFFIX, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return("Subread statistics not available")
    }

    if(file.exists(SHINYREPS_TARGET)){
        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file)
        add_factors <- colnames(targets)[!colnames(targets) %in% c("group", "sample", "file")]
    }

    # create a matrix using feature names as rownames, sample names as colnames
    l.infiles <- list.files(FOLDER, pattern=SUFFIX, full.names = T)
    l.counts <- lapply(l.infiles, function(f) {
       
      	# extract sample name and replace it with nicer sample name given in targets file (if available)
              samplename <- gsub(paste0(SUFFIX,"$"), "", basename(f))
      
      	if(file.exists(SHINYREPS_TARGET)){
                  # replace files names with nicer sample names given in targets file 
                  # if sample is missing in targets file, use reduced file name
                  samplename <- ifelse(sum(sapply(targets$sample_ext, grepl, samplename))==1,
                                       targets[sapply(targets$sample_ext, grepl, samplename),"sample"],
                                       gsub(paste0("^",SHINYREPS_PREFIX),"",samplename))
              } else {
                  if(!is.na(SHINYREPS_PREFIX)) { samplename <- gsub(paste0("^",SHINYREPS_PREFIX), "", samplename) }
              }
              df.readcount <- read.table(f, header=T, sep='\t', comment.char = '#', col.names=c("Geneid","Length",samplename))
              df.readcount$Length <- NULL
              return(df.readcount)
    })
    
    # merge counts to data frame in wide format & remove "length..." columns
    f.merge <- function(x,y) {merge(x,y,by="Geneid")}
    df.counts <- Reduce(f.merge, l.counts)
    rm(l.counts)
    
    # eradicate biotype classes that are not present (or could be summarized as "other" being <1%)
    v.countsums <- rowSums(df.counts[-1])
    v.relsums <- v.countsums / sum(v.countsums)
    df.counts$Geneid[v.relsums < 0.01] <- "other"
    
    df.counts <- aggregate(df.counts[-1],
                           by=list(Geneid =df.counts$Geneid),
                           FUN=sum)
    
    # plot
    df.counts.melt <- melt(df.counts, id.var="Geneid")
    colnames(df.counts.melt) <- c("type","sample","count")
    
    # remove possible starting "X" (in case sample name starts with a number)
    df.counts.melt$sample <- gsub("^X","",df.counts.melt$sample)

    # get reverse order for plotting
    samplenames <- gsub("^X","",colnames(df.counts)[-1])
    df.counts.melt$sample <- factor(df.counts.melt$sample, levels=samplenames[order(samplenames, decreasing=TRUE)])

    plot <- ggplot() + 
        geom_bar(data=df.counts.melt, aes(x=sample, y=count, fill=type), position="fill", stat="identity") + 
        labs(x="", y="", fill="") +
	theme(axis.text.y = element_text(size=8)) +
	guides(fill = guide_legend(reverse=TRUE)) + 
	coord_flip() 
    
    return(plot)
}


##
## DEhelper.geneBodyCov: go through geneBodyCov output dir and create a md table with
##     the gene coverage plots
##
DEhelper.geneBodyCov <- function(web=FALSE) {
    
    # logs folder
    if(!file.exists(SHINYREPS_GENEBODYCOV_LOG)) {
        return("geneBodyCov statistics not available")
    }
    
    SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){3})
    if(SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/geneBodyCov" else SHINYREPS_GENEBODYCOV_LOG
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_GENEBODYCOV_LOG, pattern=".png$")
    df <- sapply(samples, function(f) {
        paste0("![geneBodyCov img](", QC, "/", basename(f), ")")
    })
    
    # sample names
    samples <- sapply(df, function(x) {
        gsub("_geneBodyCov.png)$", "", basename(x))
    })

    if(file.exists(SHINYREPS_TARGET)){
        
        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\.*$", "",targets$file )

        # replace files names with nicer sample names given in targets file
	      # if sample is missing in targets file, use reduced file name
        samples <- sapply(samples, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                        targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                                        gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {
        if(!is.na(SHINYREPS_PREFIX)) {
            samples <- gsub(paste0("^",SHINYREPS_PREFIX), "", samples)
        }
    }

    # sort alphabetically
    samples <- samples[order(samples)]
    df <- df[names(samples)]

    # fill up additional columns if number of samples is not a multiple of SHINYREPS_PLOTS_COLUMN
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    while(length(samples) %% SHINYREPS_PLOTS_COLUMN != 0) samples <- c(samples, "")

    df      <- matrix(df     , ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    samples <- matrix(samples, ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    colnames(df.names) <- rep(" ", SHINYREPS_PLOTS_COLUMN)
    
    kable(as.data.frame(df.names), align="c", output=F, format="markdown")
}

##
## DEhelper.geneBodyCov2: go through geneBodyCov output dir and plot into one plot
##
DEhelper.geneBodyCov2 <- function(web=F) {
    
    # logs folder
    if(!file.exists(SHINYREPS_GENEBODYCOV_LOG)) {
        return("geneBodyCov statistics not available")
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/geneBodyCov" else SHINYREPS_GENEBODYCOV_LOG
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(QC, pattern=".csv$", full.names=T)
    names(samples) <- gsub("_geneBodyCov.csv", "", basename(samples))
    
    if(file.exists(SHINYREPS_TARGET)){
        
        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file )

        # replace files names with nicer sample names given in targets file
	      # if sample is missing in targets file, use reduced file name
        names(samples) <- sapply(names(samples), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,  
                                                                      targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                                      gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {
        if(!is.na(SHINYREPS_PREFIX)) {
            names(samples) <- gsub(paste0("^",SHINYREPS_PREFIX), "", names(samples))
        }
    }
    
    l <- lapply(1:length(samples), function(i) {
             df <- read.csv(samples[[i]]) 
             colnames(df) <- c("perc", names(samples)[i])
             return(df[,2, drop=F])
      })
    df <- do.call(cbind, l)
    df$perc <- 1:nrow(df) #these are the 100 bins used in the original plot
    df <- reshape2::melt(df, id.vars="perc", variable.name="sample", value.name="cov")
    #sort the samples either by the group order of the targets file or alphabetically
    #if we do not have a targets file
    if(file.exists(SHINYREPS_TARGET)){
       sample_order <- targets$sample[order(paste0(targets$group, targets$sample))]
       if(all(sample_order %in% df$sample)){
          df$sample <- factor(df$sample, levels = sample_order)
          df$group <- factor(as.character(targets$group[match(df$sample, targets$sample)]),
                             levels = sort(unique(targets$group)))
          if("replicate" %in% colnames(targets)){
             df$replicate <- factor(as.character(targets$replicate[match(df$sample, targets$sample)]), 
                                    levels = sort(unique(as.character(targets$replicate))))
          }
       }
    }else{
      # sort alphabetically
      df$sample <- factor(df$sample, levels = sort(unique(df$sample)))  
    }
    
    # plot
    plot_df <- highlight_key(df, ~sample)
    p <- ggplot(plot_df, aes(x=perc, y=cov, group=sample)) +
           geom_line(color="darkgrey") +
           labs(title="",
                x = "Gene body percentile 5' -> 3'",
                y = "Averaged normalised coverage") +
           scale_color_hue( c=50, l=70) +
           ylim(0,1) +
           theme_bw()
    gg <- ggplotly(p) 

    # since plotly might not work for everyone we should additionally create some static plots
    num_samples <- length(unique(df$sample))
    if(num_samples < 10){ #we only have 9 colors Set1
      sample.cols <- brewer.pal(num_samples, "Set1")
    }else{
      sample.cols <- colorRampPalette(brewer.pal(9, "Set1"))(num_samples)
    }

    plot_list <- list()
    plot_list[["plotly"]] <- gg

    if(!any(c("group","replicate") %in% colnames(df))){
       p_per_sample <- ggplot(df, aes(x=perc,
                                      y=cov,
                                      group=sample,
                                      color=sample)) +
           geom_line() +
           labs(title="",
                x = "Gene body percentile 5' -> 3'",
                y = "Averaged normalised coverage") +
           scale_color_manual(values=sample.cols) +
           ylim(0,1) +
           theme_bw() +
           theme(legend.title=element_blank()) 
       plot_list[["p_per_sample"]] <- p_per_sample
    }

    if(("group" %in% colnames(df)) && !("replicate" %in% colnames(df))){
      num_groups <- length(unique(df$group))
      plot_list[["num_groups"]] <- num_groups
      p_per_sample_splitByGroup <- ggplot(df, aes(x=perc,
                                                  y=cov,
                                                  group=sample,
                                                  color=sample)) +
        geom_line() +
        labs(title="",
             x = "Gene body percentile 5' -> 3'",
             y = "Averaged normalised coverage") +
        scale_color_manual(values=sample.cols) +
        ylim(0,1) +
        theme_bw() +
        theme(legend.title=element_blank())  + 
        facet_wrap(~group,ncol=2) +
        theme(legend.position="top",
              plot.title=element_blank()) +
        guides(color=guide_legend(ncol=3))
      plot_list[["p_per_sample_splitByGroup"]] <- p_per_sample_splitByGroup
    }

    if(all(c("group","replicate") %in% colnames(df))){
       num_groups <- length(unique(df$group))
       num_replicates <- length(unique(df$replicate))
       plot_list[["num_groups"]] <- num_groups
       plot_list[["num_replicates"]] <- num_replicates
       
       p_per_replicate_splitByGroup <- ggplot(df, aes(x=perc,
                                                      y=cov,
                                                      group=sample,
                                                      color=replicate)) +
           geom_line() +
           labs(x = "Gene body percentile 5' -> 3'",
                y = "Averaged normalised coverage") +
           scale_color_manual(values=define.replicate.palette(num_replicates)) +
           ylim(0,1) +
           theme_bw() + 
           theme(legend.position="top") +
           guides(color=guide_legend(ncol=6)) +
	   facet_wrap(~group,ncol=3)

       plot_list[["p_per_replicate_splitByGroup"]] <- p_per_replicate_splitByGroup

       p_per_group_splitByReplicate <- ggplot(df, aes(x=perc,
                                                          y=cov,
                                                          group=sample,
                                                          color=group)) +
           geom_line() +
           labs(x = "Gene body percentile 5' -> 3'",
                y = "Averaged normalised coverage") +
           scale_color_manual(values=define.group.palette(num_groups)) +
           ylim(0,1) +
           theme_bw() +
	   theme(legend.position="top") +
	   guides(color=guide_legend(ncol=3)) +
           facet_wrap(~replicate,ncol=3)

       plot_list[["p_per_group_splitByReplicate"]] <- p_per_group_splitByReplicate

    }
  return(plot_list)
}
##
##DEhelper.strandspecificity: get the strand specificity from the qc and display them
##
DEhelper.strandspecificity <- function(){

    # logs folder
    if(!file.exists(SHINYREPS_INFEREXPERIMENT_LOGS)) {
        return("Strand specificity statistics not available")
    }
    
    filelist <- list.files(path=SHINYREPS_INFEREXPERIMENT_LOGS, full.names=TRUE)
    strandspecifity <- lapply(filelist, read.table, sep=":", skip=3, header=FALSE, row.names=1, blank.lines.skip=TRUE)
    strandspecifity <- do.call(cbind, strandspecifity)

    samplenames <- gsub("_inferexperiment.txt", "", basename(filelist))

    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file )

        # replace files names with nicer sample names given in targets file
	      # if sample is missing in targets file, use reduced file name
        samplenames <- sapply(samplenames, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,   
                                                                targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                                gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {

        if(!is.na(SHINYREPS_PREFIX)) {
            samplenames <- gsub(paste0("^",SHINYREPS_PREFIX), "", samplenames)
        }
    }

    colnames(strandspecifity) <- samplenames 
    rownames(strandspecifity) <- c("ambiguous", "sense", "antisense")
    
    # sort alphabetically for report
    strandspecifity <- as.data.frame(t(strandspecifity)[order(rownames(t(strandspecifity))),])

    # add another column specifying the strandedness based on which strandedness is expected [no|yes|reverse]
    if (SHINYREPS_STRANDEDNESS == "yes") {
        strandspecifity$`strandedness [%]` <- 100*round(strandspecifity$sense/(strandspecifity$sense+strandspecifity$antisense),digits=4)
    } else {
        if (SHINYREPS_STRANDEDNESS == "reverse") {
            strandspecifity$`strandedness [%]` <- 100*round(strandspecifity$antisense/(strandspecifity$sense+strandspecifity$antisense),digits=4)
        }
    }

    kable(strandspecifity, output=F, format="markdown", align=c("c"))
}

##
## DEhelper.GO_Enrichment: get the GO enrichment results and display them
##
DEhelper.GO_Enrichment <- function(){

    #csv file
    if(!file.exists(SHINY_GO)){
        return("GO enrichment statistics not available")
    }
}


##
## DEhelper.Subread: parse Subread summary stats and create a md table
##
DEhelper.Subread <- function() {
    
    FOLDER <- SHINYREPS_SUBREAD
    SUFFIX <- paste0(SHINYREPS_SUBREAD_SUFFIX, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return("Subread statistics not available")
    }
    
    # create a matrix using feature names as rownames, sample names as colnames
    x <- sapply(list.files(FOLDER, pattern=SUFFIX), function(f) {
        
        f <- file(paste0(FOLDER, '/', f))
        l <- readLines(f)
        close(f)
        
	## get all interesting rows, extract count and category names
        assigned.and.unassigned <- l[grep("Assigned|Unassigned",l)]
        category.counts <- as.numeric(sapply(assigned.and.unassigned, function(i) {strsplit(i,"\t")[[1]][2]}))
        names(category.counts) <- sapply(assigned.and.unassigned, function(i) {strsplit(i,"\t")[[1]][1]})
        return(category.counts)
    })
    
    # correct column names
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
    
    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
        # replace files names with nicer sample names given in targets file
	# if sample is missing in targets file, use reduced file name
        colnames(x) <- sapply(colnames(x), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                                targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                                                gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {

        if(!is.na(SHINYREPS_PREFIX)) {
            colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX), "", colnames(x))
        }
    }

    # create md table (omitting various values that are 0 for now)
    # from x we remove the ones which are unmapped to calculate percentages
    # only for the mapped ones
    x <- x[rownames(x) != "Unassigned_Unmapped", ]
    x <- rbind(total=x, colSums(x))
    rownames(x)[nrow(x)] <- "total"
    df <- data.frame(assigned=paste0(format(x["Assigned", ], big.mark=","), " (", format((x["Assigned", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)"), 
		     unassigned_ambig=paste0(format(x["Unassigned_Ambiguity", ], big.mark=","), " (", format((x["Unassigned_Ambiguity", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)"), 
		     unassigned_multimap=paste0(format(x["Unassigned_MultiMapping", ], big.mark=","), " (", format((x["Unassigned_MultiMapping", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)"), 
		     unassigned_nofeat=paste0(format(x["Unassigned_NoFeatures", ], big.mark=","), " (", format((x["Unassigned_NoFeatures", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)"))
    rownames(df) <- colnames(x)

    # sort alphabetically for plotting
    df <- df[order(rownames(df)),]
    kable(df, align=c("r", "r", "r", "r"), output=F, format="markdown")
    
}

##
## extract the intron/exon and intergenic regions from the qualimap report
##
DEhelper.Qualimap <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_QUALIMAP_LOGS)) {
        return("Read distribution statistics not available")
    }  
    
    QC <- SHINYREPS_QUALIMAP_LOGS    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(QC, pattern="rnaseq_qc_results.txt$", recursive=T)
    if(length(samples) == 0) {
        return("Qualimap report not available")
    }
    sample_names <- gsub("_counts_qualimap", "", dirname(samples))
    
    samples <- paste0(QC, "/", samples)
    samples <- sapply(samples, function(f){
      l <- readLines(f)
      l <- l[grep("exonic|intronic|intergenic", l)] 
      l <- strsplit(gsub(" ", "", l), "=")
      names(l) <- c("exonic","intronic","intergenic")
      l <- sapply(l, function(x){
           as.numeric(gsub("\\%\\)", "", strsplit(x[[2]], "\\(")[[1]][2]))    
      })
    })
    colnames(samples) <- sample_names
    # qualimap outputs an additional category "overlapping_exon", which is part of intron or intergenic,
    # but since it is not specified, which of the two, we will not display this category here.

    #Todo fix the naming scheme, the levels and then put it in the barplot 
    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "",targets$file)

        # replace files names with nicer sample names given in targets file
        # if sample is missing in targets file, use reduced file name
        colnames(samples) <- sapply(colnames(samples), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,  
                                                                            targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                                            gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {
        if(!is.na(SHINYREPS_PREFIX)) {
            colnames(samples) <- gsub(paste0("^",SHINYREPS_PREFIX), "", colnames(samples))
        }
    }
    
    # sort according to the groups and samplenames in the targets or alphabetically
    samples <- melt(samples, value.name="perc")
    colnames(samples) <- c("class", "sample", "perc")

    # sort the samples either by the group order of the targets file or alphabetically
    if(file.exists(SHINYREPS_TARGET) & SHINYREPS_SORT_ALPHA==FALSE){
       sample_order   <- targets$sample[order(paste0(targets$group, targets$sample), decreasing=TRUE)]
       samples$sample <- factor(samples$sample, levels = sample_order)
    }else{
       samples$sample <- factor(samples$sample, levels = sort(unique(samples$sample), decreasing=TRUE))  
    }

    samples$class <- factor(samples$class, levels=c("intergenic", "intronic", "exonic" )) 

    # plot
    p <- ggplot(samples, aes(x=sample,y=perc,fill=class)) +
           geom_col(width = 0.8) +
           labs(x = "",
                y = "% of mapped reads") +
           theme_bw() +
           theme(axis.text.y = element_text(size = 8),
                 legend.position = "top") +
           guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) +
           scale_fill_brewer(palette = "Dark2") +
           coord_flip()
  return(p)
}

##
##DEhelper.insertsize: get the insertsize from the qc and display mean and sd 
##
DEhelper.insertsize <- function(){

	if (SHINYREPS_PAIRED == "yes") {
		filelist <- list.files(path=SHINYREPS_INSERTSIZE,full.names=TRUE, pattern="insertsizemetrics.tsv$")
		insertsizes <- lapply(filelist, read.table, sep="\t", header=TRUE, nrow=1)
		insertsizes <- do.call(rbind, insertsizes)
		samplenames <- basename(filelist)
		
		if(file.exists(SHINYREPS_TARGET)){
		  
		  # get target names
		  targets <- read.delim(SHINYREPS_TARGET)
		  targets$sample_ext <- gsub("\\..*$", "",targets$file)
		  
		  # replace files names with nicer sample names given in targets file
		  # if sample is missing in targets file, use reduced file name
		  samplenames <- sapply(samplenames, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,  
		                                                                      targets[sapply(targets$sample_ext, grepl, i),"sample"], 
		                                                                      gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
		} else {
		  if(!is.na(SHINYREPS_PREFIX)) {
		    samplenames <- gsub(paste0("^",SHINYREPS_PREFIX), "", samplenames)
		    samplenames <- gsub("_insertsizemetrics.tsv","", samplenames)
		  }
		}

		rownames(insertsizes) <- samplenames 
		insertsizes <- insertsizes[,c("MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE", "STANDARD_DEVIATION")]
		colnames(insertsizes) <- c("Median", "Mean", "SD")
		kable(insertsizes, output=F, align=c("l"), format="markdown")
	}
}


# Helper to plot the insertsize histogram equivalent to the one from picard
# Input is the Picard generated metrics file
DEhelper.insertsize.helper <- function(metricsFile){
  #find the start of our metrics informatioun 
  startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)
  
  firstBlankLine=0
  
  for (i in 1:length(startFinder)) {
    if (startFinder[i] == "") {
      if (firstBlankLine == 0) {
        firstBlankLine = i+1
      } else {
        secondBlankLine = i+1
        break
      }
    }
  }
  
  histogram <- read.table(metricsFile, header=TRUE, sep="\t", skip=secondBlankLine, comment.char="", quote='', check.names=FALSE)
  
  ## The histogram has a fr_count/rf_count/tandem_count for each metric "level"
  ## This code parses out the distinct levels so we can output one graph per level
  headers <- sapply(sub(".fr_count","",names(histogram),fixed=TRUE), "[[" ,1)
  headers <- sapply(sub(".rf_count","",headers,fixed=TRUE), "[[" ,1)
  headers <- sapply(sub(".tandem_count","",headers,fixed=TRUE), "[[" ,1)
  
  ## Duplicate header names could cause this to barf.  But it really shouldn't when we have "All_reads.fr_count" and 
  ## "All_reads.rf_count" for example.  Not sure why this would fail, but I care.
  if (any(duplicated(headers))) {
    levels = unique(headers[2:length(headers)]);
  } else {
    levels <- c()
    for (i in 2:length(headers)) {
      if (!(headers[i] %in% levels)) {
        levels[length(levels)+1] <- headers[i]
      }
    }
  }
  
  title_info <- basename(metricsFile)  
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET)
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    title_info <- sapply(title_info, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,   
                                                    targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                    gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
  } else {
    if(!is.na(SHINYREPS_PREFIX)) {
      title_info <- gsub(paste0("^",SHINYREPS_PREFIX), "", title_info)
    }
    title_info <- gsub("_insertsizemetrics.tsv$","", title_info)
  }

  #we get the histogram which has the names of the levels e.g. all_reads and readgroups/sample groups depending on
  #accumulation level which was used.
  #the colnames of histogram are something like all.read.fr_count, all.read.rf_count etc.
  #to get the whole shebang into a wider format we have to add the information
  hist_long <- reshape2::melt(histogram, id.var = "insert_size") %>% 
    extract(col=variable,
            into=c("group", "counttype" ),
            regex='([^\\.]+)\\.([^\\.]+)') %>% 
    dplyr::rename( amount = value)
  #we also have to add the comulative sum per group to the whole shebang
  hist_long <- hist_long %>% group_by(group) %>% 
    arrange( desc(insert_size)) %>% 
    mutate( cumulative = (cumsum(amount)/sum(amount))*max(amount))
  #1. Create one plot per group (all_reads etc)
  #2. save the plots in a list
  #3. use arrangeGrob to arrange
  hist_long <- split(hist_long, hist_long$group)
  hist_plots <- lapply(hist_long, function(hist_data){
    p <- ggplot(hist_data, aes(x    = insert_size,
                               y    = amount,
                               fill = counttype)) +
      geom_bar( stat = "identity", 
                position = "identity",
                alpha=0.7) +
      scale_fill_hue( c=50, l=40) +
      geom_line(aes(x = insert_size,
                    y = cumulative),
                linetype = "dashed",
                color    = "grey") +
      scale_y_continuous( name="# reads",
                          sec.axis = sec_axis(~./max(hist_data$amount),
                                                    name = "Cumulative fraction of reads > insert size")) +
                           xlab("insert size in bp") +
                           theme_bw() +
      theme(legend.justification=c(1,1), legend.position=c(1,1), 
            legend.background = element_rect(colour = "transparent", fill = "transparent")) + 
      facet_grid(~group)
      return(p)
})
 
 hist_plot <- arrangeGrob(grobs = hist_plots,
                           top = textGrob(title_info)) 
 return(hist_plot)

}

## 
## DEhelper.subchunkify: small function stolen from here 
##http://michaeljw.com/blog/post/subchunkify/
## to dynamically create chunks and adjust their size accordingly.
##
DEhelper.subchunkify <- function(g, fig_height=7, fig_width=5) {
  g_deparsed <- paste0(deparse(
    function() {grid.draw(g)}
  ), collapse = '')
  
  sub_chunk <- paste0("`","``{r sub_chunk_", floor(runif(1) * 10000),
                      ", fig.height=", fig_height,
                      ", fig.width=", fig_width,
                      ", echo=FALSE}","\n(", 
                      g_deparsed, ")()",
                      "\n`","``")
  
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}


##
##DEhelper.insertsize.plot: get the insertsize histograms and display them 
##
DEhelper.insertsize.plot <- function(){
  # logs folder
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),
                                     error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
  }
  if (SHINYREPS_PAIRED == "yes" &
      length(list.files(path = SHINYREPS_INSERTSIZE,
                        pattern = "insertsizemetrics.tsv$")) > 0) {
    samples <- list.files(path = SHINYREPS_INSERTSIZE,
                          full.names = TRUE,
                          pattern = "insertsizemetrics.tsv$")
    #we generate the plots
    insert_plots <- lapply(samples, 
                           DEhelper.insertsize.helper)
    return(arrangeGrob(grobs = insert_plots,
                       ncol = SHINYREPS_PLOTS_COLUMN))
  }else{
    return("No insertsize histograms available.")
  }
}


##
## DEhelper.Trackhub: display the UCSC trackhub URL
##
DEhelper.Trackhub <- function() {
    
    # output file with trackhub URL
    if(!file.exists(SHINYREPS_TRACKHUB_DONE)) {
        return("UCSC GB Trackhub URL file not available")
    }
    
    # Trackhub URL is second line of file
    url <- scan(SHINYREPS_TRACKHUB_DONE, skip=0, nlines=1, what='character')
    if (grepl("hub.txt", url)) {
        return(url)
    } else {
        return("UCSC GB Trackhub URL not available")
    }
}


##
## DEhelper.cutadapt: get trimming statistics from the Cutadapt folder and display them
## 
#' @param targetsdf targets data.frame or character with file path to targets object
#' @param colorByFactor character with column name of sample table to be used for coloring the plot. Coloring by filename if NULL. 
#' @param sampleColumnName character with column name(s) of targets table containing file names
#' @param plotfun define function to be used for plotting
#' @param labelOutliers logical, shall outlier samples be labeled
#' @param outlierIQRfactor numeric, factor is multiplied by IQR to determine outlier
#'
#' @return plot cutadapt statistics as side effect
DEhelper.cutadapt <- function(targetsdf=SHINYREPS_TARGET, colorByFactor="group", sampleColumnName =c("file"), 
                              plotfun=DEhelper.cutadapt.plot, labelOutliers=T, outlierIQRfactor=1.5
                              ){
  
  # logs folder
  if(!all(sapply(SHINYREPS_CUTADAPT_STATS, file.exists))) {
    return(paste("Cutadapt statistics not available for", names(which(!sapply(SHINYREPS_CUTADAPT_STATS, file.exists)))))
  }
  
  x <- list.files(SHINYREPS_CUTADAPT_STATS,pattern='*cutadapt.log$',full.names=TRUE) 
  
  # get Command line parameters of first file
  cutadaptpars <- system(paste("grep \"Command line parameters\"", x[1]), intern=T)
  
  paired <- grepl("(-p )|(--paired-output )", cutadaptpars) # output for R2 available?
  
  x <- sapply(x, function(f) { 
    
    trimmed.R1.perc <- trimmed.R2.perc <- trimmed.reads.perc <- tooshort.reads.perc <- NULL # initialize with NULL in case not needed
    
    if(paired) { # log lines slightly differ dependent on se or pe
      total.reads <- system(paste("grep \"Total read pairs processed\"", f, "| awk '{print $5}'"), intern=TRUE)
      total.reads <- gsub(",", "", total.reads)
      trimmed.R1.perc <- system(paste("grep \"Read 1 with adapter\"", f, "| awk '{print $6}'"), intern=TRUE)
      trimmed.R1.perc <- gsub("\\(|\\)|\\%", "", trimmed.R1.perc)
      trimmed.R2.perc <- system(paste("grep \"Read 2 with adapter\"", f, "| awk '{print $6}'"), intern=TRUE)
      trimmed.R2.perc <- gsub("\\(|\\)|\\%", "", trimmed.R2.perc)
      
    } else {
      total.reads <- system(paste("grep \"Total reads processed\"", f, "| awk '{print $4}'"), intern=TRUE)
      total.reads <- gsub(",", "", total.reads)
      trimmed.reads.perc <- system(paste("grep \"Reads with adapters\"", f, "| awk '{print $5}'"), intern=TRUE)
      trimmed.reads.perc <- gsub("\\(|\\)|\\%", "", trimmed.reads.perc)
    }
    
    tooshort.reads.perc <- system(paste("grep \"that were too short\"", f, "| awk '{print $7}'"), intern=TRUE)
    tooshort.reads.perc <- gsub("\\(|\\)|\\%", "", tooshort.reads.perc)
    
    # trimming of each adapter
    adapters <- system(paste("grep Sequence:", f, "| awk '{print $9}'"), intern=T)
    adapters.perc <- round(100*(as.numeric(adapters) / as.numeric(total.reads)),1)
    adapterprime <- gsub(";", "", system(paste("grep Sequence:", f, "| awk '{print $5}'"), intern=T))
    
    names(adapters.perc) <- gsub(" *=== *", "", system(paste("grep \"=== .*Adapter\"", f), intern=T))
    namespart1 <- gsub("First read:.*", "R1_", names(adapters.perc))
    namespart1 <- gsub("Second read:.*", "R2_", namespart1)
    namespart2 <- gsub("^.*Adapter", "Adapter", names(adapters.perc))
    names(adapters.perc) <- paste0(if(paired) {namespart1} else {""}, adapterprime, namespart2)
    
    ## add trimmed reads for each adapter here
    return(c("total reads"=total.reads, trimmed_R1=trimmed.R1.perc, trimmed_R2=trimmed.R2.perc, 
             trimmed=trimmed.reads.perc, "too short"=tooshort.reads.perc, adapters.perc))
  })
  
  # transpose dataframe
  x.df <- as.data.frame(t(x), make.names=F, stringsAsFactors = F) 
  x.df <- as.data.frame(lapply(x.df, as.numeric))
  colnames(x.df) <- rownames(x)
  
  # use cutadapt call from first log file for naming some of the unnamed adapters 
  cutadaptpars <- unlist(strsplit(cutadaptpars, split=" ")) 
  indexAdapter <- grep("(^-a$)|(--adapter)|(^-g$)|(--front)|(^-A$)|(^-G$)", cutadaptpars) # index of all adapters applied
  indexAdapterSelected <- indexAdapter[grep("[ACGT].[[:digit:]]*}", cutadaptpars[indexAdapter+1])] # select e.g. polyA, polyT
  
  # rename those adapters columns trimmed by -a commands 
  if (length(indexAdapterSelected)>0) {
    colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)] <- 
      paste0(gsub("Adapter.*$", "", colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)]), cutadaptpars[indexAdapterSelected+1])
  }
  
  #reduce length of file names 
  row.names(x.df) <- basename(colnames(x))
  x.df$filename_unmod <- factor(row.names(x.df))
  if(!is.na(SHINYREPS_PREFIX)) {
    row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))
  }
  row.names(x.df) <- gsub("\\.cutadapt\\.log$", "", row.names(x.df))
  if(nrow(x.df)>1){
    if(is.na(SHINYREPS_PREFIX)) {row.names(x.df)  <- gsub(lcPrefix(row.names(x.df) ), "", row.names(x.df) )}
    row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
  }
  
  # passing the different factors given in targetsdf to x.df which was created from cutadapt file names 
  if(!is.null(colorByFactor)) { # add information to x.df
    
    if(is.null(targetsdf)) {stop("If 'colorByFactor' is given you must also provide 'targetsdf'!")}
    
    if(!is.data.frame(targetsdf) && is.character(targetsdf) && file.exists(targetsdf)){
      targetsdf <- read.delim(targetsdf)
    } 
      
    if(length(sampleColumnName)>1) { # melt in case of multiple file name columns (as for ChIP-Seq)
      targetsdf <- targetsdf[,unique(c(colorByFactor, sampleColumnName, "sample"))]
      targetsdf <- reshape2::melt(targetsdf, measure.vars=sampleColumnName, value.name = "filename") 
      for (i in colorByFactor) {targetsdf[, i] <- paste0(targetsdf[, i], " (", targetsdf$variable, ")")}
      targetsdf[,c(colorByFactor, "filename")] <- lapply(targetsdf[,c(colorByFactor, "filename")], factor)
      
    } else {
      targetsdf$filename <- targetsdf[,sampleColumnName]
    }
    
    targetsdf$filename <- gsub("\\..*$", "", targetsdf$filename ) # shorten filename suffix
    index <- sapply(targetsdf$filename, grep, x.df$filename_unmod, ignore.case = T) # grep sample name in file names
    if(is.list(index)) {
      #targetsdf <- targetsdf[sapply(index, length) ==1, ] # remove targetsdf entries not (uniquely) found in x.df
      targetsdf <- targetsdf[sapply(index, length)!=0,] # remove targetsdf entries not found in x.df
      index <- sapply(targetsdf$filename, grep, x.df$filename_unmod, ignore.case = T) # redo grep sample name in file names
    }
    
    if(!identical(sort(unname(unlist(index))), 1:nrow(x.df))) {
      stop("There seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
    }
    
    # x.df <- data.frame(x.df[unlist(index),], targetsdf, check.names =F) ##temp
    x.df <- data.frame(x.df[unlist(t(index)),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    if("sample" %in% colnames(x.df) && !any(duplicated(x.df$sample))) { # use sample column as identifier if present and unique
      x.df$filename <- x.df$sample
      row.names(x.df) <- x.df$sample } 
    
    if(any(!colorByFactor %in% colnames(x.df))) {
      if(all(!colorByFactor %in% colnames(x.df))) {
        cat("\nNone of the column names given in colorByFactor is available. Perhaps sample names are not part of fastq file names? Using filename instead.")
        colorByFactor <- "filename"
      } else { # one plot each element of colorByFactor
        cat("\n", colorByFactor[!colorByFactor %in% colnames(x.df)], "not available. Using", colorByFactor[colorByFactor %in% colnames(x.df)], "instead.")
        colorByFactor <- colorByFactor[colorByFactor %in% colnames(x.df)]
      }
    }
  } else {
    x.df$filename <- row.names(x.df)
    colorByFactor <- "filename"
  }
  
  # melt data frame for plotting
  x.melt <- reshape2::melt(x.df, measure.vars=c(grep("trimmed", colnames(x.df), value=T), 
                                                "too short", 
                                                grep("(Adapter)|(})", colnames(x.df), value=T)), variable.name="reads")
  # everything which is not a value should be a factor
  
  # one plot for each element of colorByFactor
  violin.list <- lapply(colorByFactor, plotfun, data=x.melt, labelOutliers=labelOutliers, outlierIQRfactor=outlierIQRfactor) # "colorByFactor" is submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  DT::datatable(x.df[,c(colorByFactor, "total reads", 
                        grep("trimmed", colnames(x.df), value=T),
                        "too short", 
                        grep("(Adapter)|(})", colnames(x.df), value=T))], 
                options = list(pageLength= 20))
}


# plotting function for DEhelper.cutadapt 
DEhelper.cutadapt.plot <- function(data, color.value, labelOutliers=T, outlierIQRfactor=1.5){
  
  is_outlier <- function(x) { # function for identification of outlier
    if(IQR(x)!=0) {
      return(x < quantile(x, 0.25) - outlierIQRfactor * IQR(x) | x > quantile(x, 0.75) + outlierIQRfactor * IQR(x))
    } else {
      return(x < mean(x) - outlierIQRfactor * mean(x) | x > mean(x) + outlierIQRfactor * mean(x))
    }
  }
  
  data <- data %>%
    dplyr::group_by(reads) %>%
    dplyr::mutate(outlier=is_outlier(value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outlier=ifelse(outlier,filename,as.numeric(NA))) %>%
    as.data.frame()
  
  ylab <- "% reads"
  
  # prepare palette of appropriate length according to the different factors given in colorByFactor
  colourCount = length(unique(data[,color.value]))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  p <- ggplot(data, aes_string(x="reads",
                               y="value",
                               color=color.value ))+
    geom_quasirandom(groupOnX=TRUE) +
    geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA)  
  if(labelOutliers) {p <- p + ggrepel::geom_text_repel(data=. %>% filter(!is.na(outlier)), aes(label=filename), show.legend=F)}
  p <- p + scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
    ylab(ylab) +
    xlab("") +
    theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1)) 
  
  return(p)
}

## 
#' rMATS differential splicing analysis
#'
#' @param No Params needed

DEhelper.rmats <- function() {
    
    # source modified maser functions
    # by default, plotTranscriptsMod.R allows human chromosome annotation only. 
    # volcano plot: the axis labels are switched and the legend is messed up.
    source(file.path(SHINYREPS_MASER_SCRIPTS,"plotTranscriptsMod.R"))
    source(file.path(SHINYREPS_MASER_SCRIPTS,"volcanoMod.R"))
    SHINYREPS_MASER_FDR = as.numeric(SHINYREPS_MASER_FDR)
    SHINYREPS_MASER_DPSI = as.numeric(SHINYREPS_MASER_DPSI)

    # The function topEvents() allows to select statistically significant events given a FDR cutoff and minimum PSI change. 
    rmats_top <- maser::topEvents(rmats_filt, fdr = SHINYREPS_MASER_FDR, deltaPSI = SHINYREPS_MASER_DPSI)

    # Plot the distribution of splicing events
    plotspldist <- maser::splicingDistribution(rmats_filt, fdr = SHINYREPS_MASER_FDR, deltaPSI = SHINYREPS_MASER_DPSI)
    print(plotspldist)

    for(e in c("A3SS", "A5SS", "SE", "RI", "MXE")) {
        if(nrow(slot(rmats_top, paste0(e,"_events"))) > 0) {
            cat("\n\n", fill=T)
            cat("#### results for", e, fill=T)  
            cat("\n\n", fill=T)
           
            top <- summary(rmats_top, type=e)
	    top <- top[order(top$FDR), !colnames(top) %in% c("GeneID")]
	    colnames(top) <- plyr::mapvalues(colnames(top), from=c("ID", "PValue", "IncLevelDifference"), to=c("rMATS_ID", "pval", "IncLevelDiff"))
	    top$pval <- signif(top$pval, 3)
	    top$FDR <- signif(top$FDR, 3)
	    cat("\n\nSplicing events have genomic locations, raw counts, PSI levels and statistics (p-values) regarding their differential splicing between conditions.")
	    cat("\nTop significant", e, ":\n\n")
            DT::datatable(top, options = list(pageLength = 10), rownames = F)
	    plotPCA <- maser::pca(rmats_filt, type = e)
            plotVolcano <- volcanoMod(rmats_filt, fdr = SHINYREPS_MASER_FDR, deltaPSI = SHINYREPS_MASER_DPSI, type = e)
            gridExtra::grid.arrange(grobs=list(plotPCA, plotVolcano), ncol=1)
    
            top1 <- maser::geneEvents(rmats_filt, geneS = top$geneSymbol[1], fdr = SHINYREPS_MASER_FDR, deltaPSI = SHINYREPS_MASER_DPSI)
    
            ## Display affected transcripts and PSI levels
            plotTranscriptsMod(events=top1, type = e, event_id = top$rMATS_ID[1], gtf = ens_gtf, zoom = F, show_PSI = TRUE, title =top$geneSymbol[1]) 
    
        }
    }  
}

##		      
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    tryCatch({
        ver <- read.delim(file=SHINYREPS_TOOL_VERSIONS)
        colnames(ver) <- c("Tool name","Environment", "Version")
        kable(as.data.frame(ver),output=F, format="markdown")
    }, error=function(e) cat("tool versions not available.\n", fill=TRUE))
}
