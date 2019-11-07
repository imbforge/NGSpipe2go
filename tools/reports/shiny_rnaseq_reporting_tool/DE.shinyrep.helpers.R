##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("edgeR")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("reshape2")
library("ggrepel")
library("VennDiagram")
library("grid")
library("knitr")        # for markdown output
library("plotly")
library("GeneOverlap")
library("pheatmap")
library("viridis")
library("gridExtra")
library("dplyr")
library("tidyr")
library("forcats")
#' loadGlobalVars: read configuration from bpipe vars
#'
#' @param f - a file defining multiple variables for reporting to run. 
#'
#' @return A set of variables that are mentioned in input file, e.g.
#'         SHINYREPS_PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/butter/imb_butter_2016_10_alina_rnaseq_u2os_oe/"
#'         SHINYREPS_ORG <- "human"
#'         SHINYREPS_DB <- "hg38"
#'              
#' @description File content should be:
#'              SHINYREPS_PROJECT=/fsimb/groups/imb-bioinfocf/projects/butter/imb_butter_2016_10_alina_rnaseq_u2os_oe/
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
    p <- plotPCA(rld, intgroup=colnames(colData(dds))[1])
    print(p + 
          scale_color_manual(values=brewer.pal(9,"Set1")[1:length(levels(colData(dds)[,"group"]))]) +
          geom_text_repel(aes(label=rownames(colData(dds))), show.legend=FALSE) + 
          theme_bw())
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
    p <- plotPCA(rlog(pairwise.dds[[i]]), intgroup=colnames(colData(pairwise.dds[[i]]))[1])
    print(p + 
          scale_color_manual(values=brewer.pal(9,"Set1")[1:2]) + 
          geom_text_repel(aes(label=rownames(colData(pairwise.dds[[i]]))), show.legend=FALSE) + 
          theme_bw())
}


## DEhelper.cluster
#' Heatmap of top variant 'n' genes of the counts-per-milion table.
#'
#' @param n - amount of transcripts/rows to be plotted
#'            [default = 40]
#' @param rld - rlog transformed DESeq2 object
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
#'           DEhelper.cluster(n=10)
#'           
#DEhelper.DESeq2.cluster <- function(n=25) {
#    rows <- order(apply(assay(rld), 1, sd), decreasing=TRUE)[1:n]
#    hmcol  <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#    heatmap.2(assay(rld)[rows, ], col=hmcol, trace="none", margin=c(10, 6), scale="none")
#}
DEhelper.DESeq2.cluster.sd <- function(n=40) {
        
        # extract assay and over-write gene_id in rownames with gene_name
        assay.rld <- assay(rld)
        rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(rld)), gtf$gene_id)]

        # pick top n most variable genes
        select.highestSD <- order(apply(assay.rld,1,sd),decreasing=TRUE)[1:n]

        # set color scheme
        col <- colorRampPalette(brewer.pal(9,"GnBu"))(255)

        # extract information for legend
	if (length(add_factors)==0) {
	        legend.df <- data.frame(group=colData(rld)[,c("group")],row.names=rownames(colData(rld)))
	} else {
        	legend.df <- as.data.frame(colData(rld)[,c("group",add_factors)])
	}

        # fix group colors for legend and possible first additional factor if available
        if (length(add_factors)==0) {
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(rld)[,"group"]))],
                                                       unique(colData(rld)[,"group"])))
        } else {
                if (length(unique(colData(rld)[,add_factors[1]])) <= 8) {
                        mypalette <- brewer.pal(8,"Dark2")[1:length(unique(colData(rld)[,add_factors[1]]))]
                } else {
                        mypalette <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(colData(rld)[,add_factors[1]])))
                }
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(rld)[,"group"]))],
                                                       unique(colData(rld)[,"group"])),
                                      subject = setNames(mypalette,
                                                         unique(colData(rld)[,add_factors[1]])))
                names(legend_colors) <- c("group",add_factors[1])
        }

        # plot heatmap
        pheatmap(assay.rld[select.highestSD,],         
                 cluster_rows=TRUE,
                 cluster_cols=TRUE,
                 show_rownames=TRUE,
                 annotation_col=legend.df,
                 annotation_colors=legend_colors,
                 color=col,
                 border_color=NA,
                 main=paste("Normalized expression values of",n,"most variable genes"),
		 fontsize_row=5,
                 fontsize=6,
                 treeheight_row=20,
                 treeheight_col=20,
		 annotation_names_col=FALSE)
}


## DEhelper.DESeq2.cluster.sd.pairwise
#' Heatmap of top variant 'n' genes of the regularized log transformed counts-per-milion table
#' in a pairwise comparison of two groups
#'
#' @param i - iterator to select contrast from conts [default = 1]
#' @param n - amount of transcripts/rows to be plotted
#'            [default = 40]
#' @param pairwise.dds - DESeq2 object
#'
#' @return A plot (base plotting) object.
#'
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           DEhelper.DESeq2.cluster.sd.pairwise(i=1,n=10)
#'           
DEhelper.DESeq2.cluster.sd.pairwise <- function(i=1,n=40) {

        # extract assay and over-write gene_id in rownames with gene_name
        assay.rld <- assay(rlog(pairwise.dds[[i]]))
        rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(pairwise.dds[[i]])), gtf$gene_id)]

        # pick top n most variable genes
        select.highestSD <- order(apply(assay.rld,1,sd),decreasing=TRUE)[1:n]

        # set color scheme
        col <- colorRampPalette(brewer.pal(9,"GnBu"))(255)

        # extract information for legend
        if (length(add_factors)==0) {
                legend.df <- data.frame(group=colData(pairwise.dds[[i]])[,c("group")],row.names=rownames(colData(pairwise.dds[[i]])))
        } else {
                legend.df <- as.data.frame(colData(pairwise.dds[[i]])[,c("group",add_factors)])
        }

        # fix group colors for legend and possible first additional factor if available
        if (length(add_factors)==0) {
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(pairwise.dds[[i]])[,"group"]))],
                                                       unique(colData(pairwise.dds[[1]])[,"group"])))
        } else {
                if (length(unique(colData(pairwise.dds[[i]])[,add_factors[1]])) <= 8) {
                        mypalette <- brewer.pal(8,"Dark2")[1:length(unique(colData(pairwise.dds[[i]])[,add_factors[1]]))]
                } else {
                        mypalette <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(colData(pairwise.dds[[i]])[,add_factors[1]])))
                }
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(pairwise.dds[[i]])[,"group"]))],
                                                       unique(colData(pairwise.dds[[i]])[,"group"])),
                                      subject = setNames(mypalette,
                                                         unique(colData(pairwise.dds[[i]])[,add_factors[1]])))
                names(legend_colors) <- c("group",add_factors[1])
        }

        # plot heatmap
        pheatmap(assay.rld[select.highestSD,],
                 cluster_rows=TRUE,
                 cluster_cols=TRUE,
                 show_rownames=TRUE,
                 annotation_col=legend.df,
                 annotation_colors=legend_colors,
                 color=col,
                 border_color=NA,
                 main=paste("Normalized expression values of",n,"most variable genes"),
                 fontsize_row=5,
                 fontsize=6,
                 treeheight_row=20,
                 treeheight_col=20,
		 annotation_names_col=FALSE)
}

## DEhelper.DESeq2.cluster.mean
#' Heatmap of top 'n' highest mean genes of the regularized log transformed counts-per-milion table.
#'
#' @param n - amount of transcripts/rows to be plotted
#'            [default = 40]
#' @param rld - rlog transformed DESeq2 object
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
#'           DEhelper.DESeq2.cluster.mean(n=10)
#'           
DEhelper.DESeq2.cluster.mean <- function(n=40) {

        # extract assay and over-write gene_id in rownames with gene_name
        assay.rld <- assay(rld)
        rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(rld)), gtf$gene_id)]

        # pick top n most variable genes
        select.highestMean <- order(rowMeans(assay.rld),decreasing=TRUE)[1:n]

        # set color scheme
	col <- rev(heat.colors(255))

	# extract information for legend
        if (length(add_factors)==0) {
                legend.df <- data.frame(group=colData(rld)[,c("group")],row.names=rownames(colData(rld)))
        } else {
                legend.df <- as.data.frame(colData(rld)[,c("group",add_factors)])
        }

        # fix group colors for legend and possible first additional factor if available
        if (length(add_factors)==0) {
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(rld)[,"group"]))],
                                                       unique(colData(rld)[,"group"])))
        } else {
                if (length(unique(colData(rld)[,add_factors[1]])) <= 8) {
                        mypalette <- brewer.pal(8,"Dark2")[1:length(unique(colData(rld)[,add_factors[1]]))]
                } else {
                        mypalette <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(colData(rld)[,add_factors[1]])))
                }
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(rld)[,"group"]))],
                                                       unique(colData(rld)[,"group"])),
                                      subject = setNames(mypalette,
                                                         unique(colData(rld)[,add_factors[1]])))
                names(legend_colors) <- c("group",add_factors[1])
        }

        # plot heatmap
        pheatmap(assay.rld[select.highestMean,],
                 cluster_rows=FALSE,
                 cluster_cols=TRUE,
                 show_rownames=TRUE,
                 annotation_col=legend.df,
                 annotation_colors=legend_colors,
                 col=col,
                 border_color=NA,
                 main=paste("Normalized expression values of",n,"genes with highest mean"),
                 fontsize_row=5,
                 fontsize=6,
                 treeheight_row=20,
                 treeheight_col=20,
		 annotation_names_col=FALSE)
}

## DEhelper.DESeq2.cluster.mean.pairwise
#' Heatmap of top 'n' highest mean genes of the regularized log transformed counts-per-milion table
#' in a pairwise comparison of two groups
#'
#' @param i - iterator to select contrast from conts [default = 1]
#' @param n - amount of transcripts/rows to be plotted
#'            [default = 40]
#' @param pairwise.dds - DESeq2 object
#'
#' @return A plot (base plotting) object.
#'
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           DEhelper.DESeq2.cluster.mean.pairwise(n=10)
#'           
DEhelper.DESeq2.cluster.mean.pairwise <- function(i=1,n=40) {

        # extract assay and over-write gene_id in rownames with gene_name
        assay.rld <- assay(rlog(pairwise.dds[[i]]))
        rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(pairwise.dds[[i]])), gtf$gene_id)]

        # pick top n most variable genes
        select.highestMean <- order(rowMeans(assay.rld),decreasing=TRUE)[1:n]

        # set color scheme
	col <- rev(heat.colors(255))

        # extract information for legend
        if (length(add_factors)==0) {
                legend.df <- data.frame(group=colData(pairwise.dds[[i]])[,c("group")],row.names=rownames(colData(pairwise.dds[[i]])))
        } else {
                legend.df <- as.data.frame(colData(pairwise.dds[[i]])[,c("group",add_factors)])
        }

        # fix group colors for legend and possible first additional factor if available
        if (length(add_factors)==0) {
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(pairwise.dds[[i]])[,"group"]))],
                                                       unique(colData(pairwise.dds[[1]])[,"group"])))
        } else {
                if (length(unique(colData(pairwise.dds[[i]])[,add_factors[1]])) <= 8) {
                        mypalette <- brewer.pal(8,"Dark2")[1:length(unique(colData(pairwise.dds[[i]])[,add_factors[1]]))]
                } else {
                        mypalette <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(colData(pairwise.dds[[i]])[,add_factors[1]])))
                }
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(pairwise.dds[[i]])[,"group"]))],
                                                       unique(colData(pairwise.dds[[i]])[,"group"])),
                                      subject = setNames(mypalette,
                                                         unique(colData(pairwise.dds[[i]])[,add_factors[1]])))
                names(legend_colors) <- c("group",add_factors[1])
        }

        # plot heatmap
        pheatmap(assay.rld[select.highestMean,],
                 cluster_rows=FALSE,
                 cluster_cols=TRUE,
                 show_rownames=TRUE,
                 annotation_col=legend.df,
                 annotation_colors=legend_colors,
                 col=col,
                 border_color=NA,
                 main=paste("Normalized expression value of",n,"genes with highest mean"),
                 fontsize_row=5,
                 fontsize=6,
                 treeheight_row=20,
                 treeheight_col=20,
		 annotation_names_col=FALSE)
}


## DEhelper.corr: 
#' Heatmap of sample to sample distances.
#'
#' @param rld - rlog transformed DESeq2 object
#' @return A plot (base plotting) object.
#'
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           rld <- rlog(dds)
#'           # rld is taken from environment
#'           DEhelper.DESeq2.corr()
#'           
#DEhelper.DESeq2.corr <- function() {
#    distsRL <- dist(t(assay(rld)))
#    mat <- as.matrix(distsRL)
#    hc  <- hclust(distsRL)
#    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#    heatmap.2(mat, Rowv=as.dendrogram(hc), 
#              symm=TRUE, trace="none", 
#              col = rev(hmcol), margin=c(13, 13))
#}
DEhelper.DESeq2.corr <- function() {

        # extract assay and over-write gene_id in rownames with gene_name
        assay.rld <- assay(rld)
        rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(rld)), gtf$gene_id)]

        # extract sample to sample distances
        sampleDists <- dist(t(assay.rld))
        sampleDistMatrix <- as.matrix(sampleDists)
        
        # set color scheme
	col <- plasma(255)

        # extract information for legend
        if (length(add_factors)==0) {
                legend.df <- data.frame(group=colData(rld)[,c("group")],row.names=rownames(colData(rld)))
        } else {
                legend.df <- as.data.frame(colData(rld)[,c("group",add_factors)])
        }

        # fix group colors for legend and possible first additional factor if available
        if (length(add_factors)==0) {
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(rld)[,"group"]))],
                                                       unique(colData(rld)[,"group"])))
        } else {
                if (length(unique(colData(rld)[,add_factors[1]])) <= 8) {
                        mypalette <- brewer.pal(8,"Dark2")[1:length(unique(colData(rld)[,add_factors[1]]))]
                } else {
                        mypalette <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(colData(rld)[,add_factors[1]])))
                }
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(rld)[,"group"]))],
                                                       unique(colData(rld)[,"group"])),
                                      subject = setNames(mypalette,
                                                         unique(colData(rld)[,add_factors[1]])))
                names(legend_colors) <- c("group",add_factors[1])
        }
	
	# sample to sample distance heatmap w/ 0's on diagonal set to NA
	# since they otherwise shift the scale too much
	sampleDistMatrix.diagNA <- sampleDistMatrix
	sampleDistMatrix.diagNA[sampleDistMatrix.diagNA==0] <- NA

        # plot heatmap
        pheatmap(sampleDistMatrix.diagNA,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 annotation_col=legend.df,
                 annotation_colors=legend_colors,
                 col=col,
                 border_color=NA,
                 fontsize=6,
                 treeheight_row=20,
                 treeheight_col=20,
		 annotation_names_col=FALSE)
}

## DEhelper.corr.pairwise: 
#' Heatmap of sample to sample distances (in a pairwise comparison).
#'
#' @param i - iterator to select contrast from conts [default = 1]
#' @param pairwise.dds - DESeq2 object
#' @return A plot (base plotting) object.
#'
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           DEhelper.DESeq2.corr.pairwise()
#'           
DEhelper.DESeq2.corr.pairwise <- function(i=1) {

        # extract assay and over-write gene_id in rownames with gene_name
        assay.rld <- assay(rlog(pairwise.dds[[i]]))
        rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(pairwise.dds[[i]])), gtf$gene_id)]

        # extract sample to sample distances
        sampleDists <- dist(t(assay.rld))
        sampleDistMatrix <- as.matrix(sampleDists)

        # set color scheme
	col <- plasma(255)
        
        # extract information for legend
        if (length(add_factors)==0) {
                legend.df <- data.frame(group=colData(pairwise.dds[[i]])[,c("group")],row.names=rownames(colData(pairwise.dds[[i]])))
        } else {
                legend.df <- as.data.frame(colData(pairwise.dds[[i]])[,c("group",add_factors)])
        }

        # fix group colors for legend and possible first additional factor if available
        if (length(add_factors)==0) {
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(pairwise.dds[[i]])[,"group"]))],
                                                       unique(colData(pairwise.dds[[1]])[,"group"])))
        } else {
                if (length(unique(colData(pairwise.dds[[i]])[,add_factors[1]])) <= 8) {
                        mypalette <- brewer.pal(8,"Dark2")[1:length(unique(colData(pairwise.dds[[i]])[,add_factors[1]]))]
                } else {
                        mypalette <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(colData(pairwise.dds[[i]])[,add_factors[1]])))
                }
                legend_colors <- list(group = setNames(brewer.pal(9,"Set1")[1:length(unique(colData(pairwise.dds[[i]])[,"group"]))],
                                                       unique(colData(pairwise.dds[[i]])[,"group"])),
                                      subject = setNames(mypalette,
                                                         unique(colData(pairwise.dds[[i]])[,add_factors[1]])))
                names(legend_colors) <- c("group",add_factors[1])
        }

	# sample to sample distance heatmap w/ 0's on diagonal set to NA
	# since they otherwise shift the scale too much
	sampleDistMatrix.diagNA <- sampleDistMatrix
	sampleDistMatrix.diagNA[sampleDistMatrix.diagNA==0] <- NA

        # plot heatmap
        pheatmap(sampleDistMatrix.diagNA,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 annotation_col=legend.df,
                 annotation_colors=legend_colors,
                 col=col,
                 border_color=NA,
                 fontsize=6,
                 treeheight_row=20,
                 treeheight_col=20,
		 annotation_names_col=FALSE)
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
     plotMA(res[[i]], main=conts[i, 1], alpha=fdr)
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
            scale_size_continuous("mean norm. counts (log10)")

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
## edgeR DE analysis
##
## DEhelper.init: some time consuming tasks that can be done in advance
DEhelper.edgeR.init <- function(task) {
    
    # Prepare the DE data frame
    renderUcscGeneLinks <- function() {
        ucsc_url <- paste0("http://genome.ucsc.edu/cgi-bin/hgGene?org=", SHINYREPS_ORG, "&db=", SHINYREPS_DB, "&hgg_gene=")
        for(i in 1:length(lrt)) {
            lrt[[i]]$table$gene <<- sapply(rownames(lrt[[i]]$table), function(x) {
                paste0("<a href=\"", ucsc_url, x, "\">", x, "</a>")
            })
        }
    }
    prepareDEdataTable <- function() {
        for(i in 1:length(lrt)) {
            #lrt[[i]]$table$FDR    <<- p.adjust(lrt[[i]]$table$PValue, method="fdr")
            lrt[[i]]$table$logFC  <<- round(lrt[[i]]$table$logFC, 2)
            lrt[[i]]$table$logCPM <<- round(lrt[[i]]$table$logCPM, 2)
            lrt[[i]]$table$LR     <<- round(lrt[[i]]$table$LR, 2)
            lrt[[i]]$table$PValue <<- lrt[[i]]$table$PValue
            lrt[[i]]$table$FDR    <<- lrt[[i]]$table$FDR
            lrt[[i]]$table$gene_name <<- lrt[[i]]$table$gene_name 
        }
    }
    
    # Cluster and correlation tasks
    prepareDistanceMatrix <- function() {
        v <<- apply(m, 1, sd, na.rm=T)        # get top variant genes
        dists <<- dist(t(m))
        mat <<- as.matrix(dists)
        hmcol <<- colorRampPalette(brewer.pal(max(length(levels(group)), 3), "Oranges"))(100)
    }
    
    # dispatch tasks
    switch(task, 
           renderUcscGeneLinks=renderUcscGeneLinks(), 
           prepareDEdataTable=prepareDEdataTable(), 
           prepareDistanceMatrix=prepareDistanceMatrix())
}

## DEhelper.MDS
DEhelper.edgeR.MDS <- function() {
    edgeR::plotMDS.DGEList(y, col=brewer.pal(max(length(levels(group)), 3), "Accent")[group])
}

##
## DEhelper.var: variance along log gene count-per-milion
##
DEhelper.edgeR.var <- function() {
    edgeR::plotBCV(y)
}

## DEhelper.cluster: Heatmap of top variant 'n' genes of the counts-per-milion table
DEhelper.edgeR.cluster <- function(n=50) {
    heatmap.2(m[rev(order(v))[1:n], ], col=hmcol, trace="none", margin=c(10, 6))
}

## DEhelper.corr: Heatmap of sample to sample distances
DEhelper.edgeR.corr <- function() {
    heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
}

## DEhelper.MAplot: MA plots
DEhelper.edgeR.MAplot <- function(i=1, fdr=.05) {
    # get DE genes (p.adjust='BH', pval<.05)
    de <- decideTestsDGE(lrt[[i]], p.value=fdr)
    degenes <- rownames(y)[as.logical(de)]
    
    # MA plot
    plotSmear(lrt[[i]], de.tags=degenes, main=names(lrt)[i])    # MA plot
    abline(h=c(-1, 1), col="blue")    # indicate 2-fold changes in the MA plot
    abline(v=0, col="blue")            # indicate >1 counts-per-million
}

## DEhelper.DEgenes: show the DE results
DEhelper.edgeR.DEgenes <- function(i=1) {
    ord  <- order(-log(lrt[[i]]$table$FDR), 
                   abs(lrt[[i]]$table$logFC), 
                  decreasing=TRUE)
#    cols <- c("gene", "logFC", "logCPM", "LR", "PValue", "FDR")
    cols <- c("gene_name", "logFC", "logCPM", "PValue", "FDR")
    lrt[[i]]$table[ord, cols]
}

##
## DEhelper.STAR: parse STAR log files and create a md table
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
            x[2] <- if(length(x) < 2) "" else x[2]
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
        kable(df[, 1, drop=F], align=c("r"), output=F)
    } else {
        kable(df, align=c("r", "r"), output=F)
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
                 "Number of reads mapped to too many loci",
                 "% of reads mapped to multiple loci",
                 "% of reads mapped to too many loci",
                 "% of reads unmapped: too many mismatches",
                 "% of reads unmapped: too short",
                 "% of reads unmapped: other"), function(x) {
                     as.numeric(gsub("%", "", gsub(".+\\|\t(.+)", "\\1", l[grep(x, l)])))
                 })    
    })
    
    # set row and column names, and output the md table
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
    df_values <- as.data.frame(t(x[2:7,]))
    df_values["Unmapped reads number"] <- x[1, ] - x[2, ] - x[4, ] - x[5,]
    df_values["% of reads unmapped"] <- x[8, ] + x[9, ] + x[10, ]
    df_values$sample <- rownames(df_values)
    # we clean up the colnames a little to make them shorter and nicer
    colnames(df_values) <- gsub("of reads mapped to", "",
                                gsub(" reads number", "", 
                                     gsub("Number of reads mapped to ", "",
                                          colnames(df_values))))
    #if we have a differential expression analysis
    #we refactor the samples depending on group/subject or alternatively on the
    #amount of unique_mapping reads
    if(file.exists(SHINYREPS_TARGET)){
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub(paste0(SHINYREPS_RNATYPES_SUFFIX,"$"), "",targets$file )
        add_factors <- colnames(targets)[!colnames(targets) %in% c("group", "sample", "file")]
        #we take the additional colnames present in the targets file and sort according
        #to them the factor levels
        targets <- targets[match(df_values$sample, targets$sample_ext),]
        #sort the samples according to the group
        df_values$sample <- targets$sample
        df_values$sample <- factor(df_values$sample, 
                                   levels=df_values$sample[ order(do.call( paste0,
                                                                     targets[,c("group", add_factors)]
                                                                    ))
                                                          ])
    } else{
      #we rorder according to the % amount of unique mapped reads mapped reads
      df_values$sample <- fct_reorder(df_values$sample, df_values$`Uniquely mapped reads %`)
    }
    df_values$sample <- gsub(lcSuffix(df_values$sample), "", df_values$sample)
    df_values$sample <- gsub(lcPrefix(df_values$sample), "", df_values$sample)
    df_melt <- melt(df_values, value.name = "reads", variable.name = "mapping_stat")
    df_melt$value_info <- ifelse(grepl("%", df_melt$mapping_stat), "perc", "reads")
    
    #we create two plots one for the % and one for the amount of reads in numbers
    p_perc <- ggplot(df_melt[df_melt$value_info == "perc",],
                     aes(x = sample, y = reads, fill = mapping_stat )) +
              geom_bar(stat     = "identity",
                       position = "stack") +
           ylab("% of reads sequenced") +
           labs(fill = "Mapping Statistic") +
	         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
	               plot.title = element_text(hjust=0.5)) +
           scale_fill_brewer(palette="Dark2") +
           ggtitle("Percentage of sequenced reads") 
    p_count <- p_perc %+%
               df_melt[df_melt$value_info == "reads",] +
               ylab("# reads") +
               ggtitle("Number of reads sequenced") 
      
    rownames(df_values) <- df_values$sample
    df_values <- df_values[, colnames(df_values) != "sample"]
    #we reformat individual columns
    df_values[, grepl("%", colnames(df_values))] <- as.data.frame(
      lapply(
        df_values[, grepl("%", colnames(df_values))], function(x){
          paste(format(x, nsmall=2), "%") 
        }))
    df_values[, !grepl("%", colnames(df_values))] <- as.data.frame(
      lapply(
        df_values[, !grepl("%", colnames(df_values))], function(x){
           format(x, big.mark=",") 
        }))
                     
    return( list(p_perc = p_perc,
                 p_count = p_count,
         stat = kable(df_values, align=c("r", "r", "r", "r","r"), output=F))
    )
    
}

##
## DEhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
DEhelper.Fastqc <- function(web=TRUE) {
    
    # logs folder
    if(!file.exists(SHINYREPS_FASTQC_OUT)) {
        return("Fastqc statistics not available")
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/fastqc" else SHINYREPS_FASTQC_OUT
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.dirs(SHINYREPS_FASTQC_OUT, recursive=F)
    df <- sapply(samples, function(f) {
        c(paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"),
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_sequence_gc_content.png)"))
    })

    # set row and column names, and output the md table
    df <- as.data.frame(t(df))
    rownames(df) <- gsub(paste0("^", SHINYREPS_PREFIX), "", basename(samples))
    rownames(df) <- gsub(paste0("_fastqc$"), "", rownames(df))
    colnames(df) <- c("Read qualities", "Sequence bias", "GC content")
    kable(df, output=F, align="c")
}

##
## DEhelper.dupRadar: go through dupRadar output dir and create a md table with
##     the duplication plots
##
DEhelper.dupRadar <- function(web=TRUE) {
    
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
    
    # put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        gsub("_dupRadar.png)", "", x)
    })
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
    
    # create a matrix using feature names as rownames, sample names as colnames
    l.infiles <- list.files(FOLDER, pattern=SUFFIX, full.names = T)
    l.counts <- lapply(l.infiles, function(f) {
        # f <- list.files(FOLDER, pattern=SUFFIX, full.names = T)[1]
        samplename <- basename(f)
        samplename <- gsub(SUFFIX, "", samplename)
        samplename <- gsub(paste0("^", SHINYREPS_PREFIX), "", samplename)
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
    
    plot <- ggplot() + 
        geom_bar(data=df.counts.melt, aes(x=sample, y=count, fill=type), position="fill", stat="identity") + 
        labs(x="", y="", fill="") +
	theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) 
    
    return(plot)
}


##
## DEhelper.geneBodyCov: go through geneBodyCov output dir and create a md table with
##     the gene coverage plots
##
DEhelper.geneBodyCov <- function(web=TRUE) {
    
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
    
    # put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        gsub("_geneBodyCov.png)", "", x)
    })
    df      <- matrix(df     , ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    samples <- matrix(samples, ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    colnames(df.names) <- rep(" ", SHINYREPS_PLOTS_COLUMN)
    
    kable(as.data.frame(df.names), align="c", output=F, format="markdown")
}

##
##DEhelper.strandspecifity: get the strandspecifity from the qc and display them
##
DEhelper.strandspecificity <- function(){

    # logs folder
    if(!file.exists(SHINYREPS_INFEREXPERIMENT_LOGS)) {
        return("Strand specificity statistics not available")
    }
    
    filelist <- list.files(path=SHINYREPS_INFEREXPERIMENT_LOGS, full.names=TRUE)
    strandspecifity <- lapply(filelist, read.table, sep=":", skip=3, header=FALSE, row.names=1, blank.lines.skip=TRUE)
    strandspecifity <- do.call(cbind, strandspecifity)
    samplenames <- basename(filelist)
    if(!is.na(SHINYREPS_PREFIX)) { 
	samplenames <- gsub(SHINYREPS_PREFIX, "", samplenames)
    }
    samplenames <- gsub("_inferexperiment.txt", "", samplenames)
    colnames(strandspecifity) <- samplenames 
    rownames(strandspecifity) <- c("other", "sense", "antisense") 
    kable(t(strandspecifity), output=F, align=c("l"))
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
## DEhelper.Bustard: call the perl XML interpreter and get the MD output
##
DEhelper.Bustard <- function() {
    f  <- SHINYREPS_BUSTARD
    
    if(!file.exists(f)) {
        return("Bustard statistics not available")
    }
    
    # call the perl XSL inetrpreter
    cmd <- paste(" bustard.pl", f)
    try(ret <- system2("perl", cmd, stdout=TRUE, stderr=FALSE))
    
    # check RC
    if(!is.null(attributes(ret))) {
        return(paste("Error parsing bustard statistics. RC:", attributes(ret)$status, "in command: perl", cmd))
    }
    
    ret     # ret contains already MD code
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
        
        
        sapply(c("Assigned",
                 "Unassigned_Ambiguity",
                 "Unassigned_MultiMapping",
                 "Unassigned_NoFeatures",
                 "Unassigned_Unmapped",
                 "Unassigned_MappingQuality",
                 "Unassigned_FragmentLength",
                 "Unassigned_Chimera",
                 "Unassigned_Secondary",
                 "Unassigned_Nonjunction",
                 "Unassigned_Duplicate"), function(y) {
                    as.numeric(  gsub( ".+\t(.+)", "\\1", l[grep(y, l)] )  )
                 })    
        
    })
    
    # correct column names
    colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
    
    # create md table (omitting various values that are 0 for now)
    #from x we romeove the ones which are unmapped to calculate percentages
    #only for the mapped ones
    x <- x[rownames(x) != "Unassigned_Unmapped", ]
    x <- rbind(total=x, colSums(x))
    rownames(x)[nrow(x)] <- "total"
    df <- data.frame(assigned=paste0(format(x[1, ], big.mark=","), " (", format((x[1, ]/x["total", ])*100, digits=2, nsmall=2), "%)"), 
                     unassigned_ambiguous=paste0(format(x[2, ], big.mark=","), " (", format((x[2, ]/x["total", ])*100, digits=2, nsmall=2), "%)"), 
                     unassigned_multimap=paste0(format(x[3, ], big.mark=","), " (", format((x[3, ]/x["total", ])*100, digits=2, nsmall=2), "%)"), 
                     unassigned_nofeature=paste0(format(x[4, ], big.mark=","), " (", format((x[4, ]/x["total", ])*100, digits=2, nsmall=2), "%)"))
    rownames(df) <- colnames(x)
    kable(df, align=c("r", "r", "r", "r"), output=F)
    
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
    samples <- list.files(QC, pattern="Reads.*.png$", recursive=T, full.names=T)
    if(length(samples) == 0) {
        return("Qualimap report not available")
    }
    df <- sapply(samples, function(f) {
        paste0("![alt text](", f, ")")
    })
    
    samples <- list.files(QC, pattern="Reads.*.png$", recursive=T, full.names=F)
    # put sample names and output an md table of 4 columns
    while(length(df) %% 2 != 0) df <- c(df, "")
    samples <-gsub(paste0("^", SHINYREPS_PREFIX), "", gsub("/.*", "", dirname(samples)))
    samples <- gsub("_counts_qualimap", "", samples)
    while(length(samples) %% 2 != 0) samples <- c(samples, "")
    df      <- matrix(df     , ncol=2, byrow=T)
    samples <- matrix(samples, ncol=2, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), ncol=2, byrow=T)
    colnames(df.names) <- c(" ", " ")
    
    kable(as.data.frame(df.names), output=F, format="markdown")
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
		samplenames <- gsub(SHINYREPS_PREFIX, "", samplenames)
		samplenames <- gsub("_insertsizemetrics.tsv","", samplenames)
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
  
  title_info <- gsub("_insertsizemetrics.tsv$", "", gsub(SHINYREPS_PREFIX, "", basename(metricsFile)))
  
  #we get the histogram which ahs the names of the leves e.g. all_reads and readgroups/sample groups depending on
  #accumulation level which was used.
  #the colnames of histogram are something like all.read.fr_count, all.read.rf_count etc.
  #to get the whole shebang into a wider format we have to add the information
  hist_long <- melt(histogram, id.var = "insert_size") %>% 
    separate( variable, 
              sep = "\\.",
              into = c("group", "counttype" )) %>%
    rename( amount = value)
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
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    tryCatch({
        ver <- read.delim(file=SHINYREPS_TOOL_VERSIONS)
        colnames(ver) <- c("Tool name","Environment", "Version")
        kable(as.data.frame(ver),output=F)
    }, error=function(e) cat("tool versions not available.\n", fill=TRUE))
}
