##################################
##
## helper functions to create the plots for the Shiny report
##
##################################


#' attach_packages: load required packages
#'
#' @param pkg character vector with package names
#'
#' @return character vector with package names which have not yet been attached to workspace before yet. 
#' This vector may be used to detach packages not needed any more.

attach_packages <- function(pkg) {
  not.yet.attached.pkg <- pkg[!(pkg %in% loadedNamespaces())]
  for (p in pkg) { 
    require(p, character.only = TRUE) 
  }
return(unique(not.yet.attached.pkg))
}



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
        #kable(df[, 1, drop=F], align=c("r"), output=F)
      kable(df[, 1, drop=F]) %>% kable_styling()
    } else {
        #kable(df, align=c("r", "r"), output=F)
      kable(df) %>% kable_styling()
      }
}

##
## DEhelper.STAR: parse STAR log files and create a md table
##
DEhelper.STAR <- function(colorByFactor=NULL, targetsdf=targets, ...) {
    
    # log file
    LOG <- SHINYREPS_STAR_LOG
    SUFFIX <- paste0(SHINYREPS_STAR_SUFFIX, '$')
    if(!file.exists(LOG)) {
        return("STAR statistics not available")
    }
    
    # look for the lines containing the strings
    # and get the values associated with this strings
    x <- list.files(LOG, pattern=SUFFIX)
    # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
    x <- selectSampleSubset(x, ...)
    
    x <- sapply(x, function(f) {
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
    colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
    colnames(x) <- gsub(lcSuffix(colnames(x)), "", colnames(x)) # remove longest common suffix
    colnames(x) <- gsub(lcPrefix(colnames(x)), "", colnames(x)) # remove longest common prefix
    df <- data.frame(input_reads=format(x[1, ], big.mark=","), 
                     uniquely_mapped=paste0(format(x[2, ], big.mark=","), " (", format(x[3, ], nsmall=2), "%)"), 
                     multi_mapped=paste0(format(x[4, ], big.mark=","), " (", format(x[6, ], nsmall=2), "%)"), 
                     too_many_loci=paste0(format(x[7,],nsmall=2), "%"),
                     unmapped=paste0(format( x[8, ] + x[9, ] + x[10, ], nsmall=2), "%")
                     )
  
    df.stacked <- data.frame(filename = gsub("\\.R[12]\\.*$", "", rownames(df)),
                             input = x[1, ],
                             unique = x[2, ],
                             unique_perc = 100*(x[2, ]/x[1, ]),
                             multi_perc = 100*(x[4, ]/x[1, ]),
                             mapped_perc = 100*((x[4, ]+x[2, ]) / x[1, ]))
    
## prepare groupwise plots
  if(!is.null(colorByFactor) && nrow(df.stacked) == nrow(targetsdf)) {
    # we want to plot the input reads and the mapped and the multi mapped reads numbers into different plots with different axises separated by one feature
    # we want to plot the same thing as percentages of the total
    # we have to plot per feature and then rearrange
    # we add one plot for the color value where we plot the percentages and color them according to the amount of input reads
  
    targetsdf$samplemod <- gsub(lcSuffix(targetsdf$sample ), "", targetsdf$sample ) # shorten filename suffix
    if(!is.na(SHINYREPS_PREFIX)) {targetsdf$samplemod  <- gsub(SHINYREPS_PREFIX, "", targetsdf$samplemod)}
    targetsdf$samplemod <- gsub(lcPrefix(targetsdf$sample ), "", targetsdf$sample ) # shorten filename prefix
    
    index <- as.numeric(sapply(targetsdf$samplemod, function(x) grep(x, df.stacked$filename, ignore.case = T))) # grep for sample name in shortened file names
    if((nrow(df.stacked) != length(index)) || any(is.na(index))) {
      stop("\nThere seem to be ambiguous sample names in targets. Can't assign them uniquely to STAR logfile names")
    }
    
    if(any(!colorByFactor %in% colnames(targetsdf))) {
      if(all(!colorByFactor %in% colnames(targetsdf))) {
        cat("\nNone of the column names given in colorByFactor is available. Perhaps sample names are not part of fastq file names? Using filename instead.\n")
        colorByFactor <- "filename"
      } else { # one plot each element of colorByFactor
        cat("\n", colorByFactor[!colorByFactor %in% colnames(targetsdf)], "not available. Using", colorByFactor[colorByFactor %in% colnames(targetsdf)], "\n")
        colorByFactor <- colorByFactor[colorByFactor %in% colnames(targetsdf)]
      }
    }

    targetsdf$filename <- df.stacked$filename[index]
    df.stacked <- merge(df.stacked, targetsdf[,unique(c("filename", colorByFactor)), drop=F], by="filename")
    df.stacked <- df.stacked[order(rownames(df.stacked)),, drop=F]
    df.stacked <- df.stacked[,!apply(df.stacked,2, function(x) any(is.na(x))), drop=F] # remove NA columns from unsuccessful matching
    
  } else {
    # if colorByFactor == NULL or targets does not fit to number of files
    colorByFactor <- "filename"
  } # end groupwise plots
    
    
    # melt data frame for plotting
    df.melt  <- melt(df.stacked, id.vars=unique(c(colorByFactor, "filename")), variable.name="map_feature")
   
    map.feature.plots <- lapply(colorByFactor, function(color.value){
      p <- ggplot(df.melt[!grepl("(perc)|(unique)", df.melt$map_feature),],
                  aes_string("map_feature", "value", color=color.value)) +
        geom_quasirandom() +
        scale_color_brewer(type= "qual", palette=2)  + # FR changed palette: scale_color_brewer(palette="Paired")+
        facet_wrap(~map_feature, scales="free") +
        scale_y_log10() +
        xlab(NULL) + 
        ylab("# Reads")
      p.perc <- ggplot(df.melt[grepl("perc", df.melt$map_feature),],
                       aes_string("map_feature", "value", color=color.value)) +
        geom_quasirandom() +
        scale_color_brewer(type= "qual", palette=2)  + # FR changed palette: scale_color_brewer(palette="Paired")+
        facet_wrap(~map_feature, scales="free") +
        xlab(NULL) + 
        guides(color="none") +
        ylab("% of input reads")  
        # theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
      return(list(p.perc,p))
    })
    
    for(i in 1:length(map.feature.plots)){
      grid.arrange(grobs=map.feature.plots[[i]], layout_matrix= matrix(c(1,2), nrow=1)) 
    }
    
    DT::datatable(df, options = list(pageLength= 20))
    #return(kable(df, format="markdown", output=F) %>% kable_styling())
}

	       
	       
##
## DEhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
DEhelper.Fastqc <- function(web=TRUE, ...) {
    
    # logs folder
    if(!file.exists(SHINYREPS_FASTQC_OUT)) {
        return("Fastqc statistics not available")
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/fastqc" else SHINYREPS_FASTQC_OUT
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.dirs(SHINYREPS_FASTQC_OUT, recursive=F)
    
    # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
    samples <- selectSampleSubset(samples, ...)
    
    df <- sapply(samples, function(f) {
        c(paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"),
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_sequence_gc_content.png)"))
    })

    # set row and column names, and output the md table
    df <- as.data.frame(t(df))
    rownames(df) <- gsub(paste0("^", SHINYREPS_PREFIX), "", basename(samples))
    # rownames(df) <- gsub(lcPrefix(rownames(df)), "", rownames(df)) # remove longest common prefix
    rownames(df) <- gsub(lcSuffix(rownames(df)), "", rownames(df)) # remove longest common suffix
    colnames(df) <- c("Read qualities", "Sequence bias", "GC content")
    kable(df, output=F, align="c")
}

##
## DEhelper.dupRadar: go through dupRadar output dir and create a md table with
##     the duplication plots
##
DEhelper.dupRadar <- function(web=TRUE, samplePattern=NULL, exclude=F) {
    
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

    # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
    samples <- selectSampleSubset(samples, samplePattern, exclude)
    
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
DEhelper.RNAtypes <- function(...) {
    
    FOLDER <- SHINYREPS_RNATYPES
    SUFFIX <- paste0(SHINYREPS_RNATYPES_SUFFIX, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return("Subread statistics not available")
    }
    
    # create a matrix using feature names as rownames, sample names as colnames
    l.infiles <- list.files(FOLDER, pattern=SUFFIX, full.names = T)
    
    # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
    l.infiles <- selectSampleSubset(l.infiles, ...)
    
    l.counts <- lapply(l.infiles, function(f) {
        # f <- list.files(FOLDER, pattern=SUFFIX, full.names = T)[1]
        samplename <- basename(f)
        samplename <- gsub(paste0("^", SHINYREPS_PREFIX), "", samplename)
        samplename <- gsub(lcSuffix(basename(l.infiles)), "", samplename) # remove longest common suffix
        df.readcount <- read.table(f, header=T, sep='\t', comment.char = '#', col.names=c("Geneid","Length",samplename))
        df.readcount$Length <- NULL
        return(df.readcount)
    })
    
    # merge counts to data frame in wide format & remove "length..." columns
    f.merge <- function(x,y) {merge(x,y,by="Geneid")}
    df.counts <- Reduce(f.merge, l.counts)
    df.counts$Geneid[df.counts$Geneid==""] <- "not annotated"
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
DEhelper.geneBodyCov <- function(web=TRUE, ...) {
    
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
    
    # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
    samples <- selectSampleSubset(samples, ...)
    
    df <- sapply(samples, function(f) {
        paste0("![geneBodyCov img](", QC, "/", basename(f), ")")
    })
    
    # put sample names and output in md table of SHINYREPS_PLOTS_COLUMN columns
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        gsub(paste0(lcSuffix(samples),")"), "", x)
        # gsub("_geneBodyCov.png)", "", x) # remove longest common suffix
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
##DEhelper.strandspecificity: get the strandspecificity from the qc and display them
##
DEhelper.strandspecificity <- function(samplePattern=NULL, ...){

    # logs folder
    if(!file.exists(SHINYREPS_INFEREXPERIMENT_LOGS)) {
        return("Strand specificity statistics not available")
    }
    
    filelist <- list.files(path=SHINYREPS_INFEREXPERIMENT_LOGS, full.names=TRUE)
    # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
    filelist <- selectSampleSubset(filelist, samplePattern, ...)
    
    strandspecifity <- lapply(filelist, read.table, sep=":", skip=3, header=FALSE, row.names=1, blank.lines.skip=TRUE)
    strandspecifity <- do.call(cbind, strandspecifity)
    samplenames <- basename(filelist)
    if(!is.na(SHINYREPS_PREFIX)) { 
	samplenames <- gsub(SHINYREPS_PREFIX, "", samplenames)
    }
    samplenames <- gsub(lcSuffix(samplenames), "", samplenames)
    colnames(strandspecifity) <- samplenames 
    rownames(strandspecifity) <- c("other", "sense", "antisense") 
    kable(t(strandspecifity), output=F, align=c("l")) %>% kable_styling()
}





##
##DEhelper.cutadapt: display read trimming stats from cutadapt
##
DEhelper.cutadapt <- function(colorByFactor=NULL, targetsdf=targets, ...){
  
# x <- sapply(list.files(SHINYREPS_CUTADAPT_LOGS,pattern='imb_gcf.*.log$',full.names=TRUE), function(f) { 
x <- list.files(SHINYREPS_CUTADAPT_LOGS,pattern='*cutadapt.log$',full.names=TRUE) 
# select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL

x <- selectSampleSubset(x, ...)

x <- sapply(x, function(f) { 
  total.reads <- system(paste("grep \"Total reads processed\"", f, "| awk '{print $4}'"), intern=TRUE)
  total.reads <- gsub(",", "", total.reads)
  
  trimmed.reads.perc <- system(paste("grep \"Reads with adapters\"", f, "| awk '{print $5}'"), intern=TRUE)
  trimmed.reads.perc <- gsub("\\(|\\)|\\%", "", trimmed.reads.perc)
  
  tooshort.reads.perc <- system(paste("grep \"Reads that were too short\"", f, "| awk '{print $7}'"), intern=TRUE)
  tooshort.reads.perc <- gsub("\\(|\\)|\\%", "", tooshort.reads.perc)
  
  # trimming of each adapter
  adapters <- system(paste("grep Sequence:", f, "| awk '{print $9}'"), intern=T)
  adapters.perc <- round(100*(as.numeric(adapters) / as.numeric(total.reads)),2)
  names(adapters.perc) <- gsub(" *=== *", "", system(paste("grep \"=== Adapter\"", f), intern=T))
  
  ## add trimmed reads for each adapter here
  return(c(total.reads, trimmed.reads.perc, tooshort.reads.perc, adapters.perc))
})

# set row and column names
x.df <- as.data.frame(t(x)) 
colnames(x.df)[1:3] <- c("total.reads", "trimmed","tooshort")
x.df <- as.data.frame(lapply(x.df, as.numeric))

#reduce size of file names 
row.names(x.df) <- basename(colnames(x))
row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
if(!is.na(SHINYREPS_PREFIX)) {row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))}
x.df$filename <- factor(row.names(x.df))


# passing the different factors given in targetsdf to x.df which was created from cutadapt logfile names (if 1 cell per file)
if(!is.null(colorByFactor) && nrow(x.df) == nrow(targetsdf)) { # if targets object fits in length, add information to x.df
  

  targetsdf$samplemod <- gsub(lcSuffix(targetsdf$sample ), "", targetsdf$sample ) # shorten filename suffix
  if(!is.na(SHINYREPS_PREFIX)) {targetsdf$samplemod  <- gsub(SHINYREPS_PREFIX, "", targetsdf$samplemod)}
  targetsdf$samplemod <- gsub(lcPrefix(targetsdf$sample ), "", targetsdf$sample ) # shorten filename prefix
  
  
  #index <- as.numeric(sapply(x.df$filename, function(x) grep(x, targetsdf$samplemod, ignore.case = T))) # grep for shortened file names in sample names
  # x.df <- cbind(x.df, targetsdf[index, , drop=F ]) 
  index <- as.numeric(sapply(targetsdf$samplemod, function(x) grep(x, x.df$filename, ignore.case = T))) # grep for sample name in shortened file names
      if(nrow(x.df) != length(index) || any(is.na(index))) {
        stop("\nThere seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
      }
  
  targetsdf$filename <- x.df$filename[index]
  x.df <- merge(x.df, targetsdf, by="filename")
  x.df <- x.df[order(rownames(x.df)),, drop=F]
  x.df <- x.df[,!apply(x.df,2, function(x) any(is.na(x))), drop=F] # remove NA columns from unsuccessful matching
  rownames(x.df) <- x.df$filename
  
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
    # if colorByFactor == NULL or targets does not fit to number of files
    colorByFactor <- "filename"
  }

# melt data frame for plotting
x.melt <- melt(x.df, measure.vars=c("trimmed", "tooshort", grep("Adapter", colnames(x.df), value=T)),
               #id.vars=group.vars,
               variable="reads")
# everything which is not a value should be a factor

# now we do a violin plot of the trimmed/too_short/etc. ones and color it
# according to the different factors given in colorByFactor 
create.violin <- function(x.melt, color.value){
  ylab <- "% reads"
  p <- ggplot(x.melt, aes_string(x="reads",
                                 y="value",
                                 color=color.value ))+
    geom_quasirandom() +
    scale_color_brewer(type= "qual", palette=2)  +    # replaced color palette: scale_color_brewer(palette="Paired")
    ylab(ylab) +
    xlab("") +
    scale_y_continuous( breaks=seq(0, max(x.melt$value), 10),
                        limits = c(0, max(x.melt$value))) + 
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) 
  
  return(p)
}

# one plot for each element of colorByFactor
violin.list <- lapply(colorByFactor, create.violin, x.melt=x.melt) # "colorByFactor" is submitted as color.value

for(i in 1:length(violin.list)){
  plot(violin.list[[i]])
  }

#kable(x.df[,c("total.reads", "trimmed","tooshort")], output=F, format="markdown", align=c("l")) %>% kable_styling()
DT::datatable(x.df[,c("total.reads", "trimmed","tooshort", grep("Adapter", colnames(x.df), value=T))], options = list(pageLength= 20))
}


##
##DEhelper.umicount: display deduplication stats from UMI_tools count
## 
DEhelper.umicount <- function(colorByFactor=NULL, targetsdf=targets, ...){
  

x <- list.files(SHINYREPS_UMICOUNT_LOG,pattern='*umicount.log$',full.names=TRUE) 
# select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL

x <- selectSampleSubset(x, ...)

x <- sapply(x, function(f) { 
  input_reads_total <- system(paste("grep \"INFO Input Reads\"", f, "| awk '{print $6}'"), intern=TRUE)
  
  skipped_reads_total <- system(paste("grep \"INFO Read skipped, no tag\"", f, "| awk '{print $8}'"), intern=TRUE)

  counted_reads_total <- system(paste("grep \"INFO Number of (post deduplication) reads counted\"", f, "| awk '{print $10}'"), intern=TRUE)

  return(c(unique(input_reads_total), unique(skipped_reads_total), unique(counted_reads_total)))
})

# set row and column names
x.df <- as.data.frame(t(x)) 
colnames(x.df) <- c("input_reads_total", "skipped_reads_total","counted_reads_total")
x.df <- as.data.frame(lapply(x.df, as.numeric))
x.df$skipped_reads <- round(100* (x.df$skipped_reads_total / x.df$input_reads_total ), 2)
x.df$counted_reads <- round(100* (x.df$counted_reads_total / x.df$input_reads_total), 2)


#reduce size of file names 
row.names(x.df) <- basename(colnames(x))
row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
if(!is.na(SHINYREPS_PREFIX)) {row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))}
x.df$filename <- factor(row.names(x.df))


# passing the different factors given in targetsdf to x.df which was created from cutadapt logfile names (if 1 cell per file)
if(!is.null(colorByFactor) && nrow(x.df) == nrow(targetsdf)) { # if targets object fits in length, add information to x.df
  
  
  targetsdf$samplemod <- gsub(lcSuffix(targetsdf$sample ), "", targetsdf$sample ) # shorten filename suffix
  if(!is.na(SHINYREPS_PREFIX)) {targetsdf$samplemod  <- gsub(SHINYREPS_PREFIX, "", targetsdf$samplemod)}
  targetsdf$samplemod <- gsub(lcPrefix(targetsdf$sample ), "", targetsdf$sample ) # shorten filename prefix
  
  
  #index <- as.numeric(sapply(x.df$filename, function(x) grep(x, targetsdf$samplemod, ignore.case = T))) # grep for shortened file names in sample names
  # x.df <- cbind(x.df, targetsdf[index, , drop=F ]) 
  index <- as.numeric(sapply(targetsdf$samplemod, function(x) grep(x, x.df$filename, ignore.case = T))) # grep for sample name in shortened file names
  if(nrow(x.df) != length(index) || any(is.na(index))) {
    stop("\nThere seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
  }
  
  targetsdf$filename <- x.df$filename[index]
  x.df <- merge(x.df, targetsdf, by="filename")
  x.df <- x.df[order(rownames(x.df)),, drop=F]
  x.df <- x.df[,!apply(x.df,2, function(x) any(is.na(x))), drop=F] # remove NA columns from unsuccessful matching
  rownames(x.df) <- x.df$filename

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
  # if colorByFactor == NULL or targets does not fit to number of files
  colorByFactor <- "filename"
}

# melt data frame for plotting
x.melt <- melt(x.df, measure.vars=c("skipped_reads", "counted_reads"),
               variable="reads")
#everything which is not a value should be a factor

#now we do a violin plot of the trimmed/too_short/etc. ones and color it
# according to the different factors given in colorByFactor 

create.violin <- function(x.melt, color.value){
  ylab <- "% reads"
  p <- ggplot(x.melt, aes_string(x="reads",
                                 y="value",
                                 color=color.value ))+
    geom_quasirandom() +
    scale_color_brewer(type= "qual", palette=2)  +    # FR replaced color palette: scale_color_brewer(palette="Paired")
    ylab(ylab) +
    xlab("") +
    scale_y_continuous( breaks=seq(0, max(x.melt$value), 10),
                        limits = c(0, max(x.melt$value))) 
  return(p)
}

# one plot each element of colorByFactor
violin.list <- lapply(colorByFactor, create.violin, x.melt=x.melt) # "colorByFactor" submitted as color.value

for(i in 1:length(violin.list)){
  plot(violin.list[[i]])
}

#kable(x.df[,c("total.reads", "trimmed","tooshort")], output=F, format="markdown", align=c("l")) %>% kable_styling()
colnames(x.df)[colnames(x.df)=="skipped_reads"] <- "skipped_reads_perc"
colnames(x.df)[colnames(x.df)=="counted_reads"] <- "counted_reads_perc"
DT::datatable(x.df[,c("input_reads_total", "skipped_reads_total","skipped_reads_perc", "counted_reads_total", "counted_reads_perc")])
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
    # kable(df, align=c("r", "r", "r", "r"), output=F)
    kable(df) %>% kable_styling()
    
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
    toolList <- SHINYREPS_TOOL_VERSIONS
    ver <- read.table(file=toolList,sep="=")
    ver$V1 <- strsplit(as.character(ver$V1),"_VERSION")
    colnames(ver) <- c("Tool name","Version")

    kable(as.data.frame(ver),output=F) %>% kable_styling()
}


##
#' isoutlier function using median absolute deviation (MAD) as defined in scran
#' 
MAD <- function(x, n) {
  constant <- 1.4826 # keep consistency with sd in normally distributed data (see ?stats::mad)
  list(low =median(x) - n * constant * median(abs(x - median(x))),
       high=median(x) + n * constant * median(abs(x - median(x))))
}


##
#' plotPCAfromQCmetrics: plot PCA from qc metrics of sce object
#'
#' @param sce SingleCellExperiment object
#' @param metrics character vector with QC metric to be used from colData(sce)
#' @param anno character vector defining up to 2 cell grouping factors to be indicated in the PCA plot by color and shape
#' 
plotPCAfromQCmetrics <- function(sce, metrics, anno){
  
  if(length(anno)>2) {stop("\nno more than 2 categories allowed in 'anno'.")}
  dotcol <- anno[1]
  if(length(anno)==2) { 
    dotshape <- anno[2]} else {
      dotshape <- NULL}
  
  pca.sce <- runPCA(sce, use_coldata=TRUE, selected_variables = metrics)  
  pca <- as.data.frame(pca.sce@reducedDims)
  rownames(pca) <- rownames(colData(pca.sce)) 
  pca <- cbind(pca, as.data.frame(colData(sce)[match(rownames(pca), colnames(sce)), unique(c("cells",anno)), drop=F]))
  
  p_allGroups <-  plotReducedDim(pca.sce, use_dimred="PCA_coldata", colour_by=dotcol, shape_by=dotshape, point_size=3) +
    geom_text_repel(aes(x=PC1, y=PC2, label=cells), subset(pca, cells %in% c("0c", "10c")))  
    # scale_color_manual(dotcol, 
    #                  values=c("grey50","#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3","#FFFF33","#A65628","#F781BF"), 
    #                  guide="legend")  # color scheme as in QC PCA plot
    
  plot(p_allGroups)
  }



##
#' selectSampleSubset: select subset of samples for including in report (e.g. in case of multiple fastq files in scRNA-seq) 
#'
#' @param samples character vector with sample names
#' @param samplePattern regular expression to apply on \code{samples}
#' @param exclude logical indicating if selected samples shall be excluded or included
#'
#' @return character vector with selected sample names

selectSampleSubset <- function(samples, samplePattern=NULL, exclude=F, maxno=NULL) {
    # use all samples for samplePattern=NULL
    if(!is.null(samplePattern)) {
      samplefilenames <- basename(samples) # apply pattern to filename, not to full path
      samples <- samples[grep(samplePattern, samplefilenames, invert=exclude)]
      if(length(samples)==0) {stop("\nYou have selected no files!\n")}
    }
    if(!is.null(maxno)) {
      samples <- samples[1:min(length(samples), maxno)]
      if(maxno > length(samples)) {cat("\nSample number restricted to", maxno)}
    }
    return(samples)
  }
  

