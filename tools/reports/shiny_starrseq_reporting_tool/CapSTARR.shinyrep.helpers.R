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
#'           # rld and dds are passed via function arguments
#'           DEhelper.DESeq2.MDS(dds=dds, rld=rld)
#' 
DEhelper.DESeq2.MDS <- function(dds=NULL, rld=NULL) {

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
#'           # pairwise.dds is passed via function arguments
#'           DEhelper.DESeq2.pairwisePCA(i, pairwise.dds=pairwise.dds)
#'
DEhelper.DESeq2.pairwisePCA <- function(i=1, pairwise.dds=NULL) {

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
#'           # dds is passed via function arguments, and should be DESeq2 rld object
#'           # gtf is also passed via function arguments
#'           DEhelper.DESeq2.heatmap(dds=rld, gtf=gtf)
#'           
DEhelper.DESeq2.heatmap <- function(i=NULL, dds=NULL, gtf=NULL, logTransform=F, anno_factors = c("group","replicate"), type="distance", n=40) {
  
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
#'           # res & conts are passed via function arguments
#'           DEhelper.DESeq2.MAplot(res=res, conts=conts)
#'           
DEhelper.DESeq2.MAplot <- function(i=1, fdr=.01, res=NULL, conts=NULL) {
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
#'           # res is passed via function arguments
#'           DEhelper.DESeq2.DEgenes(i=2, res=res)
#'           
DEhelper.DESeq2.DEgenes <- function(i=1, res=NULL) {
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
#'           # res is passed via function arguments
#'           DEhelper.DESeq2.VolcanoPlot(res=res, i = 1, fdr = 0.01, top = 20, web = FALSE)
#'           
DEhelper.DESeq2.VolcanoPlot <- function(res=NULL, i=1, fdr=.01, top=25, web=TRUE) {
    # gather results
    d <- as.data.frame(DEhelper.DESeq2.DEgenes(i, res = res))
    x.limit <- max(abs(d$log2FoldChange), na.rm=T) # find the maximum spread of the x-axis
    
    # plotting
    p <- ggplot(d) +
            geom_point(mapping=aes(log2FoldChange, -log10(padj), size=log10(baseMean+1), color=padj < fdr), alpha=0.1) +
            theme_bw() +
            xlim(-x.limit, x.limit) +
            ylab("-log10 adj. p-value") +
            xlab("log2 fold change") + 
            scale_color_manual(values=c("black", "red")) +
            scale_size_continuous("mean norm. counts (log10)") +
	          guides(color="none", size=guide_legend(nrow=1)) +
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
#'           # res is passed via function arguments
#'           table_fisher <- DEhelper.DESeq2.ChrOverrepresentation(res=res, i = 1, fdr_de_gene = 0.1, fdr_fisher_test = 0.1, filter = TRUE)
#'           
DEhelper.DESeq2.ChrOverrepresentation <- function(res=NULL, i=1, fdr_de_gene=0.1, fdr_fisher_test=0.1, filter=TRUE) {
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
## DEhelper.BOWTIE2: parse bowtie2 paired end or single read mapping log files and create a md table
##
DEhelper.Bowtie2 <- function() {
  
  # log file
  LOG <- SHINYREPS_BOWTIE_LOG
  if(!file.exists(LOG)) {
    return("Bowtie statistics not available")
  }
  
  # look for the lines containing the strings
  # and get the values associated with this strings
  x <- sapply(list.files(LOG), function(f) {
    
    x <- file(paste0(LOG, "/", f))
    l <- readLines(x)
    close(x)
    
    stat.lines <- c("were (un)*paired",
                    "aligned concordantly exactly 1 time",
                    "aligned concordantly >1 times",
                    "aligned discordantly 1 time",
                    "aligned exactly 1 time",
                    "aligned >1 times",
                    "overall alignment rate")
    
    stat.names <- c("# read (pairs)",
                    "unique",
                    "multi",
                    "discordantly",
                    "single unique",
                    "single multi",
                    "overall align. rate")
    
    if(SHINYREPS_PAIRED != "yes") { # modify for single read design
      remove_lines <- grep("cordantly", stat.lines)
      stat.lines <- stat.lines[-remove_lines]
      stat.names <- stat.names[-remove_lines]
      stat.names <- gsub("single ", "", stat.names)
    }
    
    stats <- sapply(l[sapply(stat.lines, grep, x=l)], 
                    function(x) { 
                      sub("^\\s+", "", gsub(" \\(.*", "", x))
                    })	
    #recount the percentages for the stats
    stat.all <- stats[grep("overall alignment rate", names(stats))]
    stat.all <- gsub("overall alignment rate", "", stat.all)
    stats <- as.numeric(stats[!grepl("overall alignment rate", stats)])
    # and add the duplicates information
    stats.percent <- paste0("(", format((stats/stats[1])*100, digits=2), "%)")
    stats <- paste0(format(stats, big.mark=","), " ", stats.percent)
    stats <- c(stats, stat.all)
    stat.lines <- gsub("\\$", "", stat.lines)
    
    names(stats) <- stat.names
    return(stats)    
  })
  
  # set row and column names, and output the md table
  colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
  colnames(x) <- gsub(".bam.log$", "", colnames(x))
  colnames(x) <- gsub("\\..*$", "", colnames(x))
  rownames(x)[1] <- if(SHINYREPS_PAIRED != "yes") {"all reads"} else {"all pairs"}
  kable(t(x), align=c(rep("r",10)), output=F, format="html", row.names=T) %>% kableExtra::kable_styling()
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
            guides(fill="none", color="legend") +
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
    lbls <- gsub("_fastqc.zip$", "", basename(names(fastqc.stats)))
    names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", basename(names(fastqc.stats)))

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
           guides(color="none",
                  fill=guide_legend(title="",ncol=3)) +
           theme(axis.text.x = element_text(size=6,angle=90,hjust=0.5,vjust=0.5),
                 axis.text.y = element_text(size=8),
                 axis.title  = element_text(size=10),
                 plot.title  = element_text(size=12),
                 legend.text = element_text(size=7),
                 legend.position = "top") +
           scale_x_continuous(breaks=unique(as.numeric(df$position)),
                              labels=unique(df$Base))

        p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=rev(lbls)) +
          labs(y = "") +
          theme(legend.position="right") +
          guides(fill="none", color="legend") +
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
    # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", theoreticalGC=TRUE, species=SPECIES)
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
## DEhelper.Subread: parse Subread summary stats and create a md table
##
DEhelper.Subread <- function(subdir="") {
    
    FOLDER <- file.path(SHINYREPS_SUBREAD, subdir)
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
##DEhelper.insertsize: get the insertsize from the qc and display mean and sd 
##
DEhelper.insertsize <- function(subdir=""){

	if (SHINYREPS_PAIRED == "yes") {
		filelist <- list.files(path=file.path(SHINYREPS_INSERTSIZE, subdir), full.names=TRUE, pattern="insertsizemetrics.tsv$")
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
DEhelper.insertsize.plot <- function(subdir=""){
  # logs folder
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),
                                     error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
  }
  if (SHINYREPS_PAIRED == "yes" &
      length(list.files(path = file.path(SHINYREPS_INSERTSIZE, subdir),
                        pattern = "insertsizemetrics.tsv$")) > 0) {
    samples <- list.files(path = file.path(SHINYREPS_INSERTSIZE, subdir),
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
    
    trimmed.R1.perc <- trimmed.R2.perc <- trimmed.reads.perc <- trimmed.qual.perc <- tooshort.reads.perc <- toolong.reads.perc <- NULL 
    
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
    
    trimmed.qual.perc <- system(paste("grep \"Quality-trimmed\"", f, "| awk '{print $4}'"), intern=TRUE)
    trimmed.qual.perc <- gsub("\\(|\\)|\\%", "", trimmed.qual.perc)
    tooshort.reads.perc <- system(paste("grep \"that were too short\"", f, "| awk '{print $7}'"), intern=TRUE)
    tooshort.reads.perc <- gsub("\\(|\\)|\\%", "", tooshort.reads.perc)
    toolong.reads.perc <- system(paste("grep \"that were too long\"", f, "| awk '{print $7}'"), intern=TRUE)
    toolong.reads.perc <- gsub("\\(|\\)|\\%", "", toolong.reads.perc)
    
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
    return(c("total reads"=total.reads, 
             "bp quality trimmed"=trimmed.qual.perc,
             "R1 adapter trimmed"=trimmed.R1.perc, "R2 adapter trimmed"=trimmed.R2.perc, "adapter trimmed"=trimmed.reads.perc, 
             "too short"=tooshort.reads.perc, "too long"=toolong.reads.perc, adapters.perc))
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
  vars2plot <- c(grep("adapter trimmed", colnames(x.df), value=T), 
                 "too short", "too long",
                 grep("(Adapter)|(})", colnames(x.df), value=T))
  vars2plot <- vars2plot[vars2plot %in% colnames(x.df)]
  x.melt <- reshape2::melt(x.df, measure.vars=vars2plot, variable.name="reads")
  # everything which is not a value should be a factor
  
  # one plot for each element of colorByFactor
  violin.list <- lapply(colorByFactor, plotfun, data=x.melt, labelOutliers=labelOutliers, outlierIQRfactor=outlierIQRfactor) # "colorByFactor" is submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  vars4table <- c(colorByFactor, "total reads", 
                  grep("adapter trimmed", colnames(x.df), value=T),
                  "bp quality trimmed", "too short", "too long",
                  grep("(Adapter)|(})", colnames(x.df), value=T))
  vars4table <- vars4table[vars4table %in% colnames(x.df)]
  vars4table.colnames <- vars4table
  vars4table.colnames[!vars4table.colnames %in% c(colorByFactor, "total reads")] <- paste("%", vars4table.colnames[!vars4table.colnames %in% c(colorByFactor, "total reads")])
  DT::datatable(x.df[,vars4table], 
                options = list(pageLength= 20),
                colnames=vars4table.colnames)
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
##ChIPhelper.insertsize: get the insertsize from the qc and display mean and sd
##For ChIP type data, such as normal STARR-seq or log fold change version of capture STARR-seq
##
ChIPhelper.insertsize <- function(subdir="", 
                                  sampleColumnName =c("IPname", "INPUTname"), 
                                  fileColumnName =c("IP", "INPUT"), ...){
  
  if (SHINYREPS_PAIRED == "yes") {
    filelist <- list.files(path=file.path(SHINYREPS_INSERTSIZE, subdir), full.names=TRUE, pattern="insertsizemetrics.tsv$")
    filelist <- selectSampleSubset(filelist, ...)
    insertsizes <- lapply(filelist, read.table, sep="\t", header=TRUE, nrow=1)
    insertsizes <- do.call(rbind, insertsizes)
    samplenames <- basename(filelist)
    
    if(file.exists(SHINYREPS_TARGET)){
      
      # get target names
      targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors = F)
      
      if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
        targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
        targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
        targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
        for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
        targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
      } else {
        targets$file <- targets[,fileColumnName]
        targets$sample <- targets[,sampleColumnName]
      }
      
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
    knitr::kable(insertsizes, output=F, align=c("l"), format="html") %>% kableExtra::kable_styling()
  }
}


##
## Helper to plot the insertsize histogram equivalent to the one from picard
## Input is the Picard generated metrics file
##
ChIPhelper.insertsize.helper <- function(metricsFile, 
                                         sampleColumnName =c("IPname", "INPUTname"), 
                                         fileColumnName =c("IP", "INPUT")){
  #find the start of our metrics information 
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
    targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors = F)
    
    if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
      targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
      targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
    } else {
      targets$file <- targets[,fileColumnName]
      targets$sample <- targets[,sampleColumnName]
    }
    
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
  
  #we get the histogram which ahs the names of the leves e.g. all_reads and readgroups/sample groups depending on
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
## ChIPhelper.subchunkify: small function stolen from here 
##http://michaeljw.com/blog/post/subchunkify/
## to dynamically create chunks and adjust their size accordingly.
##
ChIPhelper.subchunkify <- function(g, fig_height=7, fig_width=5) {
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
##ChIPhelper.insertsize.plot: get the insertsize histograms and display them 
##
ChIPhelper.insertsize.plot <- function(subdir="", ...){
  # logs folder
  
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2 | SHINYREPS_PLOTS_COLUMN > 3) {
    SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
  }
  
  if (SHINYREPS_PAIRED == "yes" &
      length(selectSampleSubset(list.files(path = file.path(SHINYREPS_INSERTSIZE, subdir),
                                           pattern = "insertsizemetrics.tsv$"), ...)) > 0) {
    samples <- list.files(path = file.path(SHINYREPS_INSERTSIZE, subdir),
                          full.names = TRUE,
                          pattern = "insertsizemetrics.tsv$")
    samples <- selectSampleSubset(samples, ...)
    
    #we generate the plots
    insert_plots <- lapply(samples, ChIPhelper.insertsize.helper)
    return(arrangeGrob(grobs = insert_plots,
                       ncol = SHINYREPS_PLOTS_COLUMN))
  }else{
    return("No insertsize histograms available.")
  }
}


#' ChIPhelper.Fastqc.custom: prepare customized Fastqc summary plots
#' 
#' @param summarizedPlots logical, if TRUE, data from all samples is summarized in a single plot.
#' @param subdir character with sub-directory to append to the target directory.
#' @param metrics character vector with FastQC plot types to be included. Any combination of "Summary", "BaseQuals", "SeqQuals", "SeqContent", "GcContent", "DupLevels", "Overrep", "AdapterContent".
#' @param sampleColumnName character vector with column names of targets file indicating sample names.
#' @param fileColumnName character vector with column names of targets file indicating sample file names (must have order corresponding to sampleColumnName).
#'
#' @return list of ggplots made from FastQC data
#' 
ChIPhelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, subdir="", 
                                     metrics=c("BaseQuals", "SeqContent", "GcContent", "DupLevels"),
                                     sampleColumnName =c("IPname", "INPUTname"), 
                                     fileColumnName =c("IP", "INPUT")) {
  
  # logs folder
  if(!file.exists(SHINYREPS_FASTQC_OUT)) {
    return("Fastqc statistics not available")
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
  # read fastqc results in the appropriate format
  f <- list.files(QC, pattern="\\.zip$",full.names=T)
  fastqc.stats <- ngsReports::FastqcDataList(f)
  
  qclist <- list() # initialize return object
  qclist[["no.of.samples"]] <- length(f)
  
  # create proper name vectoir as labels
  lbls <- gsub("_fastqc.zip$", "", basename(names(fastqc.stats)))
  names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", basename(names(fastqc.stats)))
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors = F)
    
    # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format
    if(length(fileColumnName)>1) {
      targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
      targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
    } else {
      targets$file <- targets[,fileColumnName]
      targets$sample <- targets[,sampleColumnName]
    }
    
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
         guides(color="none",
                fill=guide_legend(title="",ncol=3)) +
         theme(axis.text.x = element_text(size=6,angle=90,hjust=0.5,vjust=0.5),
               axis.text.y = element_text(size=8),
               axis.title  = element_text(size=10),
               plot.title  = element_text(size=12),
               legend.text = element_text(size=7),
               legend.position = "top") +
         scale_x_continuous(breaks=unique(as.numeric(df$position)),
                            labels=unique(df$Base))
      
      p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=rev(lbls)) +
        labs(y = "") +
        theme(legend.position="right") +
        guides(fill="none", color="legend") +
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
  # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", theoreticalGC=TRUE, species=SPECIES)
  # the default value for SPECIES is "Hsapiens", thus, if you don't specify it, human will be used as a default
  p.gc <- ngsReports::plotGcContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=FALSE) 
  if(!summarizedPlots) {
    p.gc <- p.gc + guides(color=guide_legend(title="",ncol=4)) + 
      theme(legend.position = "top", legend.text = element_text(size=8)) 
  }
  
  return(list(no.of.samples=length(f), p.qual=p.qual, p.content=p.content, p.gc=p.gc))
}


##
## ChIPhelper.fastqscreen: add FastqScreen data to and plot it as a barplot
##
#' ChIPhelper.fastqscreen: summarizes FastQScreen results, creates summarized barplots, only relevant contanimants shown
#'
#' @param perc.to.plot - a numeric vector of length 1 setting the percent cutoff of relevant contaminants, if any sample
#'                       shows more than perc.to.plot, contaminant will be shown in plot
#' @param sampleColumnName - character vector with column names of targets file indicating sample names
#' @param fileColumnName - character vector with column names of targets file indicating sample file names (must have order corresponding to sampleColumnName)
#'
#' @return a list including a plot, the number of samples, and the number of plotted contaminants
#'
ChIPhelper.fastqscreen <- function(perc.to.plot = 1,
                                   sampleColumnName =c("IPname", "INPUTname"), 
                                   fileColumnName =c("IP", "INPUT")) {
  
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
    
    # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format
    if(length(fileColumnName)>1) {
      targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
      targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
    } else {
      targets$file <- targets[,fileColumnName]
      targets$sample <- targets[,sampleColumnName]
    }
    
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
## ChIPhelper.dupRadar: go through dupRadar output dir and create a md table with
##     the duplication plots
##
ChIPhelper.dupRadar <- function(web=FALSE,
                                sampleColumnName =c("IPname", "INPUTname"), 
                                fileColumnName =c("IP", "INPUT")) {
    
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
        
        # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format
        if(length(fileColumnName)>1) {
          targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
          targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
          targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
          for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
          targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
        } else {
          targets$file <- targets[,fileColumnName]
          targets$sample <- targets[,sampleColumnName]
        }
        
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
## ChIPhelper.Subread: parse Subread summary stats and create a md table
##
ChIPhelper.Subread <- function(subdir="",
                               sampleColumnName =c("IPname", "INPUTname"), 
                               fileColumnName =c("IP", "INPUT")) {
    
    FOLDER <- file.path(SHINYREPS_SUBREAD, subdir)
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
        
        # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format
        if(length(fileColumnName)>1) {
          targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
          targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
          targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
          for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
          targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
        } else {
          targets$file <- targets[,fileColumnName]
          targets$sample <- targets[,sampleColumnName]
        }
        
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
#' selectSampleSubset: select subset of samples for including in report (e.g. in case of multiple fastq files in scRNA-seq) 
#'
#' @param samples character vector with sample names. 
#' @param samplePattern regular expression to filter (include) on \code{samples}. No filtering if NULL or NA.
#' @param excludePattern regular expression to filter (exclude) on \code{samples}. No filtering if NULL or NA.
#' @param grepInBasename logical. If \code{TRUE} apply pattern to filename, not to full path.
#'
#' @return character vector with selected sample names
selectSampleSubset <- function(samples, samplePattern=NULL, excludePattern=NULL, grepInBasename=T, maxno=NULL) {
  
  if(!gtools::invalid(samplePattern)) {
    x <- if(grepInBasename) {basename(samples)} else {samples} # if TRUE apply pattern to filename, not to full path
    if(!any(grepl(samplePattern, x))) {
      warning(paste("\nYour pattern", samplePattern, "was not found in sample names and is ignored.\n")) 
    } else {
      samples <- samples[grep(samplePattern, x, invert=F)]
    }
    if(length(samples)==0) {stop("\nYou have selected no files!\n")}
  }
  
  if(!gtools::invalid(excludePattern)) {
    x <- if(grepInBasename) {basename(samples)} else {samples} # if TRUE apply pattern to filename, not to full path
    samples <- samples[grep(excludePattern, x, invert=T)]
    if(length(samples)==0) {stop("\nYou have selected no files!\n")}
  }
  
  if(!gtools::invalid(maxno) && is.numeric(maxno)) {
    samples <- samples[1:min(length(samples), maxno)]
    if(maxno > length(samples)) {cat("\nSample number restricted to", maxno)}
  }
  return(samples)
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



