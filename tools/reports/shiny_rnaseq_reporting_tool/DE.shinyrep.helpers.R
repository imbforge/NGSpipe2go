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

##
## loadGlobalVars: read configuration from bpipe vars
##
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
# shorten: if a text string is longer than certain length, shorten by showing the first and last characters
shorten <- function(x, max.len=40, ini=20, end=15) {
    l <- nchar(x)
    if(l > max.len) paste(substr(x, 1, ini), substr(x, (l-end), l), sep="...") else x
}

##
## DESeq2 DE analysis
##
## DEhelper.MDS
DEhelper.DESeq2.MDS <- function() {
    p <- plotPCA(rld, intgroup=colnames(colData(dds))[1])
    print(p + geom_text_repel(aes(label=rownames(colData(dds)))) + theme_bw())
}

## DEhelper.cluster: Heatmap of top variant 'n' genes of the counts-per-milion table
DEhelper.DESeq2.cluster <- function(n=25) {
    rows <- order(apply(assay(rld), 1, sd), decreasing=TRUE)[1:n]
    hmcol  <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(assay(rld)[rows, ], col=hmcol, trace="none", margin=c(10, 6), scale="none")
}

## DEhelper.corr: Heatmap of sample to sample distances
DEhelper.DESeq2.corr <- function() {
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    hc  <- hclust(distsRL)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(mat, Rowv=as.dendrogram(hc), 
              symm=TRUE, trace="none", 
              col = rev(hmcol), margin=c(13, 13))
}

## DEhelper.MAplot: MA plots
DEhelper.DESeq2.MAplot <- function(i=1, fdr=.01) {
     plotMA(res[[i]], main=conts[i, 1])    
#    x <- mapply(function(res, cont) {
#        plotMA(res, main=cont)
#        invisible(0)
#    }, res, conts[, 1], SIMPLIFY=FALSE)
}

## DEhelper.DEgenes: show the DE results
DEhelper.DESeq2.DEgenes <- function(i=1) {
    ord  <- order(-log(res[[i]]$padj), 
                   abs(res[[i]]$log2FoldChange), 
                  decreasing=TRUE)
    res[[i]][ord, ]
}

## DEhelper.DESeq2.VolcanoPlot: Volcano plots from DEseq2 results
DEhelper.DESeq2.VolcanoPlot <- function(i=1, fdr=.01, top=25, web=TRUE) {
    # gather results
    d <- as.data.frame(DEhelper.DESeq2.DEgenes(i))
    x.limit <- max(abs(d$log2FoldChange), na.rm=T) # find the maximum spread of the x-axis
    
    # plotting
    p <- ggplot(d) +
            geom_point(mapping=aes(log2FoldChange, -log10(padj), size=log2(baseMean + 1), color=padj < fdr), alpha=.1) +
            theme_bw() +
            xlim(-x.limit, x.limit) +
            ylab("-log\u2081\u2080 adj. p-value") +
            xlab("log\u2082 fold change") + 
            scale_color_manual(values=c("black", "red"), guide=FALSE) +
            scale_size_continuous("mean counts (log\u2082)")

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
    x <- sapply(list.files(LOG, pattern=SUFFIX), function(f) {
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
    df <- data.frame(input_reads=format(x[1, ], big.mark=","), 
                     uniquely_mapped=paste0(format(x[2, ], big.mark=","), " (", format(x[3, ], nsmall=2), "%)"), 
                     multi_mapped=paste0(format(x[4, ], big.mark=","), " (", format(x[6, ], nsmall=2), "%)"), 
                     unmapped=paste0(format(x[7, ] + x[8, ] + x[9, ] + x[10, ], nsmall=2), "%"))
    kable(df, align=c("r", "r", "r", "r"), output=F)
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

    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
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
    
    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
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
    samplenames <- gsub(SHINYREPS_PREFIX, "", samplenames)
    samplenames <- gsub("_inferexperiment.txt", "", samplenames)
    colnames(strandspecifity) <- samplenames 
    rownames(strandspecifity) <- c("sense", "antisense", "other") 
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
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    table.content <- data.frame(  tool=c(
                                          "FastQC", 
                                          "STAR", 
                                          "Samtools", 
                                          "Subread", 
                                          "Picard", 
                                          "R/DEseq2"
                                          ), 
                                  version=c(
                                          Toolhelper.VersionReporter("FastQC",   SHINYREPS_FASTQC_LOG ), 
                                          Toolhelper.VersionReporter("STAR",     SHINYREPS_STAR_LOG ), 
                                          Toolhelper.VersionReporter("Samtools", SHINYREPS_BAMINDEX_LOG), 
                                          Toolhelper.VersionReporter("Subread",  SHINYREPS_SUBREAD_LOG), 
                                          Toolhelper.VersionReporter("Picard",   SHINYREPS_BAM2BW_LOG), 
                                          Toolhelper.VersionReporter("R/DEseq2", SHINYREPS_DESEQ_LOGS)
                                          )
                               )
    kable(table.content)
}

## version version of one tool, knowing its log folder
Toolhelper.VersionReporter <- function(tool, logfolder) {
    
    LOG <- logfolder
    SUFFIX <- paste0(".log", "$")
    
    # logs folder
    if(!file.exists(LOG)) {
        return(paste0(tool, " version not available"))
    }
    
    x <- lapply( list.files(LOG, pattern=SUFFIX, full.names=TRUE), function(f){
        # read all lines
        l <- readLines(f)
        # need to check Version number in one line lower than "VERSION INFO"
        # e.g. FastQC v0.11.3
        l.version <- l[ grep("^VERSION INFO", l) + 1 ]
        
        return(l.version)
        
        } )
    
    # x is a list of always the same content
    r <- tryCatch(
        {
            if (is.null(x[[1]][1])) {
                return("no version tag")
            } else {
                return(x[[1]][1])
            }
        }, 
        warning = function(w) {
            return("no version tag")
        }, 
        error = function(e) {
            return("no version tag")
        }, 
        finally = {}
    )
    
}
