#####################################
##
## What: DE_DESeq2.R
## Who : Sergi Sayols
## When: 29-01-2015
##
## Script to perform DE different conditions, based on DESeq2 Negative Binomial model
## and the counts from htseq-count.
## Only supports pairwise comparisons, but any formula could be provided, including blocking factors.
##
## Args:
## -----
## targets=targets.txt      # file describing the targets.
##                          # Must fit the format expected in DESeqDataSetFromHTSeqCount
## contrasts=contrasts.txt  # file describing the contrasts
## mmatrix=~condition       # model matrix. Wrapper constrained to always use an intercept in the model
## gtf=gene_model.gtf       # gene model in gtf format, for rpkm calculation
## filter=TRUE              # perform automatic independent filtering of lowly expressed genes to maximise power
## prefix=RE                # prefix to remove from the sample name
## suffix=RE                # suffix to remove from the sample name (usually _readcounts.tsv)
## cwd=.                    # current working directory where the files .tsv files are located
## out=DE.DESeq2            # prefix filename for output
##
## IMPORTANT: This is a simplified wrapper to DESeq2 which is only able to do Wald tests on
##            simple experiments. It's meant only for pairwise comparisons in non-multifactor
##            designs. Why? to avoid messing our standard pipeline with extra parms, reusing
##            all the edgeR parms from the pipeline. With this simplification come along all
##            the limitations stated here.
##
######################################
options(stringsAsFactors=FALSE)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(pheatmap)
library(viridis)

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {

    if(length(i <- grep(string,args,fixed=T)) == 1)
        return(do.call(convert,list(gsub(string,"",args[i]))))
    
    if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
ftargets     <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
fcontrasts   <- parseArgs(args,"contrasts=","contrasts.txt") # file describing the contrasts
mmatrix      <- parseArgs(args,"mmatrix=","~condition")      # model matrix
gene.model   <- parseArgs(args,"gtf=","")       # gtf gene model
filter.genes <- parseArgs(args,"filter=",TRUE,convert="as.logical") # automatic independent filtering
pre          <- parseArgs(args,"prefix=","")    # prefix to remove from the sample name
suf          <- parseArgs(args,"suffix=","_readcounts.tsv")    # suffix to remove from the sample name
cwd          <- parseArgs(args,"cwd=","./")     # current working directory
out          <- parseArgs(args,"out=","DE.DESeq2") # output filename
pattern      <- parseArgs(args,"pattern=","\\.readcounts.tsv") # output filename
FC           <- parseArgs(args, "FC=", 1 , convert="as.numeric") # FC filter non log2 
FDR           <- parseArgs(args, "FDR=", 0.01 , convert="as.numeric") # filter FDR 

runstr <- "Rscript DE.DESeq2.R [targets=targets.txt] [contrasts=contrasts.txt] [mmatrix=~condition] [gtf=] [filter=TRUE] [prefix=RE] [suffix=RE] [cwd=.] [base=] [out=DE.DESeq2] [pattern=RE] [FC=1] [FDR=0.01]"
if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(fcontrasts)) stop(paste("File",fcontrasts,"does NOT exist. Run with:\n",runstr))
if(!file.exists(cwd))        stop(paste("Dir",cwd,"does NOT exist. Run with:\n",runstr))
if(is.na(filter.genes))      stop(paste("Filter genes has to be either TRUE or FALSE. Run with:\n",runstr))
if(!file.exists(gene.model)) stop(paste("GTF File:", gene.model, " does NOT exist. Run with: \n", runstr))

##
## create the design and contrasts matrix
##
# load targets
targets <- read.delim(ftargets,head=T,colClasses="character",comment.char="#")
if(!all(c("group", "file", "sample") %in% colnames(targets))) stop("targets file must have at least 3 columns and fit the format expected in DESeqDatcondition")

# reorder the targets file
add_factors <- colnames(targets)[!colnames(targets) %in% c("group", "sample", "file")]
targets <- targets[, c("sample", "file", "group", add_factors)]

# grep sample identifier in count file names
  countfiles <- list.files(cwd)
  countfiles <- countfiles[grep(pattern, countfiles)] # filter for valid count files
  
  index_targetsfile <- sapply(targets$file, grep,  countfiles) # grep targets in countfiles

  ## check matching ambiguity
      if(class(index_targetsfile)=="list") { # list means either zero or multiple matches
        if(any(x <- sapply(index_targetsfile, length)>1)) {
          stop(paste("\nA targets.txt entry matches multiple count file names\ntargets.txt: "),
               paste(targets$file[x], collapse=", "),
               "\ncount file names: ", paste(countfiles[unlist(index_targetsfile[x])], collapse=", "))
        }
        warning(paste("Entries in targets.txt which are not found in list of count files are removed from targets.txt:", 
                      targets$file[sapply(index_targetsfile, length)==0], collapse=", "))
        targets <- targets[!(sapply(index_targetsfile, length)==0),] # remove target entries
        index_targetsfile <- sapply(targets$file, grep,  countfiles) # recreate index vector after removal of targets.txt entries
      }
      
      if(any(x <- duplicated(index_targetsfile))) { # check for multiple target entries matching the same file name
        stop(paste("\nMultiple targets.txt entries match to the same count file name\ntargets.txt: ", 
                   paste(targets$file[x | duplicated(index_targetsfile, fromLast=T)], collapse=", "),
                   "\ncount file names:", paste(countfiles[unlist(index_targetsfile[x])], collapse=", "), "\n"))
      }
      
      if (length(unique(index_targetsfile)) < length(countfiles)) { # check for count file names not included
        warning(paste("Count file names not included in targets.txt are ignored: ",  paste(countfiles[-index_targetsfile], collapse=", ")))
      }
  
  targets$file <- countfiles[index_targetsfile] # replace entries in targets$file by count file names
  


# load contrasts
conts <- read.delim(fcontrasts,head=F,comment.char="#")

##
## calculate gene lengths
##
gtf <- import.gff(gene.model, format="gtf", feature.type="exon")
gtf.flat <- unlist(reduce(split(gtf, elementMetadata(gtf)$gene_id)))
gene.lengths <- tapply(width(gtf.flat), names(gtf.flat), sum)

# make sure we have a gene_name column (needed as output in the report later)
if(! "gene_name" %in% colnames(mcols(gtf))) {
    if("gene_id" %in% colnames(mcols(gtf))) {
        gtf$gene_name <- gtf$gene_id
    } else {
        gtf$gene_name <- NA
    }
}

# create a txdb object to collect the genes coordinates for later usage
txdb  <- makeTxDbFromGRanges(gtf)
genes <- as.data.frame(genes(txdb))

##
## DESeq analysis: right now it only allows simple linear models with pairwise comparisons
##
pairwise.dds.and.res <- lapply(conts[,1],function(cont) {

    # parse the formula in cont, get contrasts from resultNames(dds) and create a vector with coefficients
    cont.name <- gsub("(.+)=(.+)","\\1",cont)
    cont.form <- gsub("(.+)=(.+)","\\2",cont)
    factors   <- gsub("(^\\s+|\\s+$)", "", unlist(strsplit(cont.form,"\\W")))
    factors   <- factors[factors != ""]
    if(length(factors) != 2) {
        warning(paste(cont,"cannot deal with designs other than pairwise comparisons!"))
        return(NA)
    }

    # read input HTseq counts
    this_targets <- targets[targets$group %in% factors,]
    dds <- DESeqDataSetFromHTSeqCount(sampleTable=this_targets,
                                      directory=cwd, 
                                      design=as.formula(mmatrix))

    # relevel to make sure the comparison is done in the right direction (B vs A and not A vs B)
    dds$group <- relevel(dds$group, factors[2])

    # create DESeq object
    dds <- DESeq(dds)
    #we call it with the specified threshold for FC and for FDR
    res <- results(dds, 
                   independentFiltering=filter.genes,
                   format="DataFrame",
                   lfcThreshold =log2(FC),
                   alpha = FDR)

    # calculate quantification (RPKM if gene model provided, rlog transformed values otherwise)
    quantification <- apply(fpm(dds), 2, function(x, y) 1e3 * x / y, gene.lengths[rownames(fpm(dds))])
    
    # add comment "robustRPKM" to columns, such that it's clear what the value represents
    colnames(quantification) <- paste0(colnames(quantification),".robustRPKM")

    # extract the gene_name and genomic coordinates of each gene
    res$gene_name <- gtf$gene_name[match(rownames(res), gtf$gene_id)]
    i <- match(rownames(res), genes$gene_id)
    res$chr    <- genes$seqnames[i]
    res$start  <- genes$start[i]
    res$end    <- genes$end[i]
    res$strand <- genes$strand[i]
    
    # write the results
    x <- merge(res[, c("gene_name", "chr", "start", "end", "strand", "baseMean", "log2FoldChange", "padj")], quantification, by=0)
    x <- x[order(x$padj),]

    #separate the data from x into the tested genes, the upregulated genes and the downregulated genes
    tested_genes <- x[!is.na(x$padj),]
    #create the output list 
    x_info <- list(up     = tested_genes[tested_genes$log2FoldChange > log2(FC) & tested_genes$padj < FDR, ],
                   down   = tested_genes[tested_genes$log2FoldChange < log2(FC) & tested_genes$padj < FDR, ],
                   tested = tested_genes,
                   all_genes = x)
    
    colnames(x)[which(colnames(x) %in% c("baseMean", "log2FoldChange", "padj"))] <- mcols(res)$description[match(c("baseMean", "log2FoldChange", "padj"), colnames(res))]
    colnames(x)[1] <- "gene_id"

    write.csv(x, file=paste0(out, "/", cont.name, ".csv"), row.names=F)
    #we also have to replace the columnnames within the x_info list
    x_info <- lapply(x_info, function(y){
                       colnames(y)[which(colnames(y) %in% c("baseMean", "log2FoldChange", "padj"))] <- mcols(res)$description[match(c("baseMean", "log2FoldChange", "padj"), colnames(res))]
                       colnames(y)[1] <- "gene_id"
                       return(y)
                   })
    #add a description before writing out the the excel file
    x_info[["Description"]] <-data.frame(Descripton=paste("This DESeq2 analysis was performed using a FC filter of ", FC, "and a filter for the adjusted p-value of ", FDR, "."))
    write.xlsx(x_info, file=paste0(out, "/", cont.name, ".xlsx"), row.names=F)
    
    list(dds,res)
})

## separate pairwise PCAs to a list
pairwise.dds <- lapply(pairwise.dds.and.res,function(x){return(x[[1]])})

## separate results to a list
res <- lapply(pairwise.dds.and.res,function(x){return(x[[2]])})
names(res) <-  gsub("=.*", "", conts[,1]) 


##
## Sanity check plots with all the samples together
##
pdf(paste0(out,"/DE_DESeq2.pdf"))

# make a DESeq2 object with all the samples
dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
                                  directory=cwd, 
                                  design=as.formula(mmatrix))
dds <- DESeq(dds)
rld <- rlog(dds)

## quantify all samples at the same time (necessary to get unique FPM per sample per gene, since robust FPMs
## depend on all samples in the data frame, not just on individual samples

# quantify
quantification <- as.data.frame(apply(fpm(dds,robust=TRUE), 2, function(x, y) 1e3 * x / y, gene.lengths[rownames(fpm(dds,robust=TRUE))]))

# add comment "robustRPKM" to columns, such that it's clear what the value represents
names(quantification) <- paste0(names(quantification),".robustRPKM")

# extract the gene_name and genomic coordinates of each gene
names.df <- data.frame(gene_name=gtf$gene_name[match(rownames(quantification), gtf$gene_id)],row.names=rownames(quantification))
i <- match(rownames(names.df), genes$gene_id)
names.df$chr    <- genes$seqnames[i]
names.df$start  <- genes$start[i]
names.df$end    <- genes$end[i]
names.df$strand <- genes$strand[i]

# merge location and quantification
quantification.names.df <- merge(names.df,quantification,by=0)
colnames(quantification.names.df)[1] <- "gene_id"

write.csv(quantification.names.df, file=paste0(out, "/allSamples.robustRPKM.csv"), row.names=F)
write.xlsx(quantification.names.df, file=paste0(out, "/allSamples.robustRPKM.xlsx"), row.names=F)

# extract rlog assay and change to user friendly gene identifiers
assay.rld <- assay(rld)
rownames(assay.rld) <- gtf$gene_name[match(rownames(assay(rld)), gtf$gene_id)]

# extract information for legend
#legend.df <- as.data.frame(colData(rld)[,c("group","subject")])
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

# sample to sample distance heatmap
sampleDists <- dist(t(assay.rld))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_col=legend.df,
         annotation_colors=legend_colors,
         col=plasma(255),
         border_color=NA,
         main="Sample to sample distances",
         fontsize=6,
	 treeheight_row=20,
         treeheight_col=20,
         annotation_names_col=FALSE)

# sample to sample distance heatmap w/ 0's on diagonal set to NA
# since they otherwise shift the scale too much
sampleDistMatrix.diagNA <- sampleDistMatrix
sampleDistMatrix.diagNA[sampleDistMatrix.diagNA==0] <- NA
pheatmap(sampleDistMatrix.diagNA,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_col=legend.df,
         annotation_colors=legend_colors,
         col=plasma(255),
         border_color=NA,
         main="Sample to sample distances (diagonal set to NA)",
         fontsize=6,
         treeheight_row=20,
         treeheight_col=20,
 	 annotation_names_col=FALSE)

# heatmap of the top variant genes
select.highestSD <- order(apply(assay.rld,1,sd),decreasing=TRUE)[1:40]
pheatmap(assay.rld[select.highestSD,], 
         cluster_rows=TRUE, 
         cluster_cols=TRUE, 
         show_rownames=TRUE, 
         annotation_col=legend.df,
         annotation_colors=legend_colors,
         color=colorRampPalette(brewer.pal(9,"GnBu"))(255),
         border_color=NA, 
         main="Normalized expression values of 40 most variable genes",
	 fontsize=6,
         fontsize_col=10,
         fontsize_row=6,
	 treeheight_row=20,
         treeheight_col=20,
         annotation_names_col=FALSE)

# heatmap of top mean genes
select.highestMean <- order(rowMeans(assay.rld), decreasing=TRUE)[1:40]
pheatmap(assay.rld[select.highestMean,], 
         cluster_rows=FALSE, 
         cluster_cols=TRUE,
         show_rownames=TRUE, 
         annotation_col=legend.df,
         annotation_colors=legend_colors,
         col=rev(heat.colors(255)),
         border_color=NA, 
         main="Normalized expression values of 40 genes with highest mean",
	 fontsize=6,
         fontsize_col=10,
         fontsize_row=6,
	 treeheight_row=20,
         treeheight_col=20,
         annotation_names_col=FALSE)

# PCA, group based on the first factor (column 3 in targets)
p <- plotPCA(rld, intgroup=colnames(colData(dds))[1])
#plot(p + geom_text_repel(aes(label=rownames(colData(dds)))) + theme_bw())
plot(p + 
     scale_color_manual(values=brewer.pal(9,"Set1")[1:length(levels(colData(dds)[,"group"]))]) +
     geom_text_repel(aes(label=rownames(colData(dds))), show.legend=FALSE) + 
     theme_bw())


# MA plot
x <- mapply(function(res, cont) {
    plotMA(res, main=cont, ylim=c(-3,3))
    invisible(0)
}, res, conts[, 1], SIMPLIFY=FALSE)

dev.off()
#save the sessionInformation
writeLines(capture.output(sessionInfo()),paste(out, "/DE_DESeq2_session_info.txt", sep=""))
save(dds, rld, res, pairwise.dds, conts, gtf, add_factors, file=paste0(out,"/DE_DESeq2.RData"))
