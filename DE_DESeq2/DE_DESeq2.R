#####################################
##
## What: DE_DESeq2.R
## Who : Sergi Sayols
## When: 29-01-2015
##
## Script to perform DE different conditions, based on DESeq2 Negative Binomial model
## and the counts from htseq-count.
##
## Args:
## -----
## targets=targets.txt		# file describing the targets
## contrasts=contrasts.txt  # file describing the contrasts
## mmatrix=~group			# model matrix. Wrapper constrained to always use an intercept in the model
## filter=NOTUSED			# filter invariant genes? Always TRUE (DESeq2 default)
## prefix=RE				# prefix to remove from the sample name
## suffix=RE				# suffix to remove from the sample name (usually _readcounts.tsv)
## cwd=.					# current working directory where the files .tsv files are located
## out=DE.DESeq2			# prefix filename for output
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
mmatrix      <- parseArgs(args,"mmatrix=","~group")          # model matrix (or ~0+group for multiple comparisons)
filter.genes <- parseArgs(args,"filter=",TRUE,convert="as.logical") # filter invariant genes?
pre          <- parseArgs(args,"prefix=","")    # prefix to remove from the sample name
suf          <- parseArgs(args,"suffix=","_readcounts.tsv")    # suffix to remove from the sample name
cwd          <- parseArgs(args,"cwd=","./")     # current working directory
out          <- parseArgs(args,"out=","DE.DESeq2") # output filename

runstr <- "Rscript DE.DESeq2.R [targets=targets.txt] [contrasts=contrasts.txt] [mmatrix=~0+group] [filter=TRUE] [prefix=RE] [suffix=RE] [cwd=.] [out=DE.edgeR]"
if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(fcontrasts)) stop(paste("File",fcontrasts,"does NOT exist. Run with:\n",runstr))
if(!file.exists(cwd))        stop(paste("Dir",cwd,"does NOT exist. Run with:\n",runstr))
#if(is.na(filter.genes))      stop(paste("Filter (filter invariant genes) has to be either TRUE or FALSE. Run with:\n",runstr))

##
## create the design and contrasts matrix
##
# load targets
targets <- read.delim(ftargets,head=T,colClasses="character",comment.char="#",)
if(!any(grepl("sample",colnames(targets)))) {
	targets$sample <- gsub(paste0("^",pre),"",targets[,1])	  # clean file names to construct the sample names
	targets$sample <- sub(paste0(suf,"$"),"",targets$sample)
}
conds  <- sapply(colnames(targets),grepl,mmatrix)	# total number of factors in the formula (~0+A+A:B or ~group)

# load contrasts
conts <- read.delim(fcontrasts,head=F,comment.char="#")

##
## read count tables from htseq-count
##
# check if the specified targets exists in the CWD
if(!all(x <- sapply(paste0(cwd,"/",targets$file),file.exists))) stop(paste("one or more input files do not exist in",cwd,"(x=",x,")"))

# read input HTseq counts
targets <- data.frame(sampleName=targets$sample,
					  fileName=targets$file,
					  targets[,conds,drop=F])

dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets, 
								  directory=cwd, 
								  design=as.formula(mmatrix))

if(!attr(terms(design(dds)),"intercept")) stop("cannot deal with designs without the intercept")

rld <- rlog(dds)

##
## DESeq analysis: right now it only allows simple linear models with pairwise comparisons
##
dds <- DESeq(dds)
res <- lapply(conts[,1],function(cont) {
	# parse the formula in cont, get contrasts from resultNames(dds) and create a vector with coefficients
	cont.name <- gsub("(.+)=(.+)","\\1",cont)
	cont.form <- gsub("(.+)=(.+)","\\2",cont)
	factors   <- unlist(strsplit(cont.form,"\\W"))
	factors   <- factors[factors != ""]
	if(length(factors) != 2) { # || sum(conds) != 1) {
		warning(paste(cont,"cannot deal with designs other than pairwise comparisons!"))
		return(NA)
	}

	# DEseq2 stats
	res <- results(dds,contrast=list(factors[1],factors[2]))
	res <- res[order(res$padj),]

	# write the results
	factors <- Reduce(function(x,y) gsub(y,"",x),paste0("^",names(which(conds))),init=factors)
	samples <- lapply(paste0("dds$",names(which(conds))),function(x) as.character(eval(parse(text=x))) %in% factors)
	samples <- if(length(samples) == 1) unlist(samples) else do.call("&",samples)
	x <- merge(res,assay(rld)[,samples],by=0)
	write.csv(x[order(x$padj),],file=paste(out,cont.name,"csv",sep="."),row.names=F)
	res
})

##
## Plots
##
pdf(paste0(out,".pdf"))

# sample to sample distance heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
hmcol  <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, Rowv=as.dendrogram(hc),
		  symm=TRUE, trace="none",
		  col = rev(hmcol), margin=c(13, 13))

# PCA
plotPCA(rld,intgroup=names(which(conds)))

res <- lapply(1:length(conts[,1]),function(i) {
	cont <- conts[i,]

	# parse the formula in cont, get contrasts from resultNames(dds) and create a vector with coefficients
	cont.name <- gsub("(.+)=(.+)","\\1",cont)
	cont.form <- gsub("(.+)=(.+)","\\2",cont)
	factors   <- unlist(strsplit(cont.form,"\\W"))
	factors   <- factors[factors != ""]
	factors <- Reduce(function(x,y) gsub(y,"",x),paste0("^",names(which(conds))),init=factors)
	samples <- lapply(paste0("dds$",names(which(conds))),function(x) as.character(eval(parse(text=x))) %in% factors)
	samples <- if(length(samples) == 1) unlist(samples) else do.call("&",samples)

	# MA plot
	plotMA(res[[i]],main=cont.name,ylim=c(-2,2))

	# heatmap of the top variant genes
	rows <- order(apply(counts(dds[,samples],normalized=TRUE),1,sd),decreasing=TRUE)[1:30]
	hmcol  <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
	heatmap.2(assay(rld)[rows,samples],col=hmcol,trace="none",margin=c(10,6),scale="none")

	# sample to sample distance heatmap
	distsRL <- dist(t(assay(rld)[,samples]))
	mat <- as.matrix(distsRL)
	hc <- hclust(distsRL)
	heatmap.2(mat, Rowv=as.dendrogram(hc),
			  symm=TRUE, trace="none",
			  col = rev(hmcol), margin=c(13, 13))

	# PCA
	plotPCA(rld[,samples],intgroup=names(which(conds)))
})

dev.off()
save.image(file=paste0(out,".RData"))		# for further reports
