#####################################
##
## What: DE.edgeR.R
## Who : Sergi Sayols
## When: 06-10-2014
##
## Script to perform DE different conditions, based on edgeR Negative Binomial model
## and the counts from htseq-count.
##
## Args:
## -----
## targets=targets.txt		# file describing the targets
## contrasts=contrasts.txt  # file describing the contrasts
## mmatrix=~0+group			# model matrix
## filter=TRUE				# filter invariant genes?
## prefix=RE				# prefix to remove from the sample name
## suffix=RE				# suffix to remove from the sample name (usually  _readcounts.tsv)
## cwd=.					# current working directory where the files .tsv files are located
## robust=FALSE             # robustly estimate dispersion?
## out=DE.edgeR				# prefix filename for output
##
## IMPORTANT: due to the difficulties in expressing the contrasts in complex designs, one should
## ---------  use the helper script DE.edgeR.mmatrix.R to generate beforehand a design matrix
##            based on a formula. With this script one can properly define the comparisons in 
##            the file contrasts.txt
##
######################################
options(stringsAsFactors=FALSE)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(ggplot2)

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
mmatrix      <- parseArgs(args,"mmatrix=","~0+group")        # model matrix (or ~0+group for multiple comparisons)
filter.genes <- parseArgs(args,"filter=",TRUE,convert="as.logical") # filter invariant genes?
pre          <- parseArgs(args,"prefix=","")    # prefix to remove from the sample name
suf          <- parseArgs(args,"suffix=","_readcounts.tsv")    # suffix to remove from the sample name
cwd          <- parseArgs(args,"cwd=","./")     # current working directory
robust       <- parseArgs(args,"robust=",FALSE,convert="as.logical") # robustly estimate dispersion
out          <- parseArgs(args,"out=","DE.edgeR") # output filename

runstr <- "Rscript DE.edgeR.R [targets=targets.txt] [contrasts=contrasts.txt] [mmatrix=~0+group] [filter=TRUE] [prefix=RE] [suffix=RE] [cwd=.] [robust=FALSE] [out=DE.edgeR]"
if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(fcontrasts)) stop(paste("File",fcontrasts,"does NOT exist. Run with:\n",runstr))
if(!file.exists(cwd))        stop(paste("Dir",cwd,"does NOT exist. Run with:\n",runstr))
if(is.na(filter.genes))      stop(paste("Filter (filter invariant genes) has to be either TRUE or FALSE. Run with:\n",runstr))
if(is.na(robust))            stop(paste("Robust (robustly estimate dispersion) has to be either TRUE or FALSE. Run with:\n",runstr))

##
## create the design and contrasts matrix
##
# load targets
targets <- read.delim(ftargets,head=T,colClasses="character",comment.char="#",sep="")
if(!any(grepl("^sample$",colnames(targets)))) {
	targets$sample <- gsub(paste0("^",pre),"",targets[,1])	  # clean file names to construct the sample names
	targets$sample <- sub(paste0(suf,"$"),"",targets$sample)
}

# create the model matrix from the formula provided
design <- model.matrix(as.formula(mmatrix),data=targets)
rownames(design) <- targets$sample
conds  <- sapply(colnames(targets),grepl,mmatrix)	# total number of factors in the formula (~0+A+A:B or ~group)

# clean the colnames in the design matrix that contain the targets.txt colname
#for(x in colnames(targets)[conds]) {
#	colnames(design) <- gsub(paste0("^",x),"",colnames(design))
#	colnames(design) <- gsub(paste0("(.+:)",x),"\\1",colnames(design))
#	if(length(i <- grep("Intercept",colnames(design))) > 0) {
#		colnames(design)[-i] <- make.names(colnames(design)[-i]) # make sintactically valid names after gsub
#	} else {
#		colnames(design) <- make.names(colnames(design)) # make sintactically valid names after gsub
#	}
#}

# define the group the samples belong to if there was only 1 categorical factor in the formula
if(sum(conds) == 1) {
	group <- as.factor(targets[,conds])
#	colnames(design) <- levels(group)
} else {
#	group <- as.factor(rep(1,length(targets$sample)))
#	colnames(design) <- levels(group)
	group <- as.factor(targets[,which(conds)[1]])
}

# load contrasts
x <- read.delim(fcontrasts,head=F,comment.char="#")
conts <- makeContrasts(contrasts=x[,1],levels=design)
colnames(conts) <- gsub("\\s*=.+","",x[,1])

##
## read count tables from htseq-count and subread counts transformed to htseq format
##
# check if the specified targets exists in the CWD
if(!all(x <- sapply(paste0(cwd,"/",targets$file),file.exists))) stop(paste("one or more input files do not exist in",cwd,"(x=",x,")"))

# read files
countTable <- lapply(paste0(cwd,"/",targets$file),read.table,header=F,row.names=1)

# check, if all files are sorted the same way
# this is inefficient code and needs to be improved
tmp.rownames <- lapply( countTable, row.names )
check.identity <- c()
for ( i in seq(length(tmp.rownames)) ) {
	check.identity <- c( check.identity, identical( tmp.rownames[[1]], tmp.rownames[[i]]) ) 
}
if ( ! all(check.identity) ) {
	print( targets$file[!check.identity] )
	# die!!!
	stop(paste("File sorting seems to be different in the mentioned files, compared to: ", targets$file[[1]], sep=''))
}
rm(tmp.rownames)
rm(check.identity)

# continue reading files
countTable <- Reduce(cbind,countTable)
colnames(countTable) <- targets$sample

# Filter out non-gene counts from htseq-count. Versions prior to 0.5.4 didn't include the doble underscore
countTable <- countTable[!(rownames(countTable) == "ambiguous"),]
countTable <- countTable[!(rownames(countTable) == "no_feature"),]
countTable <- countTable[!(rownames(countTable) == "alignment_not_unique"),]
countTable <- countTable[!(rownames(countTable) == "too_low_aQual"),]
countTable <- countTable[!(rownames(countTable) == "not_aligned"),]
countTable <- countTable[grep("^__",rownames(countTable),invert=T),]

##
## Create digital gene expression, estimate dispersions and calculate linear model for DE
##
# DGE and compute effective library size estimation by TMM normalization
if(sum(conds) == 1) {
	y <- DGEList(counts=countTable,group=group)
} else {
	y <- DGEList(counts=countTable)
}
y <- calcNormFactors(y)

# filter uniformative genes
if(filter.genes) {	
	m <- 1e6 * t(t(y$counts) / y$samples$lib.size)	# per library size and million reads normalized counts matrix
	y <- y[rowSums(m > 1) >= 2,] 	# select genes with RPM > 1 in more than 1 sample
}

# dispersions
y <- estimateGLMCommonDisp(y,design)	# overall dispersion for the dataset
if(robust) {
	y <- estimateGLMRobustDisp(y,design)	# Compute a robust estimate of the negative binomial dispersion
} else {
	y <- estimateGLMTrendedDisp(y,design)	# gene-wise dispersion estimates
	y <- estimateGLMTagwiseDisp(y,design)	# tagwise dipersion, strongly recommended in multifactor designs
}

# call DE (fit model and likelihood ratio test)
fit <- glmFit(y,design)	# fit model
lrt <- lapply(colnames(conts),function(x) glmLRT(fit,contrast=conts[,x]))
names(lrt) <- colnames(conts)

##
## sanity checks: TO BE REMOVED WHEN DOING SHINY AND/OR MARKDOWN REPORTS
##
pdf(paste0(out,".pdf"))
m <- cpm(y,prior.count=2,log=TRUE) # get scaled counts per million
v <- apply(m,1,sd,na.rm=T)		

# MDS and variance estimation
plotMDS.DGEList(y,cex=.5,col=brewer.pal(length(levels(group)),"Set1")[group])
plotBCV(y)	# variance along log gene count-per-milion

# Heatmap of top variant 50 genes of the counts-per-milion table
hmcol <- colorRampPalette(brewer.pal(9,"Oranges"))(100)
heatmap.2(m[rev(order(v))[1:50],],col=hmcol,trace="none",margin=c(10,6))

# Heatmap of sample to sample distances
dists <- dist(t(m))
mat <- as.matrix(dists)
heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(13,13))

# MA plots
lapply(1:length(lrt),function(i) {
	
	# get DE genes (p.adjust='BH', pval<.05)
	de <- decideTestsDGE(lrt[[i]])
	degenes <- rownames(y)[as.logical(de)]

	# MA plot
	plotSmear(lrt[[i]],de.tags=degenes,main=names(lrt)[[i]])	# MA plot
	abline(h=c(-1,1),col="blue")	# indicate 2-fold changes in the MA plot
	abline(v=0,col="blue")			# indicate >1 counts-per-million

	# save the DE results into a file (whole gene list in a csv file, 50 top in a tex file)
	groups  <- rownames(conts)[conts[,i] != 0]	# the groups involved in this contrast
	samples <- rownames(design)[apply(design[,groups],1,sum) > 0]	# samples belonging to the involved groups
	x <- merge(m[,samples],lrt[[i]]$table,by=0)
	x$FDR <- p.adjust(x$PValue,method="fdr")
	write.csv(x[order(x$FDR),],file=paste(out,names(lrt)[i],"csv",sep="."),row.names=F)

	invisible(0)
})

dev.off()
save.image(file=paste0(out,".RData"))		# for further reports
