#!/fsimb/groups/imb-bioinfocf/common-tools/dependencies/R/default/bin/Rscript
#####################################
##
## What: DE.edgeR.R
## Who : Sergi Sayols
## When: 06-10-2014
##
## Helper script to DE.edgeR.R to generate a design matrix based on a formula. The idea of
## this program is to help choosing beforehand the contrasts that will be tested in the
## contrasts.txt file
##
## Args:
## -----
## targets=targets.txt		# file describing the targets
## mmatrix=~group			# model matrix
## out=mmatrix.txt			# prefix filename for output (blank for stdout)
##
######################################
options(stringsAsFactors=FALSE)
library(edgeR,quietly=TRUE)

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {

	if(length(i <- grep(string,args,fixed=T)) == 1)
		return(do.call(convert,list(gsub(string,"",args[i]))))
	
	if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
ftargets   <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
mmatrix    <- parseArgs(args,"mmatrix=","~0+group")          # model matrix
out        <- parseArgs(args,"out=","")                    # output file

runstr <- "Rscript DE.edgeR.mmatrix.R [targets=targets.txt] [mmatrix=~0+group] [out=output]"
if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))

##
## create the design and contrasts matrix
##
# load targets
targets <- read.delim(ftargets,head=T,colClasses = "character",sep="",comment.char="#")
if(!any(grepl("sample",colnames(targets)))) {
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

if(out == "") {
	print(design)
} else {
	write.table(design,file=output)
}
