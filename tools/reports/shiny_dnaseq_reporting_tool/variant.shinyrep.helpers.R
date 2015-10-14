##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
# library("edgeR")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("knitr")		# for markdown output

##
## loadGlobalVars: read configuration from bpipe vars
##
loadGlobalVars <- function(f="shinyReports.txt") {

	# read in the conf file
	conf <- readLines(f)
	conf <- conf[grep("^SHINYREPS_",conf)]
	
	# create the vars
	sapply(conf,function(x) {
		x <- unlist(strsplit(x,"=",fixed=T))
		assign(x[1],x[2],envir=.GlobalEnv)
	})
	
	invisible(0)
}

##
## DEhelper.init: some time consuming tasks that can be done in advance
##
VARhelper.init <- function(task) {
	
	# Prepare the DE data frame
	renderUcscGeneLinks <- function() {
		ucsc_url <- paste0("http://genome.ucsc.edu/cgi-bin/hgGene?org=",SHINYREPS_ORG,"&db=",SHINYREPS_DB,"&hgg_gene=")
		for(i in 1:length(lrt)) {
			lrt[[i]]$table$gene <<- sapply(rownames(lrt[[i]]$table),function(x) {
				paste0("<a href=\"",ucsc_url,x,"\">",x,"</a>")
			})
		}
	}
	prepareDEdataTable <- function() {
		for(i in 1:length(lrt)) {
			lrt[[i]]$table$FDR    <<- p.adjust(lrt[[i]]$table$PValue,method="fdr")
			lrt[[i]]$table$logFC  <<- round(lrt[[i]]$table$logFC,2)
			lrt[[i]]$table$logCPM <<- round(lrt[[i]]$table$logCPM,2)
			lrt[[i]]$table$LR     <<- round(lrt[[i]]$table$LR,2)
			lrt[[i]]$table$PValue <<- round(lrt[[i]]$table$PValue,4)
			lrt[[i]]$table$FDR    <<- round(lrt[[i]]$table$FDR,4)
		}
	}
	
	# Cluster and correlation tasks
	prepareDistanceMatrix <- function() {
		v <<- apply(m,1,sd,na.rm=T)		# get top variant genes
		dists <<- dist(t(m))
		mat <<- as.matrix(dists)
		hmcol <<- colorRampPalette(brewer.pal(max(length(levels(group)),3),"Oranges"))(100)
	}
	
	# dispatch tasks
	switch(task,
		   renderUcscGeneLinks=renderUcscGeneLinks(),
		   prepareDEdataTable=prepareDEdataTable(),
		   prepareDistanceMatrix=prepareDistanceMatrix())
}

##
## VARhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
VARhelper.Fastqc <- function(web=TRUE) {
	
	# logs folder
	if(!file.exists(SHINYREPS_FASTQC_LOG)) {
		return("Fastqc statistics not available")
	}
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/fastqc" else SHINYREPS_FASTQC_LOG
	
	# construct the image url from the folder ents (skip current dir .)
	samples <- list.dirs(SHINYREPS_FASTQC_LOG,recursive=F)
	df <- sapply(samples,function(f) {
		c(paste0("![alt text](",QC,"/",basename(f),"/Images/duplication_levels.png)"), 
		  paste0("![alt text](",QC,"/",basename(f),"/Images/per_base_quality.png)"), 
		  paste0("![alt text](",QC,"/",basename(f),"/Images/per_base_sequence_content.png)"))
	})

	# set row and column names, and output the md table
	df <- as.data.frame(t(df))
	rownames(df) <- gsub(paste0("^",SHINYREPS_PREFIX),"",basename(samples))
	colnames(df) <- c("Duplication","Read qualities","Sequence bias")
	kable(df,output=F)
}

##
## VARhelper.BWA: parse BWA mem and samtools flagstat output
##
VARhelper.BWA <- function() {

	# log file, which was copied from .bpipe folder
	# contains the runtime STDOUT of BWA and the samtools flagstat STDOUT
	LOG <- SHINYREPS_BWA_LOG
	SUFFIX <- paste0(SHINYREPS_BWA_SUFFIX, '$')
	
	if(!file.exists(LOG)) {
		return("BWA statistics not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, full.names=TRUE),function(f) { # list all files and feed them into function one by one
		l <- readLines(f) # read file content to l
		#close(f)
		
		sapply(c("in total",                             #1
				 "secondary",                            #2
				 "mapped \\(",                           #3
				 "paired in sequencing",                 #4
				 "read1",                                #5
				 "read2",                                #6
				 "properly paired",                      #7
				 "with itself and mate mapped",          #8
				 "singletons",                           #9
				 "with mate mapped to a different chr$", #10
				 "with mate mapped to a different chr \\(mapQ>=5\\)"),function(y) {   #11
				 	as.numeric(gsub("(^\\d+).+","\\1",l[grep(y,l)])) # grep returns line number, then get the respective line ([]) and extract the first number out of it (gsub and replace the whole line with it)
				 })	
	})
	
	# transform x from list to matrix (in extreme cases also with only one column)
	x <- do.call(cbind, x)
	# set row and column names, and output the md table
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	df <- data.frame(total=x[1,],
					 mapped=paste0( x[3,] ," (", round(x[3,] / x[1,] * 100, digits=2), "%)" ),
					 proper_pair=paste0( x[7,] ," (", round(x[7,] / x[1,] * 100, digits=2), "%)" ) ,
					 secondary_alignments=paste0( x[2,]," (", round(x[2,] / x[1,] * 100, digits=2), "%)" ),
					 unmapped=paste0( x[3,] - x[4,]," (", round((x[3,] - x[4,]) / x[4,] * 100, digits=2), "%)" ),
					 different_chromosome=paste0( x[10,], ", ", x[11,], " (mapQ>=5)" )
					 )
	kable(df,align=c("r","r","r","r","r","r"),output=F)
}

##
## VARhelper.GATKug: parse GATK UnifiedGenotyper output for omitted reads
##
VARhelper.GATKug <- function() {

	# log file, which was copied from .bpipe folder
	# contains the runtime STDERR of GATK Unified Genotyper
	LOG <- SHINYREPS_GATKug_LOG
	SUFFIX <- paste0(SHINYREPS_GATKug_SUFFIX, '$')
	
	if(!file.exists(LOG)) {
		print(LOG)
		return("GATK Unified Genotyper statistics not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, full.names=TRUE),function(f) { # list all files and feed them into function one by one
		l <- readLines(f) # read file content to l
		#close(f)
		
		a <- sapply(c(# "reads were filtered out during the traversal out",  #1 # this probably has to be done seperately
				 "failing BadCigarFilter",                            #1
				 "failing BadMateFilter",                             #2
				 "failing DuplicateReadFilter",                       #3
				 "failing FailsVendorQualityCheckFilter",             #4
				 "failing MalformedReadFilter",                       #5
				 "failing MappingQualityUnavailableFilter",           #6
				 "failing NotPrimaryAlignmentFilter",                 #7
				 "failing UnmappedReadFilter"),function(y) {          #8
				 	as.numeric(gsub(".+?(\\d+) reads.+","\\1",l[grep(y,l)])) # grep returns line number, then get the respective line ([]) and extract the first number out of it (gsub and replace the whole line with it)
				 })
		
		l.tmp <- l[grep("reads were filtered out during the traversal out",l)]
		b <- gsub(".+? - (\\d+).+?(\\d+).*","\\1;\\2", l.tmp) #9 #10
		
		return( c(a, as.numeric( strsplit(b, ';')[[1]] )) )	
		
	})
	
	# transform x from list to matrix (in extreme cases also with only one column)
	x <- do.call(cbind, x)
	# set row and column names, and output the md table
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	df <- data.frame("total reads"=x[10,],
					 "total filtered"=paste0( x[9,], " (", round(x[9,] / x[10,] * 100, digits=2), "%)" ),
					 "CIGAR"=paste0( x[1,], " (", round(x[1,] / x[10,] * 100, digits=2), "%)" ),
					 "BadMate"=paste0( x[2,], " (", round(x[2,] / x[10,] * 100, digits=2), "%)" ),
					 "Duplicate"=paste0( x[3,], " (", round(x[3,] / x[10,] * 100, digits=2), "%)" ),
					 "Malformed read"=paste0( x[5,], " (", round(x[5,] / x[10,] * 100, digits=2), "%)" ),
					 "no MappingQuality"=paste0( x[6,], " (", round(x[6,] / x[10,] * 100, digits=2), "%)" ),
					 "not Primary"=paste0( x[7,], " (", round(x[7,] / x[10,] * 100, digits=2), "%)" ),
					 "unmapped"=paste0( x[8,], " (", round(x[8,] / x[10,] * 100, digits=2), "%)" )
					 )
	
	kable(df,align=c("r","r","r","r","r","r","r","r","r"),output=F)

	
	
}

##
## VARhelper.GATKhc: parse GATK HaplotypeCaller output for omitted reads
##
VARhelper.GATKhc <- function() {

	# log file, which was copied from .bpipe folder
	# contains the runtime STDERR of GATK Unified Genotyper
	LOG <- SHINYREPS_GATKhc_LOG
	SUFFIX <- paste0(SHINYREPS_GATKhc_SUFFIX, '$')
	
	if(!file.exists(LOG)) {
		return("GATK Haplotype Caller statistics not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, full.names=TRUE),function(f) { # list all files and feed them into function one by one
		l <- readLines(f) # read file content to l
		#close(f)
		
		a <- sapply(c(# "reads were filtered out during the traversal out",  #1 # this probably has to be done seperately
				 "failing BadCigarFilter",                            #1
				 "failing DuplicateReadFilter",                       #2
				 "failing FailsVendorQualityCheckFilter",             #3
				 "failing HCMappingQualityFilter",                    #4
				 "failing MalformedReadFilter",                       #5
				 "failing MappingQualityUnavailableFilter",           #6
				 "failing NotPrimaryAlignmentFilter",                 #7
				 "failing UnmappedReadFilter"),function(y) {          #8
				 	as.numeric(gsub(".+?(\\d+) reads.+","\\1",l[grep(y,l)])) # grep returns line number, then get the respective line ([]) and extract the first number out of it (gsub and replace the whole line with it)
				 })
		
		l.tmp <- l[grep("reads were filtered out during the traversal out",l)]
		b <- gsub(".+? - (\\d+).+?(\\d+).*","\\1;\\2", l.tmp) #9 #10
		
		return( c(a, as.numeric( strsplit(b, ';')[[1]] )) )	
		
	})
	
	# transform x from list to matrix (in extreme cases also with only one column)
	x <- do.call(cbind, x)
	# set row and column names, and output the md table
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	df <- data.frame("total reads"=x[10,],
					 "total filtered"=paste0( x[9,], " (", round(x[9,] / x[10,] * 100, digits=2), "%)" ),
					 "CIGAR"=paste0( x[1,], " (", round(x[1,] / x[10,] * 100, digits=2), "%)" ),
					 "Duplicate"=paste0( x[2,], " (", round(x[2,] / x[10,] * 100, digits=2), "%)" ),
					 "MappingQuality"=paste0( x[4,], " (", round(x[4,] / x[10,] * 100, digits=2), "%)" ),
					 "Malformed read"=paste0( x[5,], " (", round(x[5,] / x[10,] * 100, digits=2), "%)" ),
					 "no MappingQuality"=paste0( x[6,], " (", round(x[6,] / x[10,] * 100, digits=2), "%)" ),
					 "not Primary"=paste0( x[7,], " (", round(x[7,] / x[10,] * 100, digits=2), "%)" ),
					 "unmapped"=paste0( x[8,], " (", round(x[8,] / x[10,] * 100, digits=2), "%)" )
					 )
	kable(df,align=c("r","r","r","r","r","r","r","r","r"),output=F)
}

##
## VARhelper.GATKvarianteval: parse GATK VariantEvaluation output for variant call statistics
##
VARhelper.GATKvarianteval <- function() {

	# log file, which locates to qc folder
	# contains the output of GATK VariantEval
	LOG <- SHINYREPS_GATKvarianteval
	
	if(!file.exists(LOG)) {
		return("GATK Variant Evaluation statistics not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, full.names=TRUE),function(f) { # list all files and feed them into function one by one
		l <- readLines(f) # read file content to l
		#close(f)
		
		parse.results <- list()
		
		# parse the CompOverlap table containing all variants (SNP & InDel)
		# header: "CompOverlap  CompRod  EvalRod  JexlExpression  Novelty  nEvalVariants  novelSites  nVariantsAtComp  compRate  nConcordant  concordantRate"
		# this yields a vector with one entry
		l.tmp <- l[grep("CompOverlap",l)]
		parse.results[[1]] <- sapply(c(
										"all ",
										"known ",
										"novel " # the additional space is necessary, because otherwise the header line is found as well.
									  ), function(y) {
													#l.grepped <- l.tmp[grep(y,l.tmp)] # select the line of interest
													#regindex <- regexpr("\\d+", l.grepped) # get the first occurence of a number, which is the total number of the respective variants
													as.numeric(gsub(".+?(\\d+).*" ,"\\1", l.tmp[grep(y,l.tmp)]))	# perform the last two lines in one go!
													 }
									 )
		
		# header: CountVariants  CompRod  EvalRod  JexlExpression  Novelty  nProcessedLoci  nCalledLoci  nRefLoci  nVariantLoci  variantRate  variantRatePerBp   nSNPs  nMNPs  nInsertions  nDeletions  nComplex  nSymbolic  nMixed  nNoCalls  nHets  nHomRef  nHomVar  nSingletons  nHomDerived  heterozygosity  heterozygosityPerBp  hetHomRatio  indelRate  indelRatePerBp  insertionDeletionRatio
		# this yields a matrix
		l.tmp <- l[grep("CountVariants",l)]
		parse.results[[2]] <- sapply( c(
										"all ",
										"known ",
										"novel "
									   ), function(y) {
													   # extract nSNPs, nMNPs, nInsertions, nDeletions, nHets, nHomVar, (calculate on our own HetHomRatio & InsDelRatio)
													   as.numeric( unlist( strsplit( l.tmp[grep(y,l.tmp)], "\\s+", perl=T ) )[c(12:15,20,22)] ) # strsplit returns list, which has to be unlisted to make it a vector. Then single elements can be addressed and extracted.
													   #as.numeric(gsub(".+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).*" ,"\\1", l.tmp[grep(y,l.tmp)]))
													  }
									)
		
		l.tmp <- l[grep("MultiallelicSummary",l)]
		parse.results[[3]] <- sapply( c(
										"all ",
										"known ",
										"novel "
									   ), function(y) {
													   # only extract number of MultiAllelicSNP
													   as.numeric( unlist( strsplit( l.tmp[grep(y,l.tmp)], "\\s+", perl=T ) )[8] )
													   #as.numeric(gsub(".+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).*" ,"\\1", l.tmp[grep(y,l.tmp)]))
													  }
									)
		
		l.tmp <- l[grep("TiTvVariantEvaluator",l)]
		parse.results[[4]] <- sapply( c(
										"all ",
										"known ",
										"novel "
									   ), function(y) {
													   # extract nTi and nTv from sample and database
													   as.numeric( unlist( strsplit( l.tmp[grep(y,l.tmp)], "\\s+", perl=T ) )[c(6,7,9,10)] )
													   #as.numeric(gsub(".+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).+?(\\d+).*" ,"\\1", l.tmp[grep(y,l.tmp)]))
													  }
									)
		
		# will consist CompOverlap(all, known, novel)[VariantCount], CountVariants (all, known, novel)[nSNPs, nMNPs, nInsertions, nDeletions, nHets, nHomVar], MultiallelicSummary(all, known, novel)[MultiAllelicSNP], TiTvVariantEvaluator(all, known, novel)[nTi, nTv, TiDB, TvDB]
		output <- c(
					 parse.results[[1]],     # CompOverlap(all, known, novel)[VariantCount] # 1,2,3
					 parse.results[[2]][,1], # CountVariants (all)[nSNPs, nMNPs, nInsertions, nDeletions, nHets, nHomVar] # 4,5,6,7,8,9
					 parse.results[[2]][,2], # CountVariants (known)[nSNPs, nMNPs, nInsertions, nDeletions, nHets, nHomVar] # 10,11,12,13,14,15
					 parse.results[[2]][,3], # CountVariants (novel)[nSNPs, nMNPs, nInsertions, nDeletions, nHets, nHomVar] # 16,17,18,19,20,21
					 parse.results[[3]],     # MultiallelicSummary(all, known, novel)[MultiAllelicSNP] # 22,23,24
					 parse.results[[4]][,1], # TiTvVariantEvaluator(all)[nTi, nTv, TiDB, TvDB # 25,26,27,28
					 parse.results[[4]][,2], # TiTvVariantEvaluator(known)[nTi, nTv, TiDB, TvDB # 29,30,31,32
					 parse.results[[4]][,3]  # TiTvVariantEvaluator(novel)[nTi, nTv, TiDB, TvDB # 33,34,35,36
					 )
		return( output )	
		
	})
	
	
	# transform x from list to matrix (in extreme cases also with only one column)
	x <- do.call(cbind, x)
	
	# set row and column names, and output the md table
	#rownames(x) <- list.files(LOG)
	#rownames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	#rownames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	#
	df <- data.frame("total counts"=paste(x[1,], x[2,],  x[3,], sep=", "),
					 "SNP"=         paste(x[4,] ,x[10,], x[16,],sep=", "),
					 "MNP"=         paste(x[5,], x[11,], x[17,],sep=", "),
					 "Insertion"=   paste(x[6,], x[12,], x[18,],sep=", "),
					 "Deletion"=    paste(x[7,], x[13,], x[19,],sep=", "),
					 #"het hom ratio"=     paste(paste(x[8,] , x[9,],  sep=';'), paste(x[14,] , x[15,], sep=';'), paste(x[20,] , x[21,], sep=';') ,sep=", "),
					 #"ins del ratio"=     paste(paste(x[6,] , x[7,],  sep=';'), paste(x[12,] , x[13,], sep=';'), paste(x[18,] , x[19,], sep=';'), sep=", "),
					 #"Ti Tv ratio"=       paste(paste(x[25,], x[26,], sep=';'), paste(x[29,] , x[30,], sep=';'), paste(x[33,] , x[34,], sep=';'),sep=", ")
					 "het hom ratio"=     paste(round(x[8,] / x[9,], digits=2),   round(x[14,] / x[15,], digits=2), round(x[20,] / x[21,], digits=2),sep=", "),
					 "ins del ratio"=     paste(round(x[6,] / x[7,], digits=2),   round(x[12,] / x[13,], digits=2), round(x[18,] / x[19,], digits=2),sep=", "),
					 "Ti Tv ratio"=       paste(round(x[25,] / x[26,], digits=2), round(x[29,] / x[30,], digits=2), round(x[33,] / x[34,], digits=2),sep=", ")
					 
					 
					 #"CIGAR"=paste0( x[1,], " (", round(x[1,] / x[10,] * 100, digits=2), "%)" ),
					 #"Duplicate"=paste0( x[2,], " (", round(x[2,] / x[10,] * 100, digits=2), "%)" ),
					 #"MappingQuality"=paste0( x[4,], " (", round(x[4,] / x[10,] * 100, digits=2), "%)" ),
					 #"Malformed read"=paste0( x[5,], " (", round(x[5,] / x[10,] * 100, digits=2), "%)" ),
					 #"no MappingQuality"=paste0( x[6,], " (", round(x[6,] / x[10,] * 100, digits=2), "%)" ),
					 #"not Primary"=paste0( x[7,], " (", round(x[7,] / x[10,] * 100, digits=2), "%)" ),
					 #"unmapped"=paste0( x[8,], " (", round(x[8,] / x[10,] * 100, digits=2), "%)" )
					 )
	
	kable(df,align=c("r","r","r","r","r","r","r","r"),output=F)
}


##
## VARhelper.CoveragePlot: produce a plot that is aimed to improve interaction between Genotype, ReadDepth, GenotypeQuality & dbSNP re-ocurrence
##

	# read file
	# vcfData <- read.table(file="results/NA12877.HC.vcf.gz", stringsAsFactors=F)
	# parse
	# Genotype data: unlist(strsplit(vcfData[,10], ":"))[c(1,3,4)]
	# position data: paste(vcfData[,1], vcfData[,2], sep='_')
	# known/novel  : ifelse(vcfData[, 3] == ".", 'novel', 'known')

VARhelper.CoveragePlot <- function() {
	
	#```{r echo=F,results='asis',error=F,warning=F,message=F}
	#cat(VARhelper.CoveragePlot(),sep="\n")
	#```
	
	
	library(ggplot2)
	
	# vcf result file from Haplotype caller
	# need to extract variant properties and compile list
	LOG <- SHINYREPS_RES
	SUFFIX <- paste0(SHINYREPS_RES_GATKhc_SUFFIX, '$')
	
	if(!file.exists(LOG)) {
		return("GATK Haplotype Caller results not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, pattern=SUFFIX, full.names=TRUE),function(f) { # list all files and feed them into function one by one
		l <- read.table(file=f, stringsAsFactors=F, strip.white=T) # read file content to l
		
		# parse
		m <- apply( l, 1, function(l.line){
			
			#trim white spaces that seem to be retained in an apply, but stripped in test cases
			l.line[2] <- gsub("^\\s+|\\s+$", "", l.line[2])
			
			tmp.list <- list()
			tmp.list[[1]] <- paste(l.line[1], l.line[2], sep='_') # chr_position
			tmp.list[[2]] <- ifelse(l.line[3] == ".", 'novel', 'known') # known to dbSNP?
			tmp.list[[3]] <- unlist( strsplit(l.line[10], ":") )[c(1,3,4)] # genotype, read depth & genotype quality
			tmp.list[[4]] <- l.line[1]
			tmp.list[[5]] <- l.line[2]
			
			return( c(tmp.list[[1]], tmp.list[[2]], tmp.list[[3]], tmp.list[[4]], tmp.list[[5]]) ) # return vector
			
			} )
		
		m <- t(m)
		
		return( data.frame("name"  = basename(f),
						   "chr"   = m[,1],
						   "dbSNP" = m[,2],
						   "GT"    = m[,3],
						   "DP"    = as.numeric(m[,4]),
						   "GQ"    = as.numeric(m[,5])
						   )  )
		
	})
	
	# plot
	n <- lapply(x, function(y){
		
		sample.name <- sub(SUFFIX, "", y$name)[1] # y$name is a vector of length data.frame
		p <- ggplot(data=y, aes(GQ, DP, colour=dbSNP, shape=GT) ) + geom_point() + scale_color_manual(values=c("#0000cc","#dd0000")) + scale_shape_manual(values=c(4,1,20,3)) + labs(x="Genotype Quality", y="Read Coverage", title=sample.name)
		print(p)
		
		return() # explicitly return NULL, which will arrive in n
		
		})
	
	return() # return NULL, which will emerge in report environment
	
}

##
## extract tool versions
##

VARhelper.ToolVersions <- function() {
  #tools <- c("FastQC", "STAR", "HTseq", "Subread", "DupRadar", "Samtools", "BedTools", "Picard", "R")
  #variables <- list(SHINYREPS_TOOL_FASTQC, SHINYREPS_TOOL_STAR, SHINYREPS_TOOL_HTSEQ, SHINYREPS_TOOL_SUBREAD, SHINYREPS_TOOL_DUPRADAR, SHINYREPS_TOOL_SAMTOOLS, SHINYREPS_TOOL_BEDTOOLS, SHINYREPS_TOOL_PICARD, SHINYREPS_TOOL_R)
  ## get the last element in path, which is the tool's version (for the tools listed)
  #versions <- sapply(variables, function(x) {
  #  y <- strsplit(x, '/')[[1]]
  #  tail(y, n=1)
  #})
  #
  ## correct the samtools version (second, but last element)
  #tmp_x <- strsplit(SHINYREPS_TOOL_SAMTOOLS, '/')[[1]]
  #versions[6] <- head(tail(tmp_x, n=2), n=1) 
  #
  #df <- data.frame(tool_name=tools, tool_version=versions)
  #
  #kable(df,align=c("l","l"),output=F)
  
}
