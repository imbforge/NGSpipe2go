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
## VARhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
VARhelper.Fastqc <- function(web=TRUE) {
	
	# output folder
	if(!file.exists(SHINYREPS_FASTQC_OUT)) {
		return("Fastqc statistics not available")
	}
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/fastqc" else SHINYREPS_FASTQC_OUT
	
	# construct the image url from the folder ents (skip current dir .)
	samples <- list.dirs(SHINYREPS_FASTQC_OUT,recursive=F)
	df <- sapply(samples,function(f) {
		c(paste0("![FastQC dup img](",QC,"/",basename(f),"/Images/duplication_levels.png)"), 
		  paste0("![FastQC qual img](",QC,"/",basename(f),"/Images/per_base_quality.png)"), 
		  paste0("![FastQC seq img](",QC,"/",basename(f),"/Images/per_base_sequence_content.png)"))
	})

	# set row and column names, and output the md table
	df <- as.data.frame(t(df))
	rownames(df) <- gsub(paste0("^",SHINYREPS_PREFIX),"",basename(samples))
	colnames(df) <- c("Duplication","Read qualities","Sequence bias")
	kable(df,output=F, align="c")
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
	x <- lapply(list.files(LOG, pattern=SUFFIX, full.names=TRUE),function(f) { # list all files and feed them into function one by one
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
	colnames(x) <- list.files(LOG, pattern=SUFFIX)
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	df <- data.frame(#"sample name"       = colnames(x),
					 "total reads"       = format( as.numeric(x[1,]),big.mark="," ),
					 "mapped"            = paste0( format( x[3,], big.mark=",") ," (", round(x[3,] / x[1,] * 100, digits=2), "%)" ),
					 "proper pair"       = paste0( format( x[7,], big.mark=",") ," (", round(x[7,] / x[1,] * 100, digits=2), "%)" ) ,
					 "secondary alignments" = paste0( format( x[2,], big.mark=",")," (", round(x[2,] / x[1,] * 100, digits=2), "%)" ),
					 "unmapped"             = paste0( format( x[3,] - x[4,], big.mark=",")," (", round((x[3,] - x[4,]) / x[4,] * 100, digits=2), "%)" ),
					 "different chromosome" = paste0( format( x[10,], big.mark=","), ", ", format( x[11,], big.mark=","), " (mapQ>=5)" )
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
		return("GATK Unified Genotyper statistics not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, pattern=SUFFIX, full.names=TRUE),function(f) { # list all files and feed them into function one by one
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
	colnames(x) <- list.files(LOG, pattern=SUFFIX)
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	df <- data.frame(#"sample name"    = colnames(x),
					 "total reads"    = format( x[10,], big.mark=","),
					 "total filtered" = paste0( format( x[9,], big.mark=","), " (", round(x[9,] / x[10,] * 100, digits=2), "%)" ),
					 "CIGAR"          = paste0( format( x[1,], big.mark=","), " (", round(x[1,] / x[10,] * 100, digits=2), "%)" ),
					 "BadMate"        = paste0( format( x[2,], big.mark=","), " (", round(x[2,] / x[10,] * 100, digits=2), "%)" ),
					 "Duplicate"      = paste0( format( x[3,], big.mark=","), " (", round(x[3,] / x[10,] * 100, digits=2), "%)" ),
					 "Malformed read" = paste0( format( x[5,], big.mark=","), " (", round(x[5,] / x[10,] * 100, digits=2), "%)" ),
					 "no MappingQuality"=paste0( format( x[6,], big.mark=","), " (", round(x[6,] / x[10,] * 100, digits=2), "%)" ),
					 "not Primary"    = paste0( format( x[7,], big.mark=","), " (", round(x[7,] / x[10,] * 100, digits=2), "%)" ),
					 "unmapped"       = paste0( format( x[8,], big.mark=","), " (", round(x[8,] / x[10,] * 100, digits=2), "%)" )
					 )
	rownames(df) <- colnames(x)
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
	x <- lapply(list.files(LOG, pattern=SUFFIX, full.names=TRUE),function(f) { # list all files and feed them into function one by one
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
	colnames(x) <- list.files(LOG, pattern=SUFFIX)
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	
	df <- data.frame(#"sample name"    = colnames(x),
					 "total reads"    = format( x[10,], big.mark=","),
					 "total filtered" = paste0( format( x[9,], big.mark=","), " (", round(x[9,] / x[10,] * 100, digits=2), "%)" ),
					 "CIGAR"          = paste0( format( x[1,], big.mark=","), " (", round(x[1,] / x[10,] * 100, digits=2), "%)" ),
					 "Duplicate"      = paste0( format( x[2,], big.mark=","), " (", round(x[2,] / x[10,] * 100, digits=2), "%)" ),
					 "MappingQuality" = paste0( format( x[4,], big.mark=","), " (", round(x[4,] / x[10,] * 100, digits=2), "%)" ),
					 "Malformed read" = paste0( format( x[5,], big.mark=","), " (", round(x[5,] / x[10,] * 100, digits=2), "%)" ),
					 "no MappingQuality"=paste0( format( x[6,], big.mark=","), " (", round(x[6,] / x[10,] * 100, digits=2), "%)" ),
					 "not Primary"    = paste0( format( x[7,], big.mark=","), " (", round(x[7,] / x[10,] * 100, digits=2), "%)" ),
					 "unmapped"       = paste0( format( x[8,], big.mark=","), " (", round(x[8,] / x[10,] * 100, digits=2), "%)" )
					 )
	rownames(df) <- colnames(x)
	kable(df,align=c("r","r","r","r","r","r","r","r","r"),output=F)
}

##
## VARhelper.GATKvarianteval: parse GATK VariantEvaluation output for variant call statistics
##
VARhelper.GATKvarianteval <- function() {

	# log file, which locates to qc folder
	# contains the output of GATK VariantEval
	LOG <- SHINYREPS_GATKvarianteval
	SUFFIX <- SHINYREPS_GATKvarianteval_SUFFIX
	
	if(!file.exists(LOG)) {
		return("GATK Variant Evaluation statistics not available")
	}
	
	# look for the lines containing the strings
	# and get the values associated with this strings
	# produce a list by lapply to be robust in projects containing only one file
	x <- lapply(list.files(LOG, pattern=SUFFIX, full.names=TRUE),function(f) { # list all files and feed them into function one by one
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
	colnames(x) <- list.files(LOG, pattern=SUFFIX)
	#colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	
	
	df <- data.frame(#"sample name"= colnames(x),format( x[3,], big.mark=",")
					 "total counts"=paste(format( x[1,], big.mark=","), format( x[2,], big.mark=","),  format( x[3,], big.mark=","), sep="; "),
					 "SNP"=         paste(format( x[4,], big.mark=","), format( x[10,], big.mark=","), format( x[16,], big.mark=","),sep="; "),
					 "MNP"=         paste(format( x[5,], big.mark=","), format( x[11,], big.mark=","), format( x[17,], big.mark=","),sep="; "),
					 "Insertion"=   paste(format( x[6,], big.mark=","), format( x[12,], big.mark=","), format( x[18,], big.mark=","),sep="; "),
					 "Deletion"=    paste(format( x[7,], big.mark=","), format( x[13,], big.mark=","), format( x[19,], big.mark=","),sep="; "),
					 "het hom ratio"=     paste(round(x[8,] / x[9,], digits=2),   round(x[14,] / x[15,], digits=2), round(x[20,] / x[21,], digits=2),sep="; "),
					 "ins del ratio"=     paste(round(x[6,] / x[7,], digits=2),   round(x[12,] / x[13,], digits=2), round(x[18,] / x[19,], digits=2),sep="; "),
					 "Ti Tv ratio"=       paste(round(x[25,] / x[26,], digits=2), round(x[29,] / x[30,], digits=2), round(x[33,] / x[34,], digits=2),sep="; ")
					 )
	rownames(df) <- colnames(x)
	kable(df,align=c("r","r","r","r","r","r","r","r"),output=F, )
}


##
## VARhelper.CoveragePlot: produce a plot that is aimed to improve interaction between Genotype, ReadDepth, GenotypeQuality & dbSNP re-ocurrence
##

VARhelper.CoveragePlot <- function() {
	
	# read file
	# vcfData <- read.table(file="results/NA12877.HC.vcf.gz", stringsAsFactors=F)
	# parse
	# Genotype data: unlist(strsplit(vcfData[,10], ":"))[c(1,3,4)]
	# position data: paste(vcfData[,1], vcfData[,2], sep='_')
	# known/novel  : ifelse(vcfData[, 3] == ".", 'novel', 'known')
	
	
	library(ggplot2)
	
	# vcf result file from Haplotype caller
	# need to extract variant properties and compile list
	LOG <- SHINYREPS_RES
	SUFFIX <- paste0(SHINYREPS_GATKhc_SUFFIX, '$')
	
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

Toolhelper.VersionReporter <- function(tool, logfolder) {
	
	LOG <- logfolder
	SUFFIX <- paste0(".log","$")
	
	# logs folder
	if(!file.exists(LOG)) {
		return(paste0(tool, " version not available"))
	}
	
	x <- lapply( list.files(LOG, pattern=SUFFIX, full.names=TRUE), function(f){
		# read all lines
		l <- readLines(f)
		# need to check Version number in one line lower than "VERSION INFO"
		# e.g. FastQC v0.11.3
		l.version <- l[ grep("^VERSION INFO",l) + 1 ]
		
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


#Toolhelper.VersionBWA <- function() {
#	
#	# bwa call without parameters
#	#[...]
#	#Program: bwa (alignment via Burrows-Wheeler transformation)
#	#Version: 0.7.12-r1039
#	#Contact: Heng Li <lh3@sanger.ac.uk>
#	#[...]
#	
#	LOG <- SHINYREPS_BWA_LOG
#	
#	# logs folder
#	if(!file.exists(LOG)) {
#		return("BWA version not available")
#	}
#	
#	x <- lapply( list.files(LOG, full.names=TRUE), function(f){
#		# read all lines
#		l <- readLines(f)
#		# e.g. Version: 0.7.12-r1039
#		l.version <- l[ grep("^VERSION INFO",l) + 1 ]
#		
#		return(l.version)
#		
#		} )
#	
#	# x is a list of always the same content
#	r <- tryCatch(
#		{
#			if (is.null(x[[1]][1])) {
#				return("no version tag")
#			} else {
#				return(x[[1]][1])
#			}
#		},
#		warning = function(w) {
#			return("no version tag")
#		},
#		error = function(e) {
#			return("no version tag")
#		},
#		finally = {}
#	)
#	
#	
#}
#
#
#Toolhelper.VersionFastQC <- function() {
#	
#	# fastqc --version
#	#
#	#FastQC v0.11.3
#	#
#	
#	LOG <- SHINYREPS_FASTQC_LOG
#	SUFFIX <- paste0(".log","$")
#	
#	# logs folder
#	if(!file.exists(LOG)) {
#		return("FastQC version not available")
#	}
#	
#	x <- lapply( list.files(LOG, pattern=SUFFIX, full.names=TRUE), function(f){
#		# read all lines
#		l <- readLines(f)
#		# need to check Version number in one line lower than "VERSION INFO"
#		# e.g. FastQC v0.11.3
#		l.version <- l[ grep("^VERSION INFO",l) + 1 ]
#		
#		return(l.version)
#		
#		} )
#	
#	# x is a list of always the same content
#	r <- tryCatch(
#		{
#			if (is.null(x[[1]][1])) {
#				return("no version tag")
#			} else {
#				return(x[[1]][1])
#			}
#		},
#		warning = function(w) {
#			return("no version tag")
#		},
#		error = function(e) {
#			return("no version tag")
#		},
#		finally = {}
#	)
#	
#	
#}
#
#Toolhelper.VersionMarkDups <- function() {
#	
#	# java -jar MarkDuplicates.jar --version
#	#
#	#1.92(1464)
#	#
#	
#	LOG <- SHINYREPS_MARKDUPS_LOG
#	
#	# logs folder
#	if(!file.exists(LOG)) {
#		return("Picard MarkDuplicates version not available")
#	}
#	
#	x <- lapply( list.files(LOG, full.names=TRUE), function(f){
#		# read all lines
#		l <- readLines(f)
#		# need to check Version number in one line lower than "VERSION INFO"
#		# e.g. FastQC v0.11.3
#		l.version <- l[ grep("^VERSION INFO",l) + 1 ]
#		# extract version number
#		
#		return(l.version)
#		
#		} )
#	
#	# x is a list of always the same content
#	r <- tryCatch(
#		{
#			if (is.null(x[[1]][1])) {
#				return("no version tag")
#			} else {
#				return(x[[1]][1])
#			}
#		},
#		warning = function(w) {
#			return("no version tag")
#		},
#		error = function(e) {
#			return("no version tag")
#		},
#		finally = {}
#	)
#	
#	
#}
#
#Toolhelper.VersionGATKug <- function() {
#	
#	# call GATK and get the version
#	#java -jar GenomeAnalysisTK.jar --version
#	#3.4-46-gbc02625
#	
#	
#	LOG <- SHINYREPS_GATKug_LOG
#	
#	# logs folder
#	if(!file.exists(LOG)) {
#		return("GATK Unified Genotyper version not available")
#	}
#	
#	x <- lapply( list.files(LOG, full.names=TRUE), function(f){
#		# read all lines
#		l <- readLines(f)
#		# need to check Version number in one line lower than "VERSION INFO"
#		# e.g.
#		# VERSION INFO
#		# 3.4-46-gbc02625
#		# \VERSION INFO
#		l.version <- l[ grep("^VERSION INFO",l) + 1 ]
#				
#		return(l.version)
#		
#		} )
#	
#	# x is a list of always the same content
#	r <- tryCatch(
#		{
#			if (is.null(x[[1]][1])) {
#				return("no version tag")
#			} else {
#				return(x[[1]][1])
#			}
#		},
#		warning = function(w) {
#			return("no version tag")
#		},
#		error = function(e) {
#			return("no version tag")
#		},
#		finally = {}
#	)
#	
#	
#}
#
#Toolhelper.VersionGATKhc <- function() {
#	
#	# call GATK and get the version
#	#java -jar GenomeAnalysisTK.jar --version
#	#3.4-46-gbc02625
#	
#	
#	LOG <- SHINYREPS_GATKhc_LOG
#	
#	# logs folder
#	if(!file.exists(LOG)) {
#		return("GATK Haplotype Caller version not available")
#	}
#	
#	x <- lapply( list.files(LOG, full.names=TRUE), function(f){
#		# read all lines
#		l <- readLines(f)
#		# need to check Version number in one line lower than "VERSION INFO"
#		# e.g.
#		# VERSION INFO
#		# 3.4-46-gbc02625
#		# \VERSION INFO
#		l.version <- l[ grep("^VERSION INFO",l) + 1 ]
#				
#		return(l.version)
#		
#		} )
#	
#	# x is a list of always the same content
#	r <- tryCatch(
#		{
#			if (is.null(x[[1]][1])) {
#				return("no version tag")
#			} else {
#				return(x[[1]][1])
#			}
#		},
#		warning = function(w) {
#			return("no version tag")
#		},
#		error = function(e) {
#			return("no version tag")
#		},
#		finally = {}
#	)
#	
#	
#}
#

