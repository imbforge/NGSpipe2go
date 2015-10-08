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
					 different_chromosome=paste0( x[10,], " (", x[11,], " (mapQ>=5))" )
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
## DEhelper.MDS
##
DEhelper.MDS <- function() {
	edgeR::plotMDS.DGEList(y,col=brewer.pal(max(length(levels(group)),3),"Accent")[group])
}

##
## DEhelper.var: variance along log gene count-per-milion
##
DEhelper.var <- function() {
	edgeR::plotBCV(y)
}

##
## DEhelper.cluster: Heatmap of top variant 'n' genes of the counts-per-milion table
##
DEhelper.cluster <- function(n=50) {
	heatmap.2(m[rev(order(v))[1:n],],col=hmcol,trace="none",margin=c(10,6))
}

##
## DEhelper.corr: Heatmap of sample to sample distances
##
DEhelper.corr <- function() {
	heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(13,13))
}

##
## DEhelper.MAplot: MA plots
##
DEhelper.MAplot <- function(i=1,fdr=.05) {
	# get DE genes (p.adjust='BH', pval<.05)
	de <- decideTestsDGE(lrt[[i]],p.value=fdr)
	degenes <- rownames(y)[as.logical(de)]
	
	# MA plot
	plotSmear(lrt[[i]],de.tags=degenes,main=names(lrt)[i])	# MA plot
	abline(h=c(-1,1),col="blue")	# indicate 2-fold changes in the MA plot
	abline(v=0,col="blue")			# indicate >1 counts-per-million
}

##
## DEhelper.DEgenes: show the DE results
##
DEhelper.DEgenes <- function(i=1) {
	ord  <- order(-log(lrt[[i]]$table$FDR),
				   abs(lrt[[i]]$table$logFC),
				  decreasing=TRUE)
	cols <- c("gene","logFC","logCPM","LR","PValue","FDR")
	lrt[[i]]$table[ord,cols]
}




##
## DEhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
DEhelper.Fastqc <- function(web=TRUE) {
	
	# logs folder
	if(!file.exists(SHINYREPS_FASTQC_LOG)) {
		return("Fastqc statistics not available")
	}
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/fastqc" else SHINYREPS_FASTQC_LOG
	
	# construct the image url from the folder contents (skip current dir .)
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
	samples <- list.files(SHINYREPS_DUPRADAR_LOG,pattern="*.png")
	df <- sapply(samples,function(f) {
		paste0("![alt text](",QC,"/",basename(f),")")
	})
	
	# put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
	while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df,"")
	samples <- sapply(df,function(x) {
		x <- sapply(x,function(x) gsub(paste0("^",SHINYREPS_PREFIX),"",basename(x)))
		gsub("_dupRadar.png)","",x)
	})
	df      <- matrix(df     ,ncol=SHINYREPS_PLOTS_COLUMN,byrow=T)
	samples <- matrix(samples,ncol=SHINYREPS_PLOTS_COLUMN,byrow=T)
	
	# add a row with the sample names
	df.names <- matrix(sapply(1:nrow(df),function(i) { c(df[i,],samples[i,]) }),
                       ncol=SHINYREPS_PLOTS_COLUMN,byrow=T)
	colnames(df.names) <- rep(" ",SHINYREPS_PLOTS_COLUMN)
	
	kable(as.data.frame(df.names),output=F)
}

##
## DEhelper.RNAtypes: go through RNAtypes output dir and create a md table with
##     the RNAtypes plots
##
DEhelper.RNAtypes <- function(web=TRUE) {
	
	# logs folder
	if(!file.exists(SHINYREPS_RNATYPES_LOG)) {
		return("RNAtypes statistics not available")
	}
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/RNAtypes" else SHINYREPS_RNATYPES_LOG
	
	# construct the image url from the folder contents (skip current dir .)
	f <- list.files(SHINYREPS_RNATYPES_LOG,pattern="RNAtypes.counts.per.png")
	df <- sapply(f,function(f) {
		paste0("![alt text](",QC,"/",basename(f),")")
	})
	
	# output an md table of 1 columns and 1 row
	df <- matrix(df,ncol=1,nrow=1)
	colnames(df) <- c(" ")

	kable(as.data.frame(df),output=F)
}

##
## DEhelper.geneBodyCov: go through dupRadar output dir and create a md table with
##     the duplication plots
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
	samples <- list.files(SHINYREPS_GENEBODYCOV_LOG,pattern="*.png")
	df <- sapply(samples,function(f) {
		paste0("![alt text](",QC,"/",basename(f),")")
	})
	
	# put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
	while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df,"")
	samples <- sapply(df,function(x) {
		x <- sapply(x,function(x) gsub(paste0("^",SHINYREPS_PREFIX),"",basename(x)))
		gsub(".geneBodyCoverage.curves.png)","",x)
	})
	df      <- matrix(df     ,ncol=SHINYREPS_PLOTS_COLUMN,byrow=T)
	samples <- matrix(samples,ncol=SHINYREPS_PLOTS_COLUMN,byrow=T)
	
	# add a row with the sample names
	df.names <- matrix(sapply(1:nrow(df),function(i) { c(df[i,],samples[i,]) }),
                       ncol=SHINYREPS_PLOTS_COLUMN,byrow=T)
	colnames(df.names) <- rep(" ",SHINYREPS_PLOTS_COLUMN)
	
	kable(as.data.frame(df.names),output=F)
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
	x <- sapply(list.files(FOLDER,pattern=SUFFIX),function(f) {
		
		f <- file(paste0(FOLDER, '/', f))
		l <- readLines(f)
		close(f)
		
		
		sapply(c("Assigned",                      #1
				 "Unassigned_Ambiguity",          #2
				 "Unassigned_MultiMapping",       #3
				 "Unassigned_NoFeatures",         #4
				 "Unassigned_Unmapped",           #5
				 "Unassigned_MappingQuality",     #6
				 "Unassigned_FragementLength",    #7
				 "Unassigned_Chimera",            #8
				 "Unassigned_Secondary",          #9
				 "Unassigned_Nonjunction",        #10
				 "Unassigned_Duplicate"),function(y) {   #11
					as.numeric(  gsub( ".+\t(.+)","\\1",l[grep(y,l)] )  )
				 })	
		
	})
	
	# correct column names
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	
	# create md table (omitting various values that are 0 for now)
	df <- data.frame(assigned=x[1,],
					 unass_ambiguous=x[2,],
					 unass_multimap=x[3,],
					 unass_nofeat=x[4,])
	kable(df,align=c("r","r","r","r"),output=F)
	
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