##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
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
## ChIPhelper.init: some time consuming tasks that can be done in advance
##
ChIPhelper.init <- function(task) {
	
	# read targets.txt
	readTargets <- function() {
		TARGETS <- paste0(SHINYREPS_PROJECT,"/",SHINYREPS_TARGETS)
		if(!file.exists(TARGETS)) {
			return("Targets file not available")
		}
	
		return(read.delim(TARGETS))
	}
	
	# read peaks from MACS2 output
	readPeaks <- function() {
	
		# check which macs output exist
		if(!file.exists(SHINYREPS_MACS2)) {
			return("MACS2 results not available")
		}
		comparisons <- paste0(targets$IPname,".vs.",targets$INPUTname,"_macs2_peaks.xls")
		exist <- sapply(paste0(SHINYREPS_MACS2,"/",comparisons),file.exists)
		targets <- targets[exist,]

		# and return the tables
		peaks <- lapply(paste0(SHINYREPS_MACS2,"/",targets$IPname,".vs.",targets$INPUTname,"_macs2_peaks.xls"),function(x) {
			x <- read.delim(x,comment.char="#")
			colnames(x) <- c("chr","start","end","length","summit","tags","-log10 pval","fold enrichment","-log10 FDR","name")
			x[order(x$chr,x$start,x$end),c(-7,-10)]
		})
		names(peaks) <- paste0(targets$IPname," vs. ",targets$INPUTname)

		return(peaks)
	}
	
	# dispatch tasks
	switch(task,
		   readTargets=readTargets(),
		   readPeaks=readPeaks())
}

##
## ChIPhelper.ComparisonsFromTargets: get te comparisons performed by MACS2 from the targets file
##
ChIPhelper.ComparisonsFromTargets <- function() {
	
	# check for targets.txt and macs2 results
	TARGETS <- paste0(SHINYREPS_PROJECT,"/",SHINYREPS_TARGETS)
	if(!file.exists(TARGETS)) {
		return("Targets file not available")
	}

	if(!file.exists(SHINYREPS_MACS2)) {
		return("MACS2 results not available")
	}
	
	# get the comparisons and clean the names
	x <- read.delim(TARGETS)
	comparisons <- paste0(x$IPname,".vs.",x$INPUTname,"_macs2_peaks.xls")
	exist <- sapply(paste0(SHINYREPS_MACS2,"/",comparisons),file.exists)
	comparisons <- gsub(".vs."," vs. ",comparisons)
	comparisons <- gsub("_macs2_peaks.xls","",comparisons)
	
	return(comparisons[exist])
}

##
## ChIPhelper.Peaks: show the peaks called by MACS2
##
ChIPhelper.Peaks <- function(i=1) {
	ord  <- order(peaks[[i]][,"-log10 FDR"],
				  peaks[[i]][,"fold enrichment"],
				  decreasing=TRUE)
	peaks[[i]][ord,]
}

##
## ChIPhelper.BOWTIE: parse bowtie log files and create a md table
##
ChIPhelper.Bowtie <- function() {
	
	# log file
	LOG <- SHINYREPS_BOWTIE_LOG
	if(!file.exists(LOG)) {
		return("Bowtie statistics not available")
	}
	
	# look for the lines aining the strings
	# and get the values associated with this strings
	x <- sapply(list.files(LOG),function(f) {
		
		x <- file(paste0(LOG,"/",f))
		l <- readLines(x)
		close(x)
		
		stats <- sapply(c("reads processed",                             #1
				 "reads with at least one reported alignment",  #2
				 "reads that failed to align",                  #3
				 "reads with alignments suppressed due to -m"), #4
				 function(x) { 
				 	gsub("^.+: ","",l[grep(x,l)])
				 })
	
		# and add the duplicates information
		f <- gsub(".bam.log","_duprm.bam.log",f)
		dups <- if(file.exists(paste0(SHINYREPS_MARKDUPS_LOG,"/",f))) {
			x <- file(paste0(SHINYREPS_MARKDUPS_LOG,"/",f))
			l <- readLines(x)
			close(x)
			gsub(".+Marking (\\d+) records as duplicates.+","\\1",l[grep("Marking \\d+ records as duplicates",l)])
		} else {
			"not available"
		}
		
		c(stats,dups)
	})

	# set row and column names, and output the md table
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(".bam.log$","",colnames(x))
	df <- data.frame(input_reads=x[1,],
					 mapped=x[2,],
					 failed=x[3,],
					 discarded=x[4,],
					 duplicates=paste0(x[5,]," (",round(100 * as.numeric(x[5,]) / as.numeric(x[1,]),2),"%)"))
	kable(df,align=c("r","r","r","r"),output=F)
}

##
## ChIPhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
ChIPhelper.Fastqc <- function(web=TRUE) {
	
	# logs folder
	if(!file.exists(SHINYREPS_FASTQC)) {
		return("Fastqc statistics not available")
	}
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/fastqc" else SHINYREPS_FASTQC
	
	# construct the image url from the folder ents (skip current dir .)
	samples <- list.dirs(SHINYREPS_FASTQC,recursive=F)
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
## ChIPhelper.IPstrength: go through IPstrength output dir and create a md table with
##     the plots
##
ChIPhelper.IPstrength<- function(web=TRUE) {
	
	# logs folder
	if(!file.exists(SHINYREPS_IPSTRENGTH)) {
		return("IPstrength statistics not available")
	}
	
    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
    }
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/ipstrength" else SHINYREPS_IPSTRENGTH
	
	# construct the image url from the folder contents (skip current dir .)
	samples <- list.files(SHINYREPS_IPSTRENGTH,pattern="*.png")
	df <- sapply(samples,function(f) {
		paste0("![alt text](",QC,"/",basename(f),")")
	})
	
	# put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
	while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df,"")
	samples <- sapply(df,function(x) {
		x <- sapply(x,function(x) gsub(paste0("^",SHINYREPS_PREFIX),"",basename(x)))
		gsub("_ipstrength.png)","",x)
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
## ChIPhelper.PhantomPeak: go through PhantomPeak output dir and create a md table with
##     the plots
##
ChIPhelper.PhantomPeak <- function(web=TRUE) {
	
	# logs folder
	if(!file.exists(SHINYREPS_PHANTOMPEAK)) {
		return("PhantomPeak statistics not available")
	}
	
    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
    }
	
	# construct the folder name, which is different for web and noweb
	QC <- if(web) "/phantompeak" else SHINYREPS_PHANTOMPEAK
	
	# construct the image url from the folder contents (skip current dir .)
	samples <- list.files(SHINYREPS_PHANTOMPEAK,pattern="*.png")
	df <- sapply(samples,function(f) {
		paste0("![alt text](",QC,"/",basename(f),")")
	})
	
	# put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
	while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df,"")
	samples <- sapply(df,function(x) {
		x <- sapply(x,function(x) gsub(paste0("^",SHINYREPS_PREFIX),"",basename(x)))
		gsub("_phantompeak.png)","",x)
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
## ChIPhelper.PBC: go through PBC output dir and create a md table with
##     the PBC stats
##
ChIPhelper.PBC <- function() {
	
	# logs folder
	if(!file.exists(SHINYREPS_PBC)) {
		return("PCR bottleneck coefficient statistics not available")
	}
	
	# construct the image url from the folder contents (skip current dir .)
	samples <- list.files(SHINYREPS_PBC,pattern="*.csv")
	df <- sapply(samples,function(f) {
		read.csv(paste0(SHINYREPS_PBC,"/",f))$PBC
	})
	
	# output md table
	df <- as.data.frame(df)
	colnames(df) <- "PBC"
	rownames(df) <- gsub("_PBC.csv","",rownames(df))
	kable(as.data.frame(df),output=F)
}

##
## ChIPhelper.Bustard: call the perl XML interpreter and get the MD output
##
ChIPhelper.Bustard <- function() {
	f  <- SHINYREPS_BUSTARD
	
	if(!file.exists(f)) {
		return("Bustard statistics not available")
	}
	
	# call the perl XSL inetrpreter
	cmd <- paste(" bustard.pl",f)
	try(ret <- system2("perl",cmd,stdout=TRUE,stderr=FALSE))
	
	# check RC
	if(!is.null(attributes(ret))) {
		return(paste("Error parsing bustard statistics. RC:",attributes(ret)$status,"in command: perl",cmd))
	}
	
	ret 	# ret contains already MD code
}
