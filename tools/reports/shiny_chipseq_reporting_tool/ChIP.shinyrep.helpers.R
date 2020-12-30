##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("knitr")        # for markdown output
library("ggbeeswarm")
library("ChIPpeakAnno")	#for peak venn diagrams
library("RColorBrewer")
library("gridExtra")
library("plyr")
library("dplyr")
library("tidyr")
library("reshape2")
library("ggplot2")
library("ngsReports")
library("DT")
library("DiffBind")
library("ComplexHeatmap")
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
## ChIPhelper.init: some time consuming tasks that can be done in advance
##
ChIPhelper.init <- function(task, subdir="", peaks_as="data.frame") {
  
  # read targets.txt
  readTargets <- function() {
    TARGETS <- paste0(SHINYREPS_PROJECT, "/", SHINYREPS_TARGETS)
    if(!file.exists(TARGETS)) {
      return("Targets file not available")
    }
    
    return(read.delim(TARGETS))
  }
  
  # read peaks from MACS2 output
  readPeaks <- function(peaksSubdir=subdir, as=peaks_as) {
    
    # check which macs output exist
    if(!file.exists(file.path(SHINYREPS_MACS2,peaksSubdir))) {
      return("MACS2 results not available")
    }
    if(file.exists(paste0(targets$IPname, ".vs.", targets$INPUTname,"_macs2_blacklist_filtered_peaks.xls"))[1]) {
      comparisons <- paste0(targets$IPname, ".vs.", targets$INPUTname, "_macs2_blacklist_filtered_peaks.xls")  
    } else {
      comparisons <- paste0(targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls")
    }
    exist <- sapply(paste0(file.path(SHINYREPS_MACS2,peaksSubdir), "/", comparisons), file.exists)
    targets <- targets[exist, ]
    
    columnNames2replace <- c(abs_summit="summit", pileup="tags", X.log10.pvalue.="-log10 pval", fold_enrichment="fold enrichment", X.log10.qvalue.="-log10 FDR")
    
    # and return the tables
    if(file.exists(paste0(file.path(SHINYREPS_MACS2,peaksSubdir), "/", targets$IPname, ".vs.", targets$INPUTname,"_macs2_blacklist_filtered_peaks.xls"))[1]) {
      peaks <- lapply(paste0(file.path(SHINYREPS_MACS2,peaksSubdir), "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_blacklist_filtered_peaks.xls"), function(x) {
        if(as=="data.frame") {
          x <- tryCatch(read.delim(x, comment.char="#"), error=function(e) as.data.frame(matrix(ncol=10)))
          colnames(x) <- plyr::revalue(colnames(x), replace=columnNames2replace , warn_missing = F)
          x[order(x$chr, x$start, x$end),  !colnames(x) %in% c("-log10 pval", "name")]
        } else {
          x <- tryCatch(ChIPseeker::readPeakFile(x))
        }
      })
    } else {
      peaks <- lapply(paste0(file.path(SHINYREPS_MACS2,peaksSubdir), "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls"), function(x) {
        if(as=="data.frame") {
          x <- tryCatch(read.delim(x, comment.char="#"), error=function(e) as.data.frame(matrix(ncol=10)))
          colnames(x) <- plyr::revalue(colnames(x), replace=columnNames2replace , warn_missing = F)
          x[order(x$chr, x$start, x$end),  !colnames(x) %in% c("-log10 pval", "name")]
        } else {
          x <- tryCatch(ChIPseeker::readPeakFile(x))
        }
      })
    }
    
    names(peaks) <- paste0(targets$IPname, " vs. ", targets$INPUTname)
    
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
  TARGETS <- paste0(SHINYREPS_PROJECT, "/", SHINYREPS_TARGETS)
  if(!file.exists(TARGETS)) {
    return("Targets file not available")
  }
  
  if(!file.exists(SHINYREPS_MACS2)) {
    return("MACS2 results not available")
  }
  
  # get the comparisons and clean the names
  x <- read.delim(TARGETS)
  
  if(file.exists(paste0(x$IPname, ".vs.", x$INPUTname,"_macs2_blacklist_filtered_peaks.xls"))[1]) {
    comparisons <- paste0(x$IPname, ".vs.", x$INPUTname, "_macs2_blacklist_filtered_peaks.xls")  
  } else {
    comparisons <- paste0(x$IPname, ".vs.", x$INPUTname, "_macs2_peaks.xls")
  }
  
  exist <- sapply(paste0(SHINYREPS_MACS2, "/", comparisons), file.exists)
  comparisons <- gsub(".vs.", " vs. ", comparisons)
  
  if(file.exists(paste0(x$IPname, ".vs.", x$INPUTname,"_macs2_blacklist_filtered_peaks.xls"))[1]) {
    comparisons <- gsub("_macs2_blacklist_filtered_peaks.xls", "", comparisons)  
  } else {
    comparisons <- gsub("_macs2_peaks.xls", "", comparisons)
  }
  
  return(comparisons[exist])
}

##
## ChIPhelper.Peaks: show the peaks called by MACS2
##
ChIPhelper.Peaks <- function(i=1) {
  ord  <- order(peaks[[i]][, "-log10 FDR"], 
                peaks[[i]][, "fold enrichment"], 
                decreasing=TRUE)
  peaks[[i]][ord, ]
}

##
## ChIPhelper.VennDiagram: shows a venn diagram per group of peaks called
##
ChIPhelper.VennDiagram <- function(subdir=""){
  groups <- unique(targets$group)
  #create granges from the peaks
  peaks <- ChIPhelper.init("readPeaks", subdir)
  peak.ranges <- lapply(peaks, function(x){
    x <- GRanges(seqnames=x$chr,
                 IRanges(x$start,
                         end=x$end),
                 strand="*"
    )
    
  })
  peak.groups <- targets$group
  for(group in groups){
    cat(paste0("#### ", group), fill=T)
    cat("\n", fill=T)
    peak <- peak.ranges[peak.groups==group]
    peaks.ov <- findOverlapsOfPeaks(peak)
    makeVennDiagram(peaks.ov,
                    margin=0.25, cat.default.pos="outer", cat.dist=0.15,
                    cat.fontface=rep("bold", length(peak)),
                    fill=brewer.pal(max(3,length(peak)), "Accent")[1:length(peak)] # brewer.pal requires min n=3
    )
    cat("\n", fill=T)
  }
  
}

##
## ChIPhelper.UpsetPlot: shows a upset plot for the peaks called in one branch 
##
ChIPhelper.UpSetPlot <- function(subdir="", setsize=20){
  #create granges from the peaks
  peaks <- ChIPhelper.init("readPeaks", subdir)
  peak.ranges <- lapply(peaks, function(x){
    x <- GRanges(seqnames=x$chr,
                 IRanges(x$start,
                         end=x$end),
                 strand="*"
    )
    
  })
  
    cat(paste0("#### Overlap of peaks per peak number:"), fill=T)
    cat("\n", fill=T)
    #this upset plot is based on the number of peaks which are overlapping (value_fun is length)
	  upset_matrix <- make_comb_mat(peak.ranges, top_n_sets=setsize, value_fun = length)
	  ht <- draw(UpSet(upset_matrix,
	                   comb_order = order(comb_size(upset_matrix), decreasing = T),
	                   comb_col = "steelblue",
	                   bg_col = "#F0F0FF",
	                   column_title=paste("# of regions for branch", subdir, "\nmax.", setsize, "sets are shown"),
	                   #right_annotation = NULL,
	                   right_annotation = upset_right_annotation(upset_matrix, 
	                                                             gp = gpar(fill = "steelblue"),
	                                                             )
	  ))
	  od = column_order(ht)
	  cs = comb_size(upset_matrix)
	  decorate_annotation("Intersection\nsize", {
	    grid.text(format(cs[od], scientific=T, digits=2),
	              x = seq_along(cs), y = unit(cs[od], "native") + unit(20, "pt"), 
	              default.units = "native", just = "bottom", gp = gpar(fontsize = 8), rot=90)
	  })
    cat("\n", fill=T)
    cat("\n", fill=T)
    cat(paste0("#### Overlap of peaks based on bp:"), fill=T)
    cat("\n", fill=T)
    cat("\n", fill=T)
    cat("\n", fill=T)
    #this upset plot is based on the number of bp which are overlapping 
	  upset_matrix <- make_comb_mat(peak.ranges, top_n_sets=setsize)
	  ht <- draw(UpSet(upset_matrix,
	                   comb_order = order(comb_size(upset_matrix), decreasing = T),
	                   comb_col = "steelblue",
	                   bg_col = "#F0F0FF",
	                   column_title=paste("# overlap in bp for branch", subdir,"\n max.", setsize, "sets are shown"),
	                   #right_annotation = NULL,
	                   right_annotation = upset_right_annotation(upset_matrix, 
	                                                             gp = gpar(fill = "steelblue"),
	                                                             )
	  ))
	  od = column_order(ht)
	  cs = comb_size(upset_matrix)
	  decorate_annotation("Intersection\nsize", {
	    grid.text(format(cs[od], scientific=T, digits=2),
	              x = seq_along(cs), y = unit(cs[od], "native") + unit(20, "pt"), 
	              default.units = "native", just = "bottom", gp = gpar(fontsize = 8), rot=90)
	  })
  
  cat("\n", fill=T)
  cat("\n", fill=T)
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
  
  # look for the lines containing the strings
  # and get the values associated with this strings
  x <- sapply(list.files(LOG), function(f) {
    
    x <- file(paste0(LOG, "/", f))
    l <- readLines(x)
    close(x)
    
    stats <- sapply(c("reads processed",                    #1
                      "reads with at least one reported alignment",  #2
                      "reads that failed to align",                  #3
                      "reads with alignments suppressed due to -m"), #4
                    function(x) { 
                      gsub("^.+: ", "", l[grep(x, l)])
                    })
    
    # and add the duplicates information
    f <- gsub(".bam.log", ".dupmarked.bam.log", f)
    dups <- if(file.exists(paste0(SHINYREPS_MARKDUPS_LOG, "/", f))) {
      x <- file(paste0(SHINYREPS_MARKDUPS_LOG, "/", f))
      l <- readLines(x)
      close(x)
      gsub(".+Marking (\\d+) records as duplicates.+", "\\1", l[grep("Marking \\d+ records as duplicates", l)])
    } else {
      "not available"
    }
    
    stats.return <- c(
      stats[1],                           # reads processed
      unlist(strsplit(stats[2], " "))[1], # reads with at least one reported alignment
      unlist(strsplit(stats[2], " "))[2], # (percentage)
      unlist(strsplit(stats[3], " "))[1], # reads that failed to align 
      unlist(strsplit(stats[3], " "))[2], # (percentage)
      unlist(strsplit(stats[4], " "))[1], # reads with alignments suppressed due to -m
      unlist(strsplit(stats[4], " "))[2]  # (percentage)
    )
    
    
    c(stats.return, dups)
  })
  
  # set row and column names, and output the md table
  colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
  colnames(x) <- gsub(".bam.log$", "", colnames(x))
  df <- data.frame(sample_names=sapply(colnames(x), shorten), 
                   input_reads=format( as.numeric(x[1, ]), big.mark=", "), 
                   mapped=paste( format( as.numeric(x[2, ]), big.mark=", "), x[3, ], sep=" "), 
                   failed=paste( format( as.numeric(x[4, ]), big.mark=", "), x[5, ], sep=" "), 
                   discarded=paste( format( as.numeric(x[6, ]), big.mark=", "), x[7, ], sep=" "), 
                   duplicates=paste0(format(as.numeric(x[8, ]), big.mark=", "), " (", round(100 * as.numeric(x[8, ]) / as.numeric(x[2, ]), 2), "%)")
  )
  kable(df, align=c("l", "r", "r", "r", "r", "r"), output=F, format="markdown", row.names=FALSE,
        col.names=c("sample names", "all reads", "mapped (% of all)", "unmapped (% of all)", "too many map. pos. (% all)", "duplicates (% of mapped)"))
}


##
## ChIPhelper.BOWTIE2: parse bowtie2 paired end or single read mapping log files and create a md table
##
ChIPhelper.Bowtie2 <- function() {
  
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
    
    if(!isPaired) { # modify for single read design
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
  rownames(x)[1] <- if(!isPaired) {"all reads"} else {"all pairs"}
  kable(t(x), align=c(rep("r",10)), output=F, format="markdown", row.names=T)
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
  samples <- list.dirs(SHINYREPS_FASTQC, recursive=F)
  df <- sapply(samples, function(f) {
    c(paste0("![fastq dup img](", QC, "/", basename(f), "/Images/duplication_levels.png)"), 
      paste0("![fastq qual img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
      paste0("![fastq sequ img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"))
  })
  
  # set row and column names, and output the md table
  df <- as.data.frame(t(df))
  x <- gsub(paste0("^", SHINYREPS_PREFIX), "", basename(samples))
  x <- gsub("_fastqc$", "", x)
  rownames(df) <- sapply(x, shorten)
  colnames(df) <- c("Duplication levels", "Read qualities", "Sequence bias")
  kable(df, output=F, align="c", format="markdown")
}

##
## ChIPhelper.ngsReports.Fastqc: joint FastQC report of all samples in the experiment
##
ChIPhelper.ngsReports.Fastqc <- function(...) {
  
  # output folder
  if(!file.exists(SHINYREPS_FASTQC)) {
    return("Fastqc statistics not available")
  }
  
  # Loading FastQC Data 
  f <- list.files(SHINYREPS_FASTQC, pattern="fastqc.zip$", full.names=TRUE)

  # select subset of samples for fastqc figures (e.g. with or without cutadapt suffix)
  f <- selectSampleSubset(f, ...)
  
  x <- ngsReports::FastqcDataList(f)
  lbls <- gsub(paste0("(^", SHINYREPS_PREFIX, "|.fastqc.zip$)"), "", names(x))
  names(lbls) <- gsub(".fastqc.zip", ".fastq.gz", names(x))
  
  print(ngsReports::plotDupLevels(x, labels=lbls))
  print(ngsReports::plotBaseQuals(x, labels=lbls))
  print(ngsReports::plotSeqContent(x, labels=lbls) +
          theme(legend.position="right") +
          guides(fill=FALSE, color="legend") +
          geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                     data=data.frame(base=c("T", "A", "C", "G")),
                     inherit.aes=FALSE, show.legend=TRUE) +
          scale_color_manual("", values=c("red", "green", "blue", "black"))
  )
}


##
## ChIPhelper.IPstrength: go through IPstrength output dir and create a md table with
##     the plots
##
ChIPhelper.IPstrength<- function(web=TRUE, subdir="") {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_IPSTRENGTH, subdir))) {
    return("IPstrength statistics not available")
  }
  
  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) file.path("/ipstrength", subdir) else file.path(SHINYREPS_IPSTRENGTH, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(file.path(SHINYREPS_IPSTRENGTH, subdir), pattern=".png$")
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![IPstrength img](", QC, "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMNS columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_ipstrength.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## ChIPhelper.peakAnnotation: go through Peak_Annotation output dir and create a md table with 
##      the coverage plots
##
ChIPhelper.peakAnnotationCoverage <- function(web=TRUE, subdir="") {
  # check if peak annotation results are available
  if(!file.exists(file.path(SHINYREPS_PEAK_ANNOTATION, subdir))){
    return("Peak annotation results not available")  
  }
  
  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 2L    # default to 2 columns
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(file.path(SHINYREPS_PEAK_ANNOTATION, subdir), pattern="_ChIPseq_Peaks_Coverageplot.png$")
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![Peak_Annotation img](", file.path(SHINYREPS_PEAK_ANNOTATION, subdir), "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMN columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_ChIPseq_Peaks_Coverageplot.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## ChIPhelper.peakAnnotationUpSet: go through Peak_Annotation output dir and create a md table with 
##      the UpSet plots
##
ChIPhelper.peakAnnotationUpSet <- function(web=TRUE, subdir="") {
  # check if peak annotation results are available
  if(!file.exists(file.path(SHINYREPS_PEAK_ANNOTATION, subdir))){
    return("Peak annotation results not available")  
  }
  
  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 2L    # default to 2 columns
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(file.path(SHINYREPS_PEAK_ANNOTATION, subdir), pattern="_ChIPseq_UpSetplot.png$")
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![Peak_Annotation img](", file.path(SHINYREPS_PEAK_ANNOTATION, subdir), "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMN columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_ChIPseq_UpSetplot.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## ChIPhelper.PhantomPeak: go through PhantomPeak output dir and create a md table with
##     the plots
##
ChIPhelper.PhantomPeak <- function(web=TRUE, subdir="", ...) {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_PHANTOMPEAK, subdir))) {
    return("PhantomPeak statistics not available")
  }
  
  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 3L    # default to 4 columns
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) file.path("/phantompeak", subdir) else file.path(SHINYREPS_PHANTOMPEAK, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(file.path(SHINYREPS_PHANTOMPEAK, subdir), pattern=".png$")
  samples <- selectSampleSubset(samples, ...)
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![PhantomPeak img](", QC, "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMN columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_phantompeak.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}



##
## ChIPhelper.PBC: go through PBC output dir and create a md table with
##     the PBC stats
##
ChIPhelper.PBC <- function(subdir="", ...) {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_PBC, subdir))) {
    return("PCR bottleneck coefficient statistics not available")
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(file.path(SHINYREPS_PBC, subdir), pattern="*.csv")
  samples <- selectSampleSubset(samples, ...)
  df <- sapply(samples, function(f) {
    read.csv(file.path(SHINYREPS_PBC, subdir, f))$PBC
  })
  
  # output md table
  df <- as.data.frame(df)
  colnames(df) <- "PBC"
  rownames(df) <- gsub("_PBC.csv", "", rownames(df))
  kable(as.data.frame(df), output=F, format="markdown")
}



##
##ChIPhelper.insertsize: get the insertsize from the qc and display mean and sd 
##
ChIPhelper.insertsize <- function(subdir="", ...){
  
  if (SHINYREPS_PAIRED == "yes") {
    filelist <- list.files(path=file.path(SHINYREPS_INSERTSIZE, subdir), full.names=TRUE, pattern="insertsizemetrics.tsv$")
    filelist <- selectSampleSubset(filelist, ...)
    insertsizes <- lapply(filelist, read.table, sep="\t", header=TRUE, nrow=1)
    insertsizes <- do.call(rbind, insertsizes)
    samplenames <- basename(filelist)
    
    if(length(samplenames)>1) {
      samplenames <- gsub(Biobase::lcPrefix(samplenames), "", samplenames) # remove longest common prefix
      samplenames <- gsub(Biobase::lcSuffix(samplenames), "", samplenames) # remove longest common suffix
    } else {
      if(!is.na(SHINYREPS_PREFIX)) {samplenames <- gsub(SHINYREPS_PREFIX, "", samplenames)}
      samplenames <- gsub("_insertsizemetrics.tsv","", samplenames)
    }
    rownames(insertsizes) <- samplenames 
    insertsizes <- insertsizes[,c("MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE", "STANDARD_DEVIATION")]
    colnames(insertsizes) <- c("Median", "Mean", "SD")
    knitr::kable(insertsizes, output=F, align=c("l"), format="markdown") %>% kableExtra::kable_styling()
  }
}


# Helper to plot the insertsize histogram equivalent to the one from picard
# Input is the Picard generated metrics file
ChIPhelper.insertsize.helper <- function(metricsFile){
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
  
  # title_info
  if(length(metricsFile)>1) {
    title_info <- gsub(Biobase::lcPrefix(basename(metricsFile)), "", basename(metricsFile)) # remove longest common prefix
    title_info <- gsub(Biobase::lcSuffix(title_info), "", title_info) # remove longest common suffix
  } else {
    title_info <- gsub("_insertsizemetrics.tsv$","", basename(metricsFile))
    if(!is.na(SHINYREPS_PREFIX)) {title_info <- gsub(SHINYREPS_PREFIX, "", title_info)}
  }
  
  
  #we get the histogram which ahs the names of the leves e.g. all_reads and readgroups/sample groups depending on
  #accumulation level which was used.
  #the colnames of histogram are something like all.read.fr_count, all.read.rf_count etc.
  #to get the whole shebang into a wider format we have to add the information
  hist_long <- melt(histogram, id.var = "insert_size") %>% 
    separate( variable, 
              sep = "\\.",
              into = c("group", "counttype" )) %>%
    rename( amount = value)
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
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),
                                     error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2) {
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


##
## ChIPhelper.Bustard: call the perl XML interpreter and get the MD output
##
ChIPhelper.Bustard <- function() {
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
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
  tryCatch({
    ver <- read.delim(file=SHINYREPS_TOOL_VERSIONS)
    colnames(ver) <- c("Tool name","Environment", "Version")
    kable(as.data.frame(ver),output=F)
  }, error=function(e) cat("tool versions not available.\n", fill=TRUE))
}

##
## ChIPhelper.BLACKLISTFILTER: go through MACS2 output dir and show BlackList 
##     Filter module file
##
ChIPhelper.BLACKLISTFILTER <- function(subdir="") {
  
  # logs folder
  if(!file.exists(sub("/([^/]*)$", paste0("/",subdir,"/\\1"), SHINYREPS_BLACKLIST_FILTER))) {
    return("Peaks overlapping blacklist regions not available")
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- read.csv(sub("/([^/]*)$", paste0("/",subdir,"/\\1"), SHINYREPS_BLACKLIST_FILTER), row.names=1)
  colnames(samples) <- c("# peaks with MACS2", "# peaks wo overlapping with BRs", "% peaks overlapping BRs")
  kable(samples, output=F, format="markdown")
}


##
## ChIPhelper.GREAT: get the GO enrichment results and display them
##
ChIPhelper.GREAT <- function(){
  
  #csv file
  if(!file.exists(SHINYREPS_GREAT)){
    return("GREAT analysis not available")
  }
}


##
## ChIPhelper.Trackhub: display the UCSC trackhub URL
##
ChIPhelper.Trackhub <- function() {
  
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
## ChIPhelper.diffbind
##
ChIPhelper.diffbind <- function(subdir="") {
  
  tryCatch({
    # read the results (if available
    res <- readRDS(file.path(SHINYREPS_DIFFBIND, subdir, "diffbind.rds"))
    db <- dba.load(dir=file.path(SHINYREPS_DIFFBIND, subdir), file='diffbind', pre="")
    
    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
      SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
    }
    
    cat(paste0("\nDataset in DiffBind analysis (", db$totalMerged, " sites in matrix):\n\n"))
    cat(knitr::kable(data.frame(db$samples[,c("SampleID", "Condition", "Replicate")], FRiP=db$SN), 
                     format="markdown"),sep="\n")
    cat("\n\n")
    
    # plots for 1st panel (independent from contrasts)
    cat("\nCorrelation heatmaps to inspect sample relatedness are generated based on the initial 
    peak caller data (occupancy) as well as based on the affinity scores, which might show a 
    slightly different clustering. The PCA plot given below is based on affinity scores as well.
    The boxplot compares the sizes of consensus peaks with the sizes of the initial sample peaks
    (data of replicates is joined for each group). If the consensus peaks are much larger than the 
    initial sample peaks it should be considered to use intersected peaks instead of union peaks
    as consensus. The overlap rate plot shows the number of peaks that appear in multiple peaksets. 
    These curves typically exhibit a roughly geometric drop-off, with the number of overlapping 
    sites halving as the overlap criterion becomes stricter by one site. When the drop-off is 
    extremely steep, this is an indication that the peaksets do not agree very well. If there are
    replicates you expect to agree, there may be a problem with the experiment.
    Finally, the Venn diagramm shows overlapping differentially bound sites for all analysed contrasts.\n\n")
    
    cat("\n\n")   
    plotsPanel1 <- file.path(SHINYREPS_DIFFBIND, subdir, 
                             c("heatmap_occupancy.png", "heatmap_affinity_scores.png", "pca_plot_all_samples.png", 
                               "peak_width_boxplot.png", "Overlap_rate_plot.png", "venn_plot_all_contrasts.png"))
    plotsPanel1 <- plotsPanel1[file.exists(plotsPanel1)]
    
    COLUMNS <- min(length(plotsPanel1), SHINYREPS_PLOTS_COLUMN)
    panel1 <- paste0("![diffbind img](", plotsPanel1, ")")
    while(length(panel1) %% COLUMNS != 0) panel1 <- c(panel1, "")
    panel1 <- matrix(panel1, ncol=COLUMNS, byrow=T)
    colnames(panel1) <- rep(" ", COLUMNS)
    cat(kable(as.data.frame(panel1), output=F, align="c", format="markdown"), sep="\n")
    cat("\n\n")
    
    
    ### for each contrast
    lapply(1:length(res), function(cont) {
      
      maxTableEntries = 20
      cat("\n\n#### Contrast:", names(res)[cont], "\n", nrow(res[[cont]]), "differential peaks at FDR 5%", 
          if(nrow(res[[cont]])>maxTableEntries){paste("(top", maxTableEntries, "entries shown)")}, "\n", fill=TRUE)
      
      if(nrow(res[[cont]])>=1) {
        
        res_table <- res[[cont]][,names(res[[cont]]) %in% c("seqnames",	"start", "end", "width", "strand", "Conc",
                                                            "Fold", "p.value", "FDR", "annotation", "geneId", "distanceToTSS")]
        
        cat(knitr::kable(res_table[1:min(nrow(res[[cont]]), maxTableEntries),], format="markdown"),sep="\n")
        # DT::datatable(res_table)
        
        cat("\nThe 'Conc' column of the result table shows the mean read concentration over all the samples 
           (the default calculation uses log2 normalized ChIP read counts with control read counts subtracted).
           The 'Fold' column indicates the difference in mean concentrations between the two groups of the contrast,
           with a positive value indicating increased binding affinity in group1 and a negative value indicating 
           increased binding affinity in group2. The additional columns give gene annotation data and confidence 
           measures for identifying these sites as differentially bound, with a raw p-value and a multiple testing 
           corrected FDR.\n")
        
        # plots for 2nd panel 
        cat("\n\n")   
        plotsPanel2 <- file.path(SHINYREPS_DIFFBIND, subdir, paste0(names(res)[cont],
                                                                    c("_pca_plot.png", "_venn_plot.png", "_volcano_plot.png")))
        plotsPanel2 <- plotsPanel2[file.exists(plotsPanel2)]
        
        COLUMNS <- min(length(plotsPanel2), SHINYREPS_PLOTS_COLUMN)
        panel2 <- paste0("![diffbind img](", plotsPanel2, ")")
        while(length(panel2) %% COLUMNS != 0) panel2 <- c(panel2, "")
        panel2 <- matrix(panel2, ncol=COLUMNS, byrow=T)
        colnames(panel2) <- rep(" ", COLUMNS)
        cat(kable(as.data.frame(panel2), output=F, align="c", format="markdown"), sep="\n")
        cat("\n\n")
        
        
        # plots for 3rd panel (always 3 elements plotted newly)
        cat("\n\n")   
        COLUMNS <- min(3, SHINYREPS_PLOTS_COLUMN)
        opar <- par(mfrow=c(ceiling(3/COLUMNS), COLUMNS))
        try(dba.plotMA(db, contrast=cont, cex.main=0.8))  # col.main="white"     
        freq <- table(res[[cont]]$seqnames)
        barplot(freq[rev(order(gsub("chr", "", names(freq))))], horiz=TRUE, las=1, xlab="number of significant peaks")
        hist(res[[cont]]$Fold, main="", xlab="log fold change", ylab="number of significant peaks", col="grey")#,
        abline(v=0, lty=2, col="blue")
        # plot(res[[cont]]$Conc, res[[cont]]$Fold, main="", xlab="log read concentration", ylab="log fold change") 
        # abline(h=0, lty=2, col="blue") # is basically a MA plot showing significant peaks only. Replaced by regular MA plot.
        par(opar)
        
        cat("\nThe PCA plot given here uses only the differentially bound sites of the respective contrast, 
            while the venn diagramm shows overlapping binding sites for the samples of this contrast.
            The dot size of the volcano plots represent the mean read concentration of each site. Significant 
            sites are highlighted in red.
            The MA plot helps to visualize the effect of the data normalization. Points in red represent sites 
            identified as differentially bound.
            The last 2 plots illustrate the distribution of significant peaks across the genome as well as
            the distrbution of fold changes.\n")
        
      }
    })
  }, error=function(e) cat("Differential binding analysis not available.\n", fill=TRUE))
}


##
## ChIPhelper.cutadapt: get cutadapt statistics from the log folder and display them
## 
#' @param targetsdf targets object
#' @param colorByFactor character with column name of sample table to be used for coloring the plot. Coloring by filename if NULL. 
#' @param sampleColumnName character with column name(s) of targets table containing file names
#'
#' @return plot cutadapt statistics as side effect
ChIPhelper.cutadapt <- function(targetsdf=targets, colorByFactor="group", sampleColumnName =c("IP", "INPUT"), ...){
  
  # logs folder
  if(!all(sapply(SHINYREPS_CUTADAPT_LOGS, file.exists))) {
    return(paste("Cutadapt statistics not available for", names(which(!sapply(SHINYREPS_CUTADAPT_LOGS, file.exists)))))
  }
  
  x <- list.files(SHINYREPS_CUTADAPT_LOGS,pattern='*cutadapt.log$',full.names=TRUE) 
  
  # get Command line parameters of first file
  cutadaptpars <- system(paste("grep \"Command line parameters\"", x[1]), intern=T)
  
  paired <- grepl("(-p )|(--paired-output )", cutadaptpars) # output for R2 available?
  
  x <- sapply(x, function(f) { 
    
    trimmed.R1.perc <- trimmed.R2.perc <- trimmed.reads.perc <- tooshort.reads.perc <- NULL # initialise with NULL in case not needed
    
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
    
    tooshort.reads.perc <- system(paste("grep \"that were too short\"", f, "| awk '{print $7}'"), intern=TRUE)
    tooshort.reads.perc <- gsub("\\(|\\)|\\%", "", tooshort.reads.perc)
    
    # trimming of each adapter
    adapters <- system(paste("grep Sequence:", f, "| awk '{print $9}'"), intern=T)
    adapters.perc <- round(100*(as.numeric(adapters) / as.numeric(total.reads)),2)
    adapterprime <- gsub(";", "", system(paste("grep Sequence:", f, "| awk '{print $5}'"), intern=T))
    
    names(adapters.perc) <- gsub(" *=== *", "", system(paste("grep \"=== .*Adapter\"", f), intern=T))
    namespart1 <- gsub("First read:.*", "R1_", names(adapters.perc))
    namespart1 <- gsub("Second read:.*", "R2_", namespart1)
    namespart2 <- gsub("^.*Adapter", "Adapter", names(adapters.perc))
    names(adapters.perc) <- paste0(if(paired) {namespart1} else {""}, adapterprime, namespart2)
    
    ## add trimmed reads for each adapter here
    return(c(total_reads=total.reads, trimmed_R1=trimmed.R1.perc, trimmed_R2=trimmed.R2.perc, 
             trimmed=trimmed.reads.perc, tooshort=tooshort.reads.perc, adapters.perc))
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
  if (length(indexAdapterSelected>0)) {
    colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)] <- 
      paste0(gsub("Adapter.*$", "", colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)]), cutadaptpars[indexAdapterSelected+1])
  }
  
  #reduce length of file names 
  row.names(x.df) <- basename(colnames(x))
  x.df$filename_unmod <- factor(row.names(x.df))
  if(nrow(x.df)>1){
    row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
    row.names(x.df)  <- gsub(lcPrefix(row.names(x.df) ), "", row.names(x.df) )
  }
  
  # passing the different factors given in targetsdf to x.df which was created from cutadapt logfile names 
  if(!is.null(colorByFactor)) { # add information to x.df
    
    if(is.null(targetsdf)) {stop("If 'colorByFactor' is given you must also provide 'targetsdf'!")}
    
    if(length(sampleColumnName>1)) { # melt in case of multiple file name columns (as for ChIP-Seq)
    targetsdf <- targetsdf[,c(colorByFactor, sampleColumnName)]
    targetsdf <- melt(targetsdf, id.vars= colorByFactor, measure.vars=sampleColumnName, value.name = "filename")
    for (i in colorByFactor) {targetsdf[, i] <- paste(targetsdf[, i], targetsdf$variable, sep="_")}
    targetsdf[,c(colorByFactor, "filename")] <- lapply(targetsdf[,c(colorByFactor, "filename")], factor)
    
    } else {
      targetsdf$filename <- targetsdf[,sampleColumnName]
    }
    
    targetsdf$filename <- gsub(lcSuffix(targetsdf$filename ), "", targetsdf$filename ) # shorten filename suffix
    targetsdf$filename <- gsub(lcPrefix(targetsdf$filename ), "", targetsdf$filename ) # shorten filename prefix
    
    index <- sapply(targetsdf$filename, grep, x.df$filename_unmod, ignore.case = T) # grep sample name in file names
    targetsdf <- targetsdf[sapply(index, length) ==1, ] # remove targetsdf entries not (uniquely) found in x.df
    
    if(!identical(sort(unname(unlist(index))), 1:nrow(x.df))) {
      stop("There seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
    }
    
    x.df <- data.frame(x.df[unlist(index),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    rownames(x.df) <- x.df$filename
    
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
  x.melt <- melt(x.df, measure.vars=c(grep("trimmed", colnames(x.df), value=T), 
                                      "tooshort", 
                                      grep("(Adapter)|(})", colnames(x.df), value=T)), variable="reads")
  # everything which is not a value should be a factor
  
  # now we do a violin plot of the trimmed/too_short/etc. ones and color it
  # according to the different factors given in colorByFactor 
  
  # prepare palette of appropriate length
  colourCount = length(unique(x.melt[,colorByFactor]))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  create.violin <- function(x.melt, color.value){
    ylab <- "% reads"
    p <- ggplot(x.melt, aes_string(x="reads",
                                   y="value",
                                   color=color.value ))+
      geom_quasirandom(groupOnX=TRUE) +
      scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
      ylab(ylab) +
      xlab("") +
      #      scale_y_continuous( breaks=seq(0, ceiling(max(x.melt$value)), 10),
      #                          limits = c(0, ceiling(max(x.melt$value)))) + 
      theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1)) 
    
    return(p)
  }
  
  # one plot for each element of colorByFactor
  violin.list <- lapply(colorByFactor, create.violin, x.melt=x.melt) # "colorByFactor" is submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  DT::datatable(x.df[,c("total_reads", 
                        grep("trimmed", colnames(x.df), value=T),
                        "tooshort", 
                        grep("(Adapter)|(})", colnames(x.df), value=T))], 
                options = list(pageLength= 20))
}






##
#' selectSampleSubset: select subset of samples for including in report (e.g. in case of multiple fastq files in scRNA-seq) 
#'
#' @param samples character vector with sample names
#' @param samplePattern regular expression to apply on \code{samples}
#' @param exclude logical indicating if selected samples shall be excluded or included
#' @param grepInBasename logical. If \code{TRUE} apply pattern to filename, not to full path.
#'
#' @return character vector with selected sample names
selectSampleSubset <- function(samples, samplePattern=NULL, exclude=F, grepInBasename=T, maxno=NULL) {
  # use all samples for samplePattern=NULL
  if(!is.null(samplePattern)) {
    
    x <- if(grepInBasename) {basename(samples)} else {samples} # if TRUE apply pattern to filename, not to full path
    samples <- samples[grep(samplePattern, x, invert=exclude)]
    
    if(length(samples)==0) {stop("\nYou have selected no files!\n")}
  }
  if(!is.null(maxno)) {
    samples <- samples[1:min(length(samples), maxno)]
    if(maxno < length(samples)) {cat("\nSample number restricted to", maxno)}
  }
  return(samples)
}

