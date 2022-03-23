##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("knitr")        # for markdown output
library("ggbeeswarm")
library("ChIPpeakAnno")	#for peak venn diagrams
library("GenomicFeatures")
library("RColorBrewer")
library("gridExtra")
library("gtools")
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
    TARGETS <- paste0(SHINYREPS_TARGET)
    if(!file.exists(TARGETS)) {
      return("Targets file not available")
    }
    
    return(read.delim(TARGETS, stringsAsFactors=F))
  }
  
  # read peaks from MACS2 output
  readPeaks <- function(peaksSubdir=subdir, as=peaks_as) {
    
    # check which macs output exist
    if(!file.exists(file.path(SHINYREPS_MACS2,peaksSubdir))) {
      return("MACS2 results not available")
    }
    
    # check is blacklist filtered peak files are available
    if(file.exists(file.path(SHINYREPS_MACS2, peaksSubdir, paste0(targets$IPname, ".vs.", targets$INPUTname,"_macs2_blacklist_filtered_peaks.xls")))[1]) {
      comparisons <- file.path(SHINYREPS_MACS2, peaksSubdir, paste0(targets$IPname, ".vs.", targets$INPUTname,"_macs2_blacklist_filtered_peaks.xls"))
    } else { # no blacklist filtered peak files available, read unfiltered peak files
      comparisons <- file.path(SHINYREPS_MACS2, peaksSubdir, paste0(targets$IPname, ".vs.", targets$INPUTname,"_macs2_peaks.xls"))
    }
    exist <- sapply(comparisons, file.exists) # check if files exist for targets entries
    targets <- targets[exist, ]
    comparisons <- comparisons[exist]
    
    # remove targets which have no peaks
    peakcount <- sapply(comparisons, function(x) {
      tryCatch({
        nrow(read.delim(x, head=TRUE, comment="#"))
      }, error=function(e) 0)
    })
    if(!all(peakcount > 0)) {
      warning("Sample(s) ", paste(basename(comparisons)[!(peakcount > 0)], collapse=", "),
              " excluded because didn't have any peaks called")
      comparisons <- comparisons[peakcount > 0 ] 
      targets <- targets[peakcount > 0, ] 
    }
    
    columnNames2replace <- c(seqnames="chr", abs_summit="summit", pileup="tags", X.log10.pvalue.="-log10 pvalue", X.log10.FDR="-log10 FDR", X.log10.qvalue.="-log10 FDR")
    
    # and return the tables
    peaks <- lapply(comparisons, function(x) {
      if(as=="data.frame") { # read in as data.frame
        x <- tryCatch(read.delim(x, comment.char="#", stringsAsFactors=F), error=function(e) as.data.frame(matrix(ncol=10)))
        colnames(x) <- dplyr::recode(colnames(x), !!!columnNames2replace)
        x[order(x$chr, x$start, x$end),  !colnames(x) %in% c("-log10 pval", "name")]
      } else { # read in as GRanges
        x <- tryCatch({
          x <- ChIPseeker::readPeakFile(x, as = "GRanges")
          if(packageVersion('ChIPseeker')>0) {BiocGenerics::start(x) <- BiocGenerics::start(x)-1} # bug in ChIPseeker: MACS xls files (1-based) are read as 0-based. Modify condition if fixed in future version.
          x
        })
      }
    })
    
    names(peaks) <- paste0(targets$IPname, " vs. ", targets$INPUTname)
    
    return(peaks)
  }
  
  # dispatch tasks
  switch(task, 
         readTargets=readTargets(), 
         readPeaks=readPeaks())
}


##
## ChIPhelper.ComparisonsFromTargets: get the comparisons performed by MACS2 from the targets file
##
ChIPhelper.ComparisonsFromTargets <- function() {
  
  # check for targets.txt and macs2 results
  TARGETS <- paste0(SHINYREPS_TARGET)
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
    tryCatch({
      peaks.ov <- findOverlapsOfPeaks(peak)
        makeVennDiagram(peaks.ov,
                        margin=0.25, cat.default.pos="outer", cat.dist=0.15,
                        cat.fontface=rep("bold", length(peak)),
                        fill=brewer.pal(max(3,length(peak)), "Accent")[1:length(peak)] # brewer.pal requires min n=3
        )
    }, error=function(e) cat("No replicates in", group, " or no peaks were called"))
    cat("\n", fill=T)
  }
}



#'
#' ChIPhelper.UpsetPlot: shows a upset plot for the peaks called in one branch of the pipeline (max 15 peaksets allowed)
#'
#' @param subdir character with sub-directory containing the peak files
#' @param Mode character defining the mode for forming the combination set (one of "distinct", "intersect", "union")
#' @param peakOverlapMode select the value function to calculate size of combination sets ("peaknumber" for number of peaks and/or "bp" for basepairs)
#' @param setsize numeric, maximal number of sets shown
#' @param targetsdf data.frame with targets data. If not NULL, combination sets are highlighted by exclusive sample groups.
#' @param addBarAnnotation logical, whether the intersection sizes are printed on top pf the column annotation
#' @param matrixlist list with upset matrices if calculated externally (allowed elements names: "matrix_peaknumber", "matrix_bp")

ChIPhelper.UpSetPlot <- function(subdir="", Mode = "distinct", peakOverlapMode=c("peaknumber", "bp"), setsize=25, targetsdf=NULL, addBarAnnotation=T, matrixlist=NULL){
  
  #create granges from the peaks
  if(is.null(matrixlist[["peak.ranges"]])) {
    peak.ranges <- ChIPhelper.init("readPeaks", subdir, peaks_as="GRanges")
  } else {
    peak.ranges <- matrixlist[["peak.ranges"]]
  }
  
  if("peaknumber" %in% peakOverlapMode) {
    cat(paste0("#### Overlap of peaks per peak number"), fill=T)
    cat("\n", fill=T)
    #this upset plot is based on the number of peaks which are overlapping (value_fun is length)
    # create upset matrix if not given:
    if(is.null(matrixlist[["matrix_peaknumber"]])) {
      upset_matrix <- make_comb_mat(peak.ranges, mode = Mode, value_fun = length)
      #we subset the matrix to only display the top sets
      upset_matrix <- upset_matrix[order(comb_size(upset_matrix), decreasing = T)[1:setsize]]
    } else {
      upset_matrix <- matrixlist[["matrix_peaknumber"]]
    }
    
    combColors <- fillColors <- "steelblue" # default color
    
    if(!is.null(targetsdf)) { # coloring setnames and respective combinations by group from targets file
      targetsdf <- targetsdf[order(targetsdf$group, targetsdf$IPname), ]
      mypalette <- define.group.palette(length(levels(factor(targetsdf$group))))
      legend_colors <- setNames(mypalette, levels(factor(targetsdf$group)))
      
      setnamesOrderByTargetsdf <- match(targetsdf$IPname, gsub(".vs.*$", "", set_name(upset_matrix)))
      fillColors <- legend_colors[targetsdf$group]
      combColors <- rep("grey40", length.out=length(comb_name(upset_matrix))) 
      for(i in targetsdf$group) { # assign color for all combinations involving a single group
        comb <- paste0(ifelse(targetsdf$group==i, ".", 0), collapse="")
        combColors[grep(comb, comb_name(upset_matrix))] <- legend_colors[i]
      }
    }
    
    ht <- draw(UpSet(upset_matrix,
                     comb_order = order(comb_size(upset_matrix), decreasing = T),
                     comb_col = combColors,
                     bg_col = "#F0F0FF",
                     set_order = if(is.null(targetsdf)) {order(set_size(upset_matrix), decreasing = TRUE)} else {setnamesOrderByTargetsdf},
                     column_title=paste("# of", Mode, "peaks for each set (max", setsize, "sets shown)"),
                     left_annotation = if(is.null(targetsdf)) {NULL} else {rowAnnotation(group = targetsdf$group, col=list(group=legend_colors))}, 
                     right_annotation = upset_right_annotation(upset_matrix, gp = gpar(fill = fillColors)),
                     column_title_gp = gpar(fontsize = 11),
                     row_names_max_width = unit(5, "cm"),
                     row_names_gp = gpar(fontsize = 10)
    ) )
    if(addBarAnnotation){
      od = column_order(ht)
      cs = comb_size(upset_matrix)
      currentvptree <- current.vpTree() # viewport name differs depending on package version
      #print(paste("\ncurrent vptree:", currentvptree))
      interactionsize_viewport <- c("intersection_size", "Intersection\nsize")
      interactionsize_viewport <- interactionsize_viewport[sapply(interactionsize_viewport, grepl, currentvptree)]
      if(length(interactionsize_viewport)==1) {
        decorate_annotation(interactionsize_viewport, {
          grid.text(format(cs[od], scientific=F),
                    x = seq_along(cs), y = unit(cs[od], "native") + unit(20, "pt"), 
                    default.units = "native", just = "bottom", gp = gpar(fontsize = 8), rot=90)        
        })
      } else {
        warning("Viewpoint not found or ambiguous. decorate_annotation is omitted.")
      }
    }
    cat("\n", fill=T)
  }
  
  if("bp" %in% peakOverlapMode) {
    cat("\n", fill=T)
    cat(paste0("#### Overlap of peaks based on bp"), fill=T)
    cat("\n", fill=T)
    cat("\n", fill=T)
    cat("\n", fill=T)
    #this upset plot is based on the number of bp which are overlapping 
    # create upset matrix if not given:
    if(is.null(matrixlist[["matrix_bp"]])) {
      upset_matrix <- make_comb_mat(peak.ranges, mode = Mode)
      #we subset the matrix to only display the top sets
      upset_matrix <- upset_matrix[order(comb_size(upset_matrix), decreasing = T)[1:setsize]]
    } else {
      upset_matrix <- matrixlist[["matrix_bp"]]
    }
    
    if(!is.null(targetsdf)) { # coloring setnames and respective combinations by group from targets file
      targetsdf <- targetsdf[order(targetsdf$group, targetsdf$IPname), ]
      mypalette <- define.group.palette(length(levels(factor(targetsdf$group))))
      legend_colors <- setNames(mypalette, levels(factor(targetsdf$group)))
      
      setnamesOrderByTargetsdf <- match(targetsdf$IPname, gsub(".vs.*$", "", set_name(upset_matrix)))
      fillColors <- legend_colors[targetsdf$group]
      combColors <- rep("grey40", length.out=length(comb_name(upset_matrix))) 
      for(i in targetsdf$group) { # assign color for all combinations involving a single group
        comb <- paste0(ifelse(targetsdf$group==i, ".", 0), collapse="")
        combColors[grep(comb, comb_name(upset_matrix))] <- legend_colors[i]
      }
    }
    
    ht <- draw(UpSet(upset_matrix,
                     comb_order = order(comb_size(upset_matrix), decreasing = T),
                     comb_col = combColors,
                     bg_col = "#F0F0FF",
                     set_order = if(is.null(targetsdf)) {order(set_size(upset_matrix), decreasing = TRUE)} else {setnamesOrderByTargetsdf},
                     column_title=paste("# of", Mode, "bp for each set (max", setsize, "sets shown)"),
                     left_annotation = if(is.null(targetsdf)) {NULL} else {rowAnnotation(group = targetsdf$group, col=list(group=legend_colors))}, 
                     right_annotation = upset_right_annotation(upset_matrix, gp = gpar(fill = fillColors)),
                     column_title_gp = gpar(fontsize = 11),
                     row_names_max_width = unit(5, "cm"),
                     row_names_gp = gpar(fontsize = 10)
    ))
    if(addBarAnnotation){
      od = column_order(ht)
      cs = comb_size(upset_matrix)
      currentvptree <- current.vpTree() # viewport name differs depending on package version
      #print(paste("\ncurrent vptree:", currentvptree))
      interactionsize_viewport <- c("intersection_size", "Intersection\nsize")
      interactionsize_viewport <- interactionsize_viewport[sapply(interactionsize_viewport, grepl, currentvptree)]
      if(length(interactionsize_viewport)==1) {
        decorate_annotation(interactionsize_viewport, {
          grid.text(format(cs[od], scientific=T, digits=2),
                    x = seq_along(cs), y = unit(cs[od], "native") + unit(20, "pt"), 
                    default.units = "native", just = "bottom", gp = gpar(fontsize = 8), rot=90)
        })
      } else {
        warning("Viewpoint not found or ambiguous. decorate_annotation is omitted.")
      }
    }  
    cat("\n", fill=T)
    cat("\n", fill=T)
  }
}


#'
#' ChIPhelper.display.plots: go through folder containing plots and create a md table with the plots
#'
#' @param plotdir character with result directory containing plots
#' @param subdir character with sub-directory for pipeline branch
#' @param plots_column numeric, number of columns to display the plots
#' @param addPlotLabel logical, add plot label from file name
ChIPhelper.display.plots<- function(web=FALSE, plotdir="", subdir="", plots_column=SHINYREPS_PLOTS_COLUMN,
                                    addPlotLabel=T) {
  
  # logs folder
  if(!file.exists(file.path(plotdir, subdir))) {
    return("No plots available")
  }
  
  plots_column <- tryCatch(as.integer(plots_column),error=function(e){4})
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) file.path(paste0("/", basename(plotdir)), subdir) else file.path(plotdir, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(QC, pattern=".png$")
  if(length(samples)==0) {return("No plots available")}
  COLUMNS <- min(length(samples), plots_column)
  df <- sapply(samples, function(f) {
    paste0("![plot img](", QC, "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMNS columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    if(addPlotLabel) {
      x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
      sapply(gsub(".png)$", "", gsub("_", " ", x)), shorten)
    } else {
      x <- sapply(x, function(x) " ")
    }
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="html") %>% kableExtra::kable_styling()
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
  kable(df, align=c("l", "r", "r", "r", "r", "r"), output=F, format="html", row.names=FALSE,
        col.names=c("sample names", "all reads", "mapped (% of all)", "unmapped (% of all)", "too many map. pos. (% all)", "duplicates (% of mapped)")) %>% kableExtra::kable_styling()
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
  colnames(x) <- gsub("\\..*$", "", colnames(x))
  rownames(x)[1] <- if(!isPaired) {"all reads"} else {"all pairs"}
  kable(t(x), align=c(rep("r",10)), output=F, format="html", row.names=T) %>% kableExtra::kable_styling()
}


#' ChIPhelper.ngsReports.Fastqc: joint FastQC report of all samples in the experiment and plot as heatmaps
#' 
#' @param metrics character vector with FastQC plot types to be included. Any combination of "Summary", "BaseQuals", "SeqQuals", "SeqContent", "GcContent", "DupLevels", "Overrep", "AdapterContent".
#'
#' @return list of ggplots made from FastQC data
#' 
ChIPhelper.ngsReports.Fastqc <- function(subdir="", 
                                         metrics=c("BaseQuals", "SeqContent", "GcContent", "DupLevels"), 
                                         ...) {
  
  # output folder
  if(!file.exists(file.path(SHINYREPS_FASTQC, subdir))) {
    return("Fastqc statistics not available")
  }
  
  qclist <- list() # initialize return object
  
  # Loading FastQC Data 
  f <- list.files(file.path(SHINYREPS_FASTQC, subdir), pattern="fastqc.zip$", full.names=TRUE)
  
  # select subset of samples for fastqc figures if desired
  f <- selectSampleSubset(f, ...)
  
  x <- ngsReports::FastqcDataList(f)
  lbls <- gsub(paste0("(^", SHINYREPS_PREFIX, "|.fastqc.zip$)"), "", names(x))
  names(lbls) <- gsub(".fastqc.zip", ".fastq.gz", names(x))
  
  if("Summary" %in% metrics) {
    qclist[["Summary"]] <- ngsReports::plotSummary(x, labels=lbls)
  }
  if("BaseQuals" %in% metrics) {
    qclist[["BaseQuals"]] <- ngsReports::plotBaseQuals(x, labels=lbls)
  }
  if("SeqQuals" %in% metrics) {
    qclist[["SeqQuals"]] <- ngsReports::plotSeqQuals(x, labels=lbls)
  }
  if("SeqContent" %in% metrics) {
    qclist[["SeqContent"]] <- ngsReports::plotSeqContent(x, labels=lbls) +
      theme(legend.position="right") +
      guides(fill=FALSE, color="legend") +
      geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                 data=data.frame(base=c("T", "A", "C", "G")),
                 inherit.aes=FALSE, show.legend=TRUE) +
      scale_color_manual("", values=c("red", "green", "blue", "black"))
  }
  if("GcContent" %in% metrics) {
    qclist[["GcContent"]] <- ngsReports::plotGcContent(x, labels=lbls, theoreticalGC=FALSE)
  }
  if("DupLevels" %in% metrics) {
    qclist[["DupLevels"]] <- ngsReports::plotDupLevels(x, labels=lbls)
  }
  if("Overrep" %in% metrics) {
    qclist[["Overrep"]] <- ngsReports::plotOverrep(x, labels=lbls)
  }
  if("AdapterContent" %in% metrics) {
    qclist[["AdapterContent"]] <- ngsReports::plotAdapterContent(x, labels=lbls)
  }
  
  return(qclist)
}


##
## ChIPhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
ChIPhelper.Fastqc <- function(web=FALSE, subdir="",
                              sampleColumnName =c("IPname", "INPUTname"), 
                              fileColumnName =c("IP", "INPUT")) {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_FASTQC, subdir))) {
    return("Fastqc statistics not available")
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC, subdir)
  
  # construct the image url from the folder ents (skip current dir .)
  samples <- list.dirs(QC, recursive=F)
  samples <- samples[sapply(samples, function(x) {file.exists(file.path(x, "fastqc_data.txt"))})] # exclude potential subdir which is also listed by list.dirs
  
  df <- sapply(samples, function(f) {
    c(paste0("![fastq dup img](", QC, "/", basename(f), "/Images/duplication_levels.png)"), 
      paste0("![fastq qual img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
      paste0("![fastq sequ img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"))
  })
  
  # set row and column names, and output the md table
  df <- as.data.frame(t(df))
  rownames(df) <-  basename(samples)
  colnames(df) <- c("Duplication levels", "Read qualities", "Sequence bias")
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors=F)
    
    if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
      targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
      targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
    } else {
      targets$file <- targets[,fileColumnName]
      targets$sample <- targets[,sampleColumnName]
    }
    
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file 
    # if sample is missing in targets file, use reduced file name
    rownames(df) <- sapply(rownames(df), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                              ifelse(sapply("R1", grepl, i), 
                                                                     paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R1"),
                                                                     ifelse(sapply("R2", grepl, i), 
                                                                            paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R2"),
                                                                            targets[sapply(targets$sample_ext, grepl, i),"sample"])),
                                                              gsub(paste0("^",SHINYREPS_PREFIX),"",i))})                                                    
  } else {
    if(!is.na(SHINYREPS_PREFIX)) {
      rownames(df) <- gsub(paste0("^",SHINYREPS_PREFIX), "", rownames(df))
    }
    rownames(df) <- gsub("_fastqc$", "", rownames(df))
    rownames(df) <- sapply(rownames(df), shorten)
  }
  
  # add a row with the sample name (as given in the rownames) before every row
  df.new <- do.call(rbind,lapply(1:nrow(df),function(i) {rbind(c("",rownames(df)[i],""),df[i,])}))
  rownames(df.new) <- NULL
  kable(df, output=F, align="c", format="html") %>% kableExtra::kable_styling() # print sample names as rownames
}


#' ChIPhelper.Fastqc.custom: prepare customized Fastqc summary plots
#' 
#' @param summarizedPlots logical, if TRUE, data from all samples is summarized in a single plot.
#' @param subdir character with sub-directory to append to the target directory.
#' @param metrics character vector with FastQC plot types to be included. Any combination of "Summary", "BaseQuals", "SeqQuals", "SeqContent", "GcContent", "DupLevels", "Overrep", "AdapterContent".
#' @param sampleColumnName character vector with column names of targets file indicating sample names.
#' @param fileColumnName character vector with column names of targets file indicating sample file names (must have order corresponding to sampleColumnName).
#'
#' @return list of ggplots made from FastQC data
#' 
ChIPhelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, subdir="", 
                                     metrics=c("BaseQuals", "SeqContent", "GcContent", "DupLevels"),
                                     sampleColumnName =c("IPname", "INPUTname"), 
                                     fileColumnName =c("IP", "INPUT")) {
  
  # logs folder
  if(!file.exists(SHINYREPS_FASTQC)) {
    return("Fastqc statistics not available")
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC, subdir)
  
  # read fastqc results in the appropriate format
  f <- list.files(QC, pattern="\\.zip$",full.names=T)
  fastqc.stats <- ngsReports::FastqcDataList(f)
  
  qclist <- list() # initialize return object
  qclist[["no.of.samples"]] <- length(f)
  
  # create proper name vectoir as labels
  lbls <- gsub("_fastqc.zip$", "", names(fastqc.stats))
  names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", names(fastqc.stats))
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors = F)
    
    if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
      targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
      targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
    } else {
      targets$file <- targets[,fileColumnName]
      targets$sample <- targets[,sampleColumnName]
    }
    
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file 
    # if sample is missing in targets file, use reduced file name
    lbls <- sapply(lbls, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                              targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                              gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    
    if(SHINYREPS_PAIRED == "yes") {
      x <- names(lbls)
      lbls <- paste0(lbls, ifelse(grepl("R1", names(lbls)), ".R1", ".R2"))
      names(lbls) <- x
    }
  } else {
    
    if(!is.na(SHINYREPS_PREFIX)) {
      lbls <- gsub(paste0("^",SHINYREPS_PREFIX), "", lbls)
    }
  }
  
  # change names also in fastqc.stats (needed for seq. quality plot)
  names(fastqc.stats) <- lbls
  
  # summary plot (independent from summarizedPlots)
  if("Summary" %in% metrics) {
    qclist[["Summary"]] <- ngsReports::plotSummary(fastqc.stats, labels=lbls)
  }
  
  
  if (summarizedPlots == TRUE) {
    
    # prepare for plotting  
    df <- reshape2::melt(lapply(fastqc.stats , function(x) x@Per_base_sequence_quality[, c("Base","Mean")]))
    names(df)[names(df)=="L1"] <- "samplename"
    
    # color code the samples as done by fastqc:
    # A warning will be issued if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
    # A failure will be raised if the lower quartile for any base is less than 5 or if the median for any base is less than 20. 
    # (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)
    cols <- c(pass    = "#5cb85c",
              warning = "#f0ad4e",
              fail    = "#d9534f")
    colorcode <- do.call(rbind,lapply(names(fastqc.stats),
                                      function(i) {
                                        min.l.quart <- min(fastqc.stats[[i]]@Per_base_sequence_quality$Lower_Quartile)
                                        min.med <- min(fastqc.stats[[i]]@Per_base_sequence_quality$Median)
                                        col.sample <- ifelse(min.l.quart>=10 & min.med>=25, 
                                                             cols["pass"],
                                                             ifelse(min.l.quart>=5 & min.med>=20,
                                                                    cols["warning"],
                                                                    cols["fail"]))
                                        return(data.frame(sample=i,
                                                          min.lower.quart=min.l.quart,
                                                          min.median=min.med,
                                                          col=col.sample))
                                      }
    ))
    
    ## only label "warning"s and "fail"ures 
    to.be.labelled <- colorcode[colorcode$col == cols["warning"] | colorcode$col == cols["fail"],]
    
    ## in case all samples "pass"ed, label "all" in legend
    if (nrow(to.be.labelled) == 0) {
      to.be.labelled <- data.frame(sample="overlay of all samples", col=cols["pass"], row.names=NULL)
    } else {
      to.be.labelled <- rbind(to.be.labelled,c("all other samples","","",col=cols["pass"]))
    }
    
    # fix position on x-axis since there can be intervals of positions summarized in fastqc
    df$position <- factor(df$Base, levels=unique(df$Base))
    xlen <- length(unique(df$Base))
    
    if("BaseQuals" %in% metrics) {
      qclist[["BaseQuals"]] <- ggplot(df, aes(x=as.numeric(position), y=value)) +
        labs(x = "position in read (bp)",
             y = "mean quality score (Phred)") +
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 20), fill = "#edc0c4", alpha = 0.3, color=NA) +
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 28), fill = "#f0e2cc", alpha = 0.3, color=NA) +
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 28, ymax = Inf), fill = "#ceebd1", alpha = 0.3, color=NA) +
        geom_line(aes(color=samplename)) +
        scale_color_manual(values = as.character(colorcode$col),
                           breaks = colorcode$sample,
                           labels = colorcode$sample) +
        geom_point(data=to.be.labelled,
                   mapping=aes(x=NaN, y=NaN, fill=sample),
                   inherit.aes=FALSE, show.legend=TRUE, size = 1.5, shape = 21, color = "white") +
        scale_fill_manual(values = as.character(to.be.labelled$col),
                          breaks = to.be.labelled$sample,
                          labels = to.be.labelled$sample) +
        geom_hline(yintercept = c(0,10,20,30,40),color="white",alpha=0.3) +
        geom_vline(xintercept = seq(0,xlen,10),color="white",alpha=0.3) +
        coord_cartesian(xlim = c(1,xlen), ylim = c(0,42)) +
        guides(color=FALSE,
               fill=guide_legend(title="",ncol=3)) +
        theme(axis.text.x = element_text(size=6,angle=90,hjust=0.5,vjust=0.5),
              axis.text.y = element_text(size=8),
              axis.title  = element_text(size=10),
              plot.title  = element_text(size=12),
              legend.text = element_text(size=7),
              legend.position = "top") +
        scale_x_continuous(breaks=unique(as.numeric(df$position)),
                           labels=unique(df$Base))
    }
    ## use rev(lbls) to plot in reverse order
    if("SeqContent" %in% metrics) {
      qclist[["SeqContent"]] <- ngsReports::plotSeqContent(fastqc.stats, labels=rev(lbls)) +
        labs(y = "") +
        theme(legend.position="right") +
        guides(fill=FALSE, color="legend") +
        geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                   data=data.frame(base=c("T", "A", "C", "G")),
                   inherit.aes=FALSE, show.legend=TRUE) +
        scale_color_manual("", values=c("red", "green", "blue", "black")) 
    }
    
  } else {
    
    if("BaseQuals" %in% metrics) {
      qclist[["BaseQuals"]] <- ngsReports::plotBaseQuals(fastqc.stats, labels=lbls, plotType="boxplot") +
        theme(axis.text.x = element_text(size=5))
    }
    if("SeqContent" %in% metrics) {
      qclist[["SeqContent"]] <- ngsReports::plotSeqContent(fastqc.stats, labels=lbls, plotType="line") +
        theme(axis.text.x = element_text(size=5), legend.position = "top")
    }
  }
  
  # GC content line plot 
  # in case you want to add a theoretical distribution to the plot, use function plotGcContent with 
  # the following settings:
  # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=TRUE, species=SPECIES)
  # the default value for SPECIES is "Hsapiens", thus, if you don't specify it, human will be used as a default
  if("GcContent" %in% metrics) {
    p.gc <- ngsReports::plotGcContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=FALSE) 
    if(!summarizedPlots) {
      p.gc <- p.gc + guides(color=guide_legend(title="",ncol=4)) + 
        theme(legend.position = "top", legend.text = element_text(size=8)) 
    }
    qclist[["GcContent"]] <- p.gc
  }
  
  if("SeqQuals" %in% metrics) {
    qclist[["SeqQuals"]] <- ngsReports::plotSeqQuals(fastqc.stats, usePlotly=summarizedPlots, plotType="line", labels=lbls)
  }
  if("DupLevels" %in% metrics) {
    qclist[["DupLevels"]] <- ngsReports::plotDupLevels(fastqc.stats, usePlotly=summarizedPlots, plotType="line", labels=lbls)
  }
  if("Overrep" %in% metrics) {
    qclist[["Overrep"]] <- ngsReports::plotOverrep(fastqc.stats, usePlotly=summarizedPlots, plotType="line", labels=lbls)
  }
  if("AdapterContent" %in% metrics) {
    qclist[["AdapterContent"]] <- ngsReports::plotAdapterContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", labels=lbls)
  }
  
  return(qclist)
}


##
## ChIPhelper.IPstrength: go through IPstrength output dir and create a md table with
##     the plots
##
ChIPhelper.IPstrength<- function(web=FALSE, subdir="") {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_IPSTRENGTH, subdir))) {
    return("IPstrength statistics not available")
  }
  
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){4})
  if(SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) file.path("/ipstrength", subdir) else file.path(SHINYREPS_IPSTRENGTH, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(QC, pattern=".png$")
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
  
  kable(as.data.frame(df.names), output=F, align="c", format="html") %>% kableExtra::kable_styling()
}

##
## ChIPhelper.peakAnnotation: go through Peak_Annotation output dir and create a md table with 
##      the coverage plots
##
ChIPhelper.peakAnnotationCoverage <- function(web=FALSE, subdir="") {
  # check if peak annotation results are available
  if(!file.exists(file.path(SHINYREPS_PEAK_ANNOTATION, subdir))){
    return("Peak annotation results not available")  
  }
  
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){2})
  if(SHINYREPS_PLOTS_COLUMN < 2 | SHINYREPS_PLOTS_COLUMN > 3) {
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
  
  kable(as.data.frame(df.names), output=F, align="c", format="html") %>% kableExtra::kable_styling()
}

##
## ChIPhelper.peakAnnotationUpSet: go through Peak_Annotation output dir and create a md table with 
##      the UpSet plots
##
ChIPhelper.peakAnnotationUpSet <- function(web=FALSE, subdir="") {
  # check if peak annotation results are available
  if(!file.exists(file.path(SHINYREPS_PEAK_ANNOTATION, subdir))){
    return("Peak annotation results not available")  
  }
  
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){2})
  if(SHINYREPS_PLOTS_COLUMN < 2 | SHINYREPS_PLOTS_COLUMN > 3) {
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
  
  kable(as.data.frame(df.names), output=F, align="c", format="html") %>% kableExtra::kable_styling()
}

##
## ChIPhelper.PhantomPeak: go through PhantomPeak output dir and create a md table with
##     the plots
##
ChIPhelper.PhantomPeak <- function(web=FALSE, subdir="", ...) {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_PHANTOMPEAK, subdir))) {
    return("PhantomPeak statistics not available")
  }
  
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2 | SHINYREPS_PLOTS_COLUMN > 3) {
    SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) file.path("/phantompeak", subdir) else file.path(SHINYREPS_PHANTOMPEAK, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(QC, pattern=".png$")
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
  
  kable(as.data.frame(df.names), output=F, align="c", format="html") %>% kableExtra::kable_styling()
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
  kable(as.data.frame(df), output=F, format="html") %>% kableExtra::kable_styling()
}



##
##ChIPhelper.insertsize: get the insertsize from the qc and display mean and sd 
##
ChIPhelper.insertsize <- function(subdir="", 
                                  sampleColumnName =c("IPname", "INPUTname"), 
                                  fileColumnName =c("IP", "INPUT"), ...){
  
  if (SHINYREPS_PAIRED == "yes") {
    filelist <- list.files(path=file.path(SHINYREPS_INSERTSIZE, subdir), full.names=TRUE, pattern="insertsizemetrics.tsv$")
    filelist <- selectSampleSubset(filelist, ...)
    insertsizes <- lapply(filelist, read.table, sep="\t", header=TRUE, nrow=1)
    insertsizes <- do.call(rbind, insertsizes)
    samplenames <- basename(filelist)
    
    if(file.exists(SHINYREPS_TARGET)){
      
      # get target names
      targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors = F)
      
      if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
        targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
        targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
        targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
        for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
        targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
      } else {
        targets$file <- targets[,fileColumnName]
        targets$sample <- targets[,sampleColumnName]
      }
      
      targets$sample_ext <- gsub("\\..*$", "",targets$file)
      
      # replace files names with nicer sample names given in targets file
      # if sample is missing in targets file, use reduced file name
      samplenames <- sapply(samplenames, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,  
                                                              targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                              gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    } else {
      if(!is.na(SHINYREPS_PREFIX)) {
        samplenames <- gsub(paste0("^",SHINYREPS_PREFIX), "", samplenames)
        samplenames <- gsub("_insertsizemetrics.tsv","", samplenames)
      }
    }
    
    rownames(insertsizes) <- samplenames 
    insertsizes <- insertsizes[,c("MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE", "STANDARD_DEVIATION")]
    colnames(insertsizes) <- c("Median", "Mean", "SD")
    knitr::kable(insertsizes, output=F, align=c("l"), format="html") %>% kableExtra::kable_styling()
  }
}


# Helper to plot the insertsize histogram equivalent to the one from picard
# Input is the Picard generated metrics file
ChIPhelper.insertsize.helper <- function(metricsFile, 
                                         sampleColumnName =c("IPname", "INPUTname"), 
                                         fileColumnName =c("IP", "INPUT")){
  #find the start of our metrics information 
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
  
  title_info <- basename(metricsFile)  
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET, stringsAsFactors = F)
    
    if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
      targets <- targets[, colnames(targets)[colnames(targets) %in% unique(c(fileColumnName, sampleColumnName))]]
      targets <- reshape2::melt(targets, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targets <- targets[!duplicated(targets$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targets$sample[targets$variable == fileColumnName[i]] <- targets[targets$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targets <- targets[, !colnames(targets) %in% sampleColumnName] # sampleColumnName not needed any more
    } else {
      targets$file <- targets[,fileColumnName]
      targets$sample <- targets[,sampleColumnName]
    }
    
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    title_info <- sapply(title_info, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,   
                                                          targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                          gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
  } else {
    if(!is.na(SHINYREPS_PREFIX)) {
      title_info <- gsub(paste0("^",SHINYREPS_PREFIX), "", title_info)
    }
    title_info <- gsub("_insertsizemetrics.tsv$","", title_info)
  }
  
  #we get the histogram which ahs the names of the leves e.g. all_reads and readgroups/sample groups depending on
  #accumulation level which was used.
  #the colnames of histogram are something like all.read.fr_count, all.read.rf_count etc.
  #to get the whole shebang into a wider format we have to add the information
  hist_long <- reshape2::melt(histogram, id.var = "insert_size") %>% 
    extract(col=variable,
            into=c("group", "counttype" ),
            regex='([^\\.]+)\\.([^\\.]+)') %>% 
    dplyr::rename( amount = value)
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
      theme(legend.justification=c(1,1), legend.position=c(1,1), 
            legend.background = element_rect(colour = "transparent", fill = "transparent")) + 
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

  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2 | SHINYREPS_PLOTS_COLUMN > 3) {
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
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
  tryCatch({
    ver <- read.delim(file=SHINYREPS_TOOL_VERSIONS)
    colnames(ver) <- c("Tool name","Environment", "Version")
    knitr::kable(as.data.frame(ver), output=F, align=c("l"), format="html") %>% kableExtra::kable_styling()
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
  kable(samples, output=F, format="html") %>% kableExtra::kable_styling()
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
    resAll <- readRDS(file.path(SHINYREPS_DIFFBIND, subdir, "diffbind.rds"))
    resContrasts <- names(resAll)
    diffbindSettings <- read.table(file.path(SHINYREPS_DIFFBIND, subdir, "diffbind_settings.txt"), header=T, sep="\t", stringsAsFactors = F)
    fdr_threshold  <- diffbindSettings[diffbindSettings$Parameter=="FDR threshold", "Value"]
    fold_threshold <- diffbindSettings[diffbindSettings$Parameter=="Fold threshold", "Value"]
    
    SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){3})
    if(SHINYREPS_PLOTS_COLUMN > 3) {
      SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
    }
    
    cat("\nOverview of applied DiffBind Settings:\n\n")
    cat(knitr::kable(diffbindSettings, format="html") %>% kableExtra::kable_styling(), sep="\n")
    
    
    ### start loop over sub_experiments
    if(any(grepl("^(SubExp_.*)_", resContrasts))) {
      sub_experiments <- unique(gsub("^(SubExp_.*_).*$", "\\1", resContrasts))
    } else {
      sub_experiments <- ""
    }
    
    for (subexpPrefix in sub_experiments) { 
      
      res <- resAll[grep(subexpPrefix, resContrasts)]
      db <- dba.load(dir=file.path(SHINYREPS_DIFFBIND, subdir), file='diffbind', pre=subexpPrefix)
      infodb <- read.table(file.path(SHINYREPS_DIFFBIND, subdir, paste0(subexpPrefix, "info_dba_object.txt")), header=T, sep="\t", stringsAsFactors = F)
      
      if(subexpPrefix !="") {cat(paste0("\n####", gsub("SubExp", "Sub experiment", gsub("_", " ", subexpPrefix)), "\n\n"))}
      cat(paste0("\nDataset in DiffBind analysis (", db$totalMerged, " sites in matrix):\n\n"))
      cat(knitr::kable(infodb, format="html") %>% kableExtra::kable_styling(), sep="\n")
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
                               paste0(subexpPrefix, c("heatmap_occupancy.png", "heatmap_affinity_scores.png", "pca_plot_all_samples.png", 
                                                      "peak_width_boxplot.png", "Overlap_rate_plot.png", "venn_plot_all_contrasts.png")))
      plotsPanel1 <- plotsPanel1[file.exists(plotsPanel1)]
      
      COLUMNS <- min(length(plotsPanel1), SHINYREPS_PLOTS_COLUMN)
      panel1 <- paste0("![diffbind img](", plotsPanel1, ")")
      while(length(panel1) %% COLUMNS != 0) panel1 <- c(panel1, "")
      panel1 <- matrix(panel1, ncol=COLUMNS, byrow=T)
      colnames(panel1) <- rep(" ", COLUMNS)
      cat(kable(as.data.frame(panel1), output=F, align="c", format="html") %>% kableExtra::kable_styling(), sep="\n")
      cat("\n\n")
      
      
      ### for each contrast
      lapply(1:length(res), function(cont) {
        
        maxTableEntries = 20
        if(subexpPrefix !="") { # adjust section hierarchy
          cat("\n\n##### Contrast:", names(res)[cont], "\n")
        } else {
          cat("\n\n#### Contrast:", names(res)[cont], "\n")
        }
        cat("\nIdentified", nrow(res[[cont]]), "differential peaks at FDR", fdr_threshold, 
            "and log fold change threshold", fold_threshold, 
            if(nrow(res[[cont]])>maxTableEntries){paste("(top", maxTableEntries, "entries shown)")}, "\n", fill=TRUE)
        
        if(nrow(res[[cont]])>=1) {
          
          res_table <- res[[cont]][,names(res[[cont]]) %in% c("seqnames",	"start", "end", "width", "strand", "Conc",
                                                              "Fold", "p.value", "FDR", "annotation", "geneId", "distanceToTSS")]
          
          cat(knitr::kable(res_table[1:min(nrow(res[[cont]]), maxTableEntries),], format="html") %>% kableExtra::kable_styling(),sep="\n")
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
          cat(kable(as.data.frame(panel2), output=F, align="c", format="html") %>% kableExtra::kable_styling(), sep="\n")
          cat("\n\n")
          
          
          # plots for 3rd panel (always 3 elements plotted newly)
          cat("\n\n")   
          COLUMNS <- min(3, SHINYREPS_PLOTS_COLUMN)
          numberOfPlots <- 3
          if(!tolower(SHINYREPS_DB) %in% registered_UCSC_genomes()$genome) {numberOfPlots <- numberOfPlots-1} # reduce number of plots if barplot not available
          opar <- par(mfrow=c(ceiling(numberOfPlots/COLUMNS), COLUMNS))
          try(dba.plotMA(db, contrast=cont, cex.main=0.8))  # col.main="white"  
          
          hist(res[[cont]]$Fold, main="", xlab="log fold change", ylab="number of significant peaks", col="grey")
          abline(v=0, lty=2, col="blue")
          
          if(tolower(SHINYREPS_DB) %in% registered_UCSC_genomes()$genome) {
            chromInfo <- getChromInfoFromUCSC(SHINYREPS_DB) # see registered_UCSC_genomes()
            colnames(chromInfo)[colnames(chromInfo)=="length"] <- "size" # in case of older GenomicFeatures versions which give colname "length"
            rownames(chromInfo) <- chromInfo$chrom
            freq <- table(res[[cont]]$seqnames)
            freq <- freq[rev(order(gsub("chr", "", names(freq))))]
            barplot(1e6 * freq/chromInfo[names(freq), "size"], horiz=TRUE, las=1, xlab="number of significant peaks per MB") # count normalized by chrom length
            # barplot(freq, horiz=TRUE, las=1, xlab="number of significant peaks") # absolute sign peaks count instead of normalized count
          }
          
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
    } # end loop subexpPrefix  
  }, error=function(e) cat("Differential binding analysis not available.\n", fill=TRUE))
}



##
## ChIPhelper.cutadapt: get trimming statistics from the Cutadapt folder and display them
## 
#' @param targetsdf targets data.frame or character with file path to targets object
#' @param colorByFactor character with column name of sample table to be used for coloring the plot (e.g. 'group' or 'sample'). Coloring by filename if NULL. 
#' @param fileColumnName character with column name(s) of targets table containing file names
#' @param sampleColumnName character with column name(s) of targets table containing sample names (order must correspond to fileColumnName).
#' @param plotfun define function to be used for plotting
#' @param labelOutliers logical, shall outlier samples be labeled
#' @param outlierIQRfactor numeric, factor is multiplied by IQR to determine outlier
#'
#' @return plot cutadapt statistics as side effect
ChIPhelper.cutadapt <- function(targetsdf=SHINYREPS_TARGET, colorByFactor="group", 
                                sampleColumnName =c("IPname", "INPUTname"), fileColumnName =c("IP", "INPUT"),
                                plotfun=ChIPhelper.cutadapt.plot, labelOutliers=T, outlierIQRfactor=1.5, ...)
{
  
  # logs folder
  if(!all(sapply(SHINYREPS_CUTADAPT_STATS, file.exists))) {
    return(paste("Cutadapt statistics not available for", names(which(!sapply(SHINYREPS_CUTADAPT_STATS, file.exists)))))
  }
  
  x <- list.files(SHINYREPS_CUTADAPT_STATS,pattern='*cutadapt.log$',full.names=TRUE) 
  
  # select subset of samples if desired
  x <- selectSampleSubset(x, ...)
  
  # get Command line parameters of first file
  cutadaptpars <- system(paste("grep \"Command line parameters\"", x[1]), intern=T)
  
  paired <- grepl("(-p )|(--paired-output )", cutadaptpars) # output for R2 available?
  
  x <- sapply(x, function(f) { 
    
    # initialize with NULL in case not needed
    trimmed.R1.perc <- trimmed.R2.perc <- trimmed.reads.perc <- trimmed.qual.perc <- tooshort.reads.perc <- toolong.reads.perc <- NULL 
    
    if(paired) { # log lines slightly differ dependent on se or pe
      total.reads <- system(paste("grep \"Total read pairs processed\"", f, "| awk '{print $5}'"), intern=TRUE)
      total.reads <- gsub(",", "", total.reads)
      trimmed.R1.perc <- system(paste("grep \"Read 1 with adapter\"", f, "| awk '{print $6}'"), intern=TRUE)
      trimmed.R1.perc <- gsub("\\(|\\)|\\%", "", trimmed.R1.perc)
      trimmed.R2.perc <- system(paste("grep \"Read 2 with adapter\"", f, "| awk '{print $6}'"), intern=TRUE)
      trimmed.R2.perc <- gsub("\\(|\\)|\\%", "", trimmed.R2.perc)
      # trimmed.qual.R1 <- system(paste("grep --after-context=2 \"Quality-trimmed\"", f, "| grep \"Read 1\" | awk '{print $3}'"), intern=TRUE) # not needed per read
      # trimmed.qual.R2 <- system(paste("grep --after-context=2 \"Quality-trimmed\"", f, "| grep \"Read 2\" | awk '{print $3}'"), intern=TRUE) # not needed per read
      
    } else {
      total.reads <- system(paste("grep \"Total reads processed\"", f, "| awk '{print $4}'"), intern=TRUE)
      total.reads <- gsub(",", "", total.reads)
      trimmed.reads.perc <- system(paste("grep \"Reads with adapters\"", f, "| awk '{print $5}'"), intern=TRUE)
      trimmed.reads.perc <- gsub("\\(|\\)|\\%", "", trimmed.reads.perc)
    }
    
    trimmed.qual.perc <- system(paste("grep \"Quality-trimmed\"", f, "| awk '{print $4}'"), intern=TRUE)
    trimmed.qual.perc <- gsub("\\(|\\)|\\%", "", trimmed.qual.perc)
    tooshort.reads.perc <- system(paste("grep \"that were too short\"", f, "| awk '{print $7}'"), intern=TRUE)
    tooshort.reads.perc <- gsub("\\(|\\)|\\%", "", tooshort.reads.perc)
    toolong.reads.perc <- system(paste("grep \"that were too long\"", f, "| awk '{print $7}'"), intern=TRUE)
    toolong.reads.perc <- gsub("\\(|\\)|\\%", "", toolong.reads.perc)
    
    # trimming of each adapter
    adapters <- system(paste("grep Sequence:", f, "| awk '{print $9}'"), intern=T)
    adapters.perc <- round(100*(as.numeric(adapters) / as.numeric(total.reads)),1)
    adapterprime <- gsub(";", "", system(paste("grep Sequence:", f, "| awk '{print $5}'"), intern=T))
    
    names(adapters.perc) <- gsub(" *=== *", "", system(paste("grep \"=== .*Adapter\"", f), intern=T))
    namespart1 <- gsub("First read:.*", "R1_", names(adapters.perc))
    namespart1 <- gsub("Second read:.*", "R2_", namespart1)
    namespart2 <- gsub("^.*Adapter", "Adapter", names(adapters.perc))
    names(adapters.perc) <- paste0(if(paired) {namespart1} else {""}, adapterprime, namespart2)
    
    ## add trimmed reads for each adapter here
    return(c("total reads"=total.reads, 
             "bp quality trimmed"=trimmed.qual.perc,
             "R1 adapter trimmed"=trimmed.R1.perc, "R2 adapter trimmed"=trimmed.R2.perc, "adapter trimmed"=trimmed.reads.perc, 
             "too short"=tooshort.reads.perc, "too long"=toolong.reads.perc, adapters.perc))
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
  if (length(indexAdapterSelected)>0) {
    colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)] <- 
      paste0(gsub("Adapter.*$", "", colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)]), cutadaptpars[indexAdapterSelected+1])
  }
  
  #reduce length of file names 
  row.names(x.df) <- basename(colnames(x))
  x.df$filename <- factor(row.names(x.df))
  if(!is.na(SHINYREPS_PREFIX)) {
    row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))
  }
  row.names(x.df) <- gsub("\\.cutadapt\\.log$", "", row.names(x.df))
  if(nrow(x.df)>1){
    if(is.na(SHINYREPS_PREFIX)) {row.names(x.df)  <- gsub(lcPrefix(row.names(x.df) ), "", row.names(x.df) )}
    row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
  }
  
  # passing the different factors given in targetsdf to x.df which was created from cutadapt file names 
  if(!is.null(colorByFactor)) { # add information to x.df
    
    if(is.null(targetsdf)) {stop("If 'colorByFactor' is given you must also provide 'targetsdf'!")}
    
    if(!is.data.frame(targetsdf) && is.character(targetsdf) && file.exists(targetsdf)){
      targetsdf <- read.delim(targetsdf, stringsAsFactors = F)
    } 
    
    if(length(fileColumnName)>1) { # melt targets in case of multiple file name columns (as for ChIP-Seq) and create general targets format 
      targetsdf <- targetsdf[, colnames(targetsdf)[colnames(targetsdf) %in% unique(c(colorByFactor, fileColumnName, sampleColumnName))]]
      targetsdf <- reshape2::melt(targetsdf, measure.vars=fileColumnName, value.name = "file") # 'file' column created
      targetsdf <- targetsdf[!duplicated(targetsdf$file), ] # in case the same inputs are used for several samples
      for(i in 1:length(sampleColumnName)) {targetsdf$sample[targetsdf$variable == fileColumnName[i]] <- targetsdf[targetsdf$variable == fileColumnName[i], sampleColumnName[i]]} # 'sample' column created
      targetsdf <- targetsdf[, !colnames(targetsdf) %in% sampleColumnName] # sampleColumnName not needed any more
      for (i in colorByFactor[!colorByFactor %in% c("file", "sample")]) {targetsdf[, i] <- paste0(targetsdf[, i], " (", targetsdf$variable, ")")} # add fileColumnName to colorByFactor
      targetsdf[,unique(c(colorByFactor, "file", "sample"))] <- lapply(targetsdf[,unique(c(colorByFactor, "file", "sample"))], factor)
    } else {
      targetsdf$file <- targetsdf[,fileColumnName]
      targetsdf$sample <- targetsdf[,sampleColumnName]
    }
    
    targetsdf$file <- gsub("\\..*$", "", targetsdf$file ) # shorten file suffix
    index <- sapply(targetsdf$file, grep, x.df$filename, ignore.case = T) # grep sample name in file names
    if(is.list(index)) {
      targetsdf <- targetsdf[sapply(index, length)!=0,] # remove targetsdf entries not found in x.df
      index <- sapply(targetsdf$file, grep, x.df$filename, ignore.case = T) # redo grep sample name in file names
    }
    
    if(!identical(sort(unname(unlist(index))), 1:nrow(x.df))) {
      stop("There seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
    }
    
    x.df <- data.frame(x.df[unlist(t(index)),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    if("sample" %in% colnames(x.df) && !any(duplicated(x.df$sample))) { # use sample column as identifier if present and unique
      x.df$filename <- x.df$sample
      row.names(x.df) <- x.df$sample } 
    
    if(any(!colorByFactor %in% colnames(x.df))) { # any colorByFactor not available?
      if(all(!colorByFactor %in% colnames(x.df))) { # none colorByFactor available?
        cat("\nNone of the column names given in colorByFactor is available. Perhaps sample names are not part of fastq file names? Using filename instead.")
        colorByFactor <- "filename"
      } else { # one plot each element of colorByFactor
        cat("\n", colorByFactor[!colorByFactor %in% colnames(x.df)], "not available. Using", colorByFactor[colorByFactor %in% colnames(x.df)], "instead.")
        colorByFactor <- colorByFactor[colorByFactor %in% colnames(x.df)]
      }
    }
  } else { # if no colorByFactor defined
    x.df$filename <- row.names(x.df)
    colorByFactor <- "filename"
  }
  
  # melt data frame for plotting
  vars2plot <- c(grep("adapter trimmed", colnames(x.df), value=T), 
                 "too short", "too long",
                 grep("(Adapter)|(})", colnames(x.df), value=T))
  vars2plot <- vars2plot[vars2plot %in% colnames(x.df)]
  x.melt <- reshape2::melt(x.df, measure.vars=vars2plot, variable.name="reads")
  # everything which is not a value should be a factor
  
  # one plot for each element of colorByFactor
  violin.list <- lapply(colorByFactor, plotfun, data=x.melt, labelOutliers=labelOutliers, outlierIQRfactor=outlierIQRfactor) # "colorByFactor" is submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  vars4table <- c(colorByFactor, "total reads", 
                  grep("adapter trimmed", colnames(x.df), value=T),
                  "bp quality trimmed", "too short", "too long",
                  grep("(Adapter)|(})", colnames(x.df), value=T))
  vars4table <- vars4table[vars4table %in% colnames(x.df)]
  vars4table.colnames <- vars4table
  vars4table.colnames[!vars4table.colnames %in% c(colorByFactor, "total reads")] <- paste("%", vars4table.colnames[!vars4table.colnames %in% c(colorByFactor, "total reads")])
  DT::datatable(x.df[,vars4table], 
                options = list(pageLength= 20),
                colnames=vars4table.colnames)
}


# plotting function for ChIPhelper.cutadapt 
ChIPhelper.cutadapt.plot <- function(data, color.value, labelOutliers=T, outlierIQRfactor=1.5){
  
  is_outlier <- function(x) { # function for identification of outlier
    if(IQR(x)!=0) {
      return(x < quantile(x, 0.25) - outlierIQRfactor * IQR(x) | x > quantile(x, 0.75) + outlierIQRfactor * IQR(x))
    } else {
      return(x < mean(x) - outlierIQRfactor * mean(x) | x > mean(x) + outlierIQRfactor * mean(x))
    }
  }
  
  data <- data %>%
    dplyr::group_by(reads) %>%
    dplyr::mutate(outlier=is_outlier(value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outlier=ifelse(outlier,filename,as.numeric(NA))) %>%
    as.data.frame()
  
  ylab <- "% reads"
  
  # prepare palette of appropriate length according to the different factors given in colorByFactor
  colourCount = length(unique(data[,color.value]))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  p <- ggplot(data, aes_string(x="reads",
                               y="value",
                               color=color.value ))+
    geom_quasirandom(groupOnX=TRUE) +
    geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA)  
  if(labelOutliers) {p <- p + ggrepel::geom_text_repel(data=. %>% filter(!is.na(outlier)), aes(label=filename), show.legend=F)}
  p <- p + scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
    ylab(ylab) +
    xlab("") +
    theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1)) 
  
  return(p)
}






##
#' selectSampleSubset: select subset of samples for including in report (e.g. in case of multiple fastq files in scRNA-seq) 
#'
#' @param samples character vector with sample names. 
#' @param samplePattern regular expression to filter (include) on \code{samples}. No filtering if NULL or NA.
#' @param excludePattern regular expression to filter (exclude) on \code{samples}. No filtering if NULL or NA.
#' @param grepInBasename logical. If \code{TRUE} apply pattern to filename, not to full path.
#'
#' @return character vector with selected sample names
selectSampleSubset <- function(samples, samplePattern=NULL, excludePattern=NULL, grepInBasename=T, maxno=NULL) {
  
  if(!gtools::invalid(samplePattern)) {
    x <- if(grepInBasename) {basename(samples)} else {samples} # if TRUE apply pattern to filename, not to full path
    if(!any(grepl(samplePattern, x))) {
      warning(paste("\nYour pattern", samplePattern, "was not found in sample names and is ignored.\n")) 
    } else {
      samples <- samples[grep(samplePattern, x, invert=F)]
    }
    if(length(samples)==0) {stop("\nYou have selected no files!\n")}
  }
  
  if(!gtools::invalid(excludePattern)) {
    x <- if(grepInBasename) {basename(samples)} else {samples} # if TRUE apply pattern to filename, not to full path
    samples <- samples[grep(excludePattern, x, invert=T)]
    if(length(samples)==0) {stop("\nYou have selected no files!\n")}
  }
  
  if(!gtools::invalid(maxno) && is.numeric(maxno)) {
    samples <- samples[1:min(length(samples), maxno)]
    if(maxno > length(samples)) {cat("\nSample number restricted to", maxno)}
  }
  return(samples)
}

#'
#' Define a group color palette.
#'
#' @param num - number of groups
#'
#' @return a color vector
#'
#' @examples group.colors <- define.group.palette(length(levels(colData(dds)[,"group"])))
#'           
define.group.palette <- function(num) {
  
  if (num <=9) {
    group.palette <- brewer.pal(9,"Set1")[1:num]
  } else {
    group.palette <- colorRampPalette(brewer.pal(9,"Set1"))(num)
  }
  return(group.palette)
}

#'
#' Define a replicate color palette.
#'
#' @param num - number of replicates
#'
#' @return a color vector
#'
#' @examples replicate.colors <- define.replicate.palette(length(levels(colData(dds)[,"replicate"])))
#'           
define.replicate.palette <- function(num) {
  
  if (num <=9) {
    replicate.palette <- brewer.pal(8,"Dark2")[1:num]
  } else {
    replicate.palette <- colorRampPalette(brewer.pal(8,"Dark2"))(num)
  }
  return(replicate.palette)
}
