##################################
##
## helper functions to create the plots for the Shiny report
##
##################################


#' attach_packages: load required packages
#'
#' @param pkg character vector with package names
#'
#' @return character vector with package names which have not yet been attached to workspace before yet. 
#' This vector may be used to detach packages not needed any more.

attach_packages <- function(pkg) {
  not.yet.attached.pkg <- pkg[!(pkg %in% loadedNamespaces())]
  for (p in pkg) { 
    require(p, character.only = TRUE) 
  }
  return(unique(not.yet.attached.pkg))
}



#' loadGlobalVars: read configuration from bpipe vars
#'
#' @param f - character containing a file name defining multiple variables for reporting to run. 
#' When f is a character vector with multiple filenames, multiple (non-unique) entries like path variables
#' are assigned as vector. This is useful when do a meta analysis using target files from different project
#' folders. Most helper functions can use these vector variables to read the respective target files
#' from multiple designations.  
#'
#' @return A set of variables that are mentioned in input file, e.g.
#'         SHINYREPS_PROJECT <- "projects/example_project/"
#'              
#' @description File content should be:
#'              SHINYREPS_PROJECT=projects/example_project
#'              
loadGlobalVars <- function(f="shinyReports.txt") {
  
  conf_list <- lapply(f, function(f) {
    # read in the conf file(s)
    conf <- readLines(f)
    conf <- conf[grep("^SHINYREPS_", conf)]
    conf <-  sapply(conf, function(x) {
      x <- unlist(strsplit(x, "=", fixed=T))
      if(length(x) == 2)
        x
      else
        c(x[1], NA)
    })
    colnames(conf) <- conf[1,]
    conf <- conf[-1,, drop=F]
  })
  
  # combine conf files if more than one
  #conf <- plyr::rbind.fill.matrix(conf_list)
  conf <- dplyr::bind_rows(lapply(conf_list, as.data.frame, drop=F, stringsAsFactors =F))
  for(i in colnames(conf)) {
    assign(i, unique(conf[,i]), envir=.GlobalEnv)
  }
  
  invisible(0)
}

##
## Some generic functions
##

#'
#' shorten: if a text string is longer than certain length, shorten by showing the first and last characters
#'
#' @param x - a character vector
#' @param max.len - a numeric vector of length 1 indicating the maximal length of a string to be tolerated
#' @param ini - a numeric vector of length 1 indicating how many characters at the beginning of the string should be conserved 
#' @param end - a numeric vector of length 1 indicating how many characters at the end of the string should be conserved 
#'
#' @return a character vector
#'
#' @examples long_string <- "This_is_a_very_long_string_and_should_be_shortened"
#'           short_string <- shorten(long_string, max.len = 20, ini = 5, end = 9)
#'           print(short_string)
#'           ## [1] "This_..._shortened"
#'           
shorten <- function(x, max.len=40, ini=20, end=15) {
  l <- nchar(x)
  if(l > max.len) paste(substr(x, 1, ini), substr(x, (l-end), l), sep="...") else x
}



##
## MPShelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
MPShelper.Fastqc <- function(web=FALSE, subdir="", ...) {
  
  # logs folder
  if(!all(sapply(file.path(SHINYREPS_FASTQC_OUT, subdir), file.exists))) {
    return(paste("Fastqc statistics not available for", names(which(!sapply(file.path(SHINYREPS_FASTQC_OUT, subdir), file.exists)))))
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.dirs(QC, recursive=F, full.names = T)
  samples <- samples[sapply(samples, function(x) {file.exists(file.path(x, "fastqc_data.txt"))})] # exclude potential subdir which is also listed by list.dirs
  
  # select subset of samples for fastqc figures (e.g. merged single cell pools) or use all samples for samplePattern=NULL
  samples <- selectSampleSubset(samples, ...)

  df <- sapply(samples, function(f) {
    c(paste0("![fastqc img](", f, "/Images/per_base_quality.png)"), 
      paste0("![fastqc img](", f, "/Images/per_base_sequence_content.png)"),
      paste0("![fastqc img](", f, "/Images/per_sequence_gc_content.png)"))
  })
  
  # set row and column names, and output the md table
  df <- as.data.frame(t(df))
  rownames(df) <-  basename(samples)
  colnames(df) <- c("Read qualities", "Sequence bias", "GC content")
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET)
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
# kable(df.new, output=F, align="c", format="markdown") # print sample names in additional rows
  kable(df, output=F, align="c") # print sample names as rownames
}


##
## MPShelper.ngsReports.Fastqc: joint FastQC report of all samples in the experiment
##
MPShelper.ngsReports.Fastqc <- function(subdir="", ...) {
  
  # output folder
  if(!file.exists(file.path(SHINYREPS_FASTQC_OUT, subdir))) {
    return("Fastqc statistics not available")
  }
  
  # Loading FastQC Data 
  f <- list.files(file.path(SHINYREPS_FASTQC_OUT, subdir), pattern="fastqc.zip$", full.names=TRUE)
  
  # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
  f <- selectSampleSubset(f, ...)
  
  x <- ngsReports::FastqcDataList(f)
  lbls <- gsub(paste0("(^", SHINYREPS_PREFIX, "|.fastqc.zip$)"), "", names(x))
  names(lbls) <- gsub(".fastqc.zip", ".fastq.gz", names(x))
  
  print(ngsReports::plotBaseQuals(x, labels=lbls))
  print(ngsReports::plotSeqContent(x, labels=lbls) +
          theme(legend.position="right") +
          guides(fill=FALSE, color="legend") +
          geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                     data=data.frame(base=c("T", "A", "C", "G")),
                     inherit.aes=FALSE, show.legend=TRUE) +
          scale_color_manual("", values=c("red", "green", "blue", "black"))
  )
  print(ngsReports::plotGcContent(x, plotType="line", gcType="Genome", theoreticalGC = F, labels=lbls))  
}

##
## MPShelper.Fastqc.custom: prepare Fastqc summary plots
##
MPShelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, subdir="", ...) {
  
  # logs folder
  if(!file.exists(file.path(SHINYREPS_FASTQC_OUT, subdir))) {
    return("Fastqc statistics not available")
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
  # read fastqc results in the appropriate format
  f <- list.files(QC, pattern="\\.zip$",full.names=T)
  
  # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
  f <- selectSampleSubset(f, ...)
  
  fastqc.stats <- ngsReports::FastqcDataList(f)
  
  # create proper name vector as labels
  lbls <- gsub("_fastqc.zip$", "", names(fastqc.stats))
  names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", names(fastqc.stats))
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET)
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file 
    # if sample is missing in targets file, use reduced file name
    lbls <- sapply(lbls, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                              targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                              gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
    
    if(SHINYREPS_PAIRED == "yes") {
      x <- names(lbls)
      lbls <- paste0(lbls, ifelse(grepl("R1", names(lbls)), "_R1", "_R2"))
      names(lbls) <- x
    }
  } else {
    
    if(!is.na(SHINYREPS_PREFIX)) {
      lbls <- gsub(paste0("^",SHINYREPS_PREFIX), "", lbls)
    }
  }
  
  # change names also in fastqc.stats (needed for seq. quality plot)
  names(fastqc.stats) <- lbls
  
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
    
    p.qual <- ggplot(df, aes(x=as.numeric(position), y=value)) +
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
    
    ## use rev(lbls) to plot in reverse order
    p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=rev(lbls)) +
      labs(y = "") +
      theme(legend.position="right") +
      guides(fill=FALSE, color="legend") +
      geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                 data=data.frame(base=c("T", "A", "C", "G")),
                 inherit.aes=FALSE, show.legend=TRUE) +
      scale_color_manual("", values=c("red", "green", "blue", "black")) 
    
  } else {
    
    p.qual <- ngsReports::plotBaseQuals(fastqc.stats, labels=lbls, plotType="boxplot") +
      theme(axis.text.x = element_text(size=5))
    p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=lbls, plotType="line") +
      theme(axis.text.x = element_text(size=5),
            legend.position = "top")
  }
  
  # GC content line plot 
  # in case you want to add a theoretical distribution to the plot, use function plotGcContent with 
  # the following settings:
  # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=TRUE, species=SPECIES)
  # the default value for SPECIES is "Hsapiens", thus, if you don't specify it, human will be used as a default
  p.gc <- ngsReports::plotGcContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=FALSE) 
  if(!summarizedPlots) {
    p.gc <- p.gc + guides(color=guide_legend(title="",ncol=4)) + 
      theme(legend.position = "top", legend.text = element_text(size=8)) 
  }
  return(list(no.of.samples=length(f), p.qual=p.qual, p.content=p.content, p.gc=p.gc))
}



#' MPShelper.fastqscreen: summarizes FastQScreen results, creates summarized barplots, only relevant contanimants shown
#'
#' @param perc.to.plot - a numeric vector of length 1 setting the percent cutoff of relevant contaminants, if any sample
#'                       shows more than perc.to.plot, contaminant will be shown in plot
#'
#' @return a list including a plot, the number of samples, and the number of plotted contaminants
#'
MPShelper.fastqscreen <- function(subdir="", perc.to.plot = 1, ncol=2, ...) {
  
  # logs folder
  if(!file.exists(SHINYREPS_FASTQSCREEN_OUT)) {
    return(list(errortext="FastQScreen statistics not available",
                no.of.genomes=1,
                no.of.samples=1,
                no.of.rows=1))
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- file.path(SHINYREPS_FASTQSCREEN_OUT, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.dirs(QC, recursive=F, full.names = T)
  samples <- samples[sapply(samples, function(x) {file.exists(file.path(x, "fastqscreen.conf"))})] # exclude potential subdir which is also listed by list.dirs or recursive list.files
  samples <- list.files(samples, pattern="_screen.txt$", recursive=F, full.names=T)
  
  # select subset of samples to plot
  samples <- selectSampleSubset(samples, ...)
  
  df <- lapply(samples, function(f) {
    #we read in the file 
    screen_data <- read.delim(f, header=T, skip=1)
    #the last line is our %hit_no_genomes
    no_hit <- as.numeric(gsub("%Hit_no_genomes: ", "",screen_data[nrow(screen_data),1]))
    screen_data <- screen_data[-nrow(screen_data),]
    rownames(screen_data) <- screen_data$Genome
    screen_data <- screen_data[, !(colnames(screen_data)=="Genome")]
    #we get the number of unique/multiple hits per one genome and calculate sum (of the percentages)
    one_genome <- data.frame(perc=rowSums(screen_data[, 
                                                      grepl("one_genome.1",
                                                            colnames(screen_data))]))
    one_genome$genome <- rownames(one_genome)
    one_genome$category <- "one genome"
    multi_genome <- data.frame(perc=rowSums(screen_data[, colnames(screen_data) %in% 
                                                          c("X.One_hit_multiple_genomes.1",
                                                            "X.Multiple_hits_multiple_genomes")]))
    multi_genome$genome <- rownames(multi_genome)
    multi_genome$category <- "multiple genomes"
    
    mapping_info <- rbind(one_genome, multi_genome)
    mapping_info <- rbind(mapping_info, c(perc=no_hit, genome="no hit", category="no hit"))
    mapping_info$perc <- as.numeric(mapping_info$perc)
    mapping_info$category <- factor(mapping_info$category, levels=c("one genome","multiple genomes","no hit"))
    return(mapping_info)
  })
  
  # sample names
  samples <- gsub("_screen.txt", "", basename(samples))
  
  if(file.exists(SHINYREPS_TARGET)){
    
    # get target names
    targets <- read.delim(SHINYREPS_TARGET, sep="\t")
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    samples <- sapply(samples, function(i) { ifelse(i %in% targets$sample_ext,
                                                    targets[targets$sample_ext == i,"sample"],
                                                    ifelse(gsub(".R1|.R2","",i) %in% targets$sample_ext,
                                                           ifelse(gsub(".R1","",i) %in% targets$sample_ext,
                                                                  paste0(targets[targets$sample_ext == gsub(".R1","",i),"sample"],".R1"),
                                                                  paste0(targets[targets$sample_ext == gsub(".R2","",i),"sample"],".R2")),
                                                           gsub(paste0("^",SHINYREPS_PREFIX),"",i)))})
  } else {
    if(!is.na(SHINYREPS_PREFIX)) {
      samples <- gsub(paste0("^",SHINYREPS_PREFIX), "", samples)
    }
  }
  
  names(df) <- samples
  #createing a df out of the lists
  df <- reshape2::melt(df, value.name="perc")
  colnames(df) <- gsub("L1", "sample", colnames(df))
  
  # filter for relevant contaminants (e.g. showing >=1% (perc.to.plot) in any sample)
  max.per.genome <- aggregate(df$perc,list(df$genome),max)
  relevant.genomes <- max.per.genome[max.per.genome[,2]>=perc.to.plot,1]
  df <- df[df$genome %in% relevant.genomes,]
  
  # sort alphabetically
  df$sample <- factor(df$sample, levels=unique(df$sample)[order(unique(df$sample),decreasing=TRUE)])
  
  # replace "Mycoplasma" by "Mycoplasma species", FastQScreen itself cannot deal with space characters
  df$genome <- gsub("Mycoplasma","Mycoplasma species",df$genome)
  
  # split/wrap per genome
  df$genome <- factor(df$genome, levels=unique(df$genome))
  
  p.category.wrap <- ggplot(df, aes(x=sample, y=perc, fill=category)) +
    geom_col(position=position_stack(reverse=T),width=0.8) +
    scale_fill_manual(values=c(alpha("#4281a4",0.8),alpha("#ffa62b",0.8),"gray60")) +   
    scale_y_continuous(breaks=seq(0,100,by=10)) +
    theme_bw(base_size=10) +
    labs(x = "",
         y = "% mapped") +
    theme(axis.text.x = element_text(vjust=0.5, angle=90),
          legend.position = "top") +
    guides(fill=guide_legend(title="mapped to", ncol=3)) +
    facet_wrap(~genome,ncol=ncol) +
    coord_flip()
  
  return(list(p.category.wrap=p.category.wrap,
              no.of.genomes=length(unique(df$genome)),
              no.of.samples=length(unique(df$sample)),
              no.of.rows = ceiling(length(unique(df$genome))/ncol))
  )
}



##
## MPShelper.cutadapt: get cutadapt statistics from the log folder and display them
## 
#' @param targetsdf targets data.frame or character with file path to targets object
#' @param colorByFactor character with column name of sample table to be used for coloring the plot. Coloring by filename if NULL. 
#' @param plotfun define function to be used for plotting
#' @param labelOutliers logical, shall outlier samples be labeled
#' @param outlierIQRfactor numeric, factor is multiplied by IQR to determine outlier
#'
#' @return plot cutadapt statistics as side effect
MPShelper.cutadapt <- function(targetsdf=SHINYREPS_TARGET, colorByFactor="group", 
                               plotfun=MPShelper.cutadapt.plot, labelOutliers=T, outlierIQRfactor=1.5, ...)
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
    
    trimmed.R1.perc <- trimmed.R2.perc <- trimmed.reads.perc <- tooshort.reads.perc <- NULL # initialize with NULL in case not needed
    
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
    adapters.perc <- round(100*(as.numeric(adapters) / as.numeric(total.reads)),1)
    adapterprime <- gsub(";", "", system(paste("grep Sequence:", f, "| awk '{print $5}'"), intern=T))
    
    names(adapters.perc) <- gsub(" *=== *", "", system(paste("grep \"=== .*Adapter\"", f), intern=T))
    namespart1 <- gsub("First read:.*", "R1_", names(adapters.perc))
    namespart1 <- gsub("Second read:.*", "R2_", namespart1)
    namespart2 <- gsub("^.*Adapter", "Adapter", names(adapters.perc))
    names(adapters.perc) <- paste0(if(paired) {namespart1} else {""}, adapterprime, namespart2)
    
    ## add trimmed reads for each adapter here
    return(c("total reads"=total.reads, trimmed_R1=trimmed.R1.perc, trimmed_R2=trimmed.R2.perc, 
             trimmed=trimmed.reads.perc, "too short"=tooshort.reads.perc, adapters.perc))
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
      targetsdf <- read.delim(targetsdf)
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
  x.melt <- reshape2::melt(x.df, measure.vars=c(grep("trimmed", colnames(x.df), value=T), 
                                                "too short", 
                                                grep("(Adapter)|(})", colnames(x.df), value=T)), variable.name="reads")
  # everything which is not a value should be a factor
  
  # one plot for each element of colorByFactor
  violin.list <- lapply(colorByFactor, plotfun, data=x.melt, labelOutliers=labelOutliers, outlierIQRfactor=outlierIQRfactor) # "colorByFactor" is submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  DT::datatable(x.df[,c(colorByFactor, "total reads", 
                        grep("trimmed", colnames(x.df), value=T),
                        "too short", 
                        grep("(Adapter)|(})", colnames(x.df), value=T))], 
                options = list(pageLength= 20))
}


# plotting function for MPShelper.cutadapt 
MPShelper.cutadapt.plot <- function(data, color.value, labelOutliers=T, outlierIQRfactor=1.5){
  
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
## MPShelper.pear: get stats for PEAR readpair assembly
## 
MPShelper.pear <- function(colorByFactor=NULL, targetsdf=targets, ...){
  
  x <- list.files(SHINYREPS_PEAR_LOG,pattern='*.log$',full.names=TRUE) 
  x <- x[!grepl("fastq", basename(x))]
  
  # select subset of samples or use all samples if samplePattern=NULL
  x <- selectSampleSubset(x, ...)
  if(length(x) == 0) {
    return("No samples matched with this pattern...")
  }
  
  x <- sapply(x, function(f) { 
    input_reads <- system(paste('grep \"Assembled reads \\.\\.\"', f, "| awk '{print $6}'"), intern=TRUE)
    
    assembled_reads <- system(paste('grep \"Assembled reads \\.\\.\"', f, "| awk '{print $4}'"), intern=TRUE)
    
    discarded_reads <- system(paste("grep \"Discarded reads \\.\\.\"", f, "| awk '{print $4}'"), intern=TRUE)
    
    not_assembled_reads <- system(paste("grep \"Not assembled reads \\.\\.\"", f, "| awk '{print $5}'"), intern=TRUE)
    
    return(c(input_reads=input_reads, assembled_reads=assembled_reads,  
             discarded_reads=discarded_reads, not_assembled_reads=not_assembled_reads))
  })
  
  # set row and column names
  x.df <- as.data.frame(t(x)) 
  x.df <- as.data.frame(lapply(x.df, function(y) {as.numeric(gsub(",", "", y))}))
  
  vars4table <- colnames(x.df)
  vars2plot <- vars4table
  
  #reduce length of file names 
  row.names(x.df) <- basename(colnames(x))
  x.df$filename_unmod <- factor(row.names(x.df))
  if(!is.na(SHINYREPS_PREFIX)) {
    row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))
  }
  row.names(x.df) <- gsub("(?:\\.R1)?\\.cutadapt\\.log$", "", row.names(x.df))
  if(nrow(x.df)>1){
    if(is.na(SHINYREPS_PREFIX)) {row.names(x.df)  <- gsub(lcPrefix(row.names(x.df) ), "", row.names(x.df) )}
    row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
  }
  
  
  # passing the different factors given in targetsdf to x.df which was created from cutadapt logfile names (if 1 cell per file)
  if(!is.null(colorByFactor) && nrow(x.df) == nrow(targetsdf)) { # if targets object fits in length, add information to x.df
    
    index <- sapply(targetsdf$file, grep, x.df$filename_unmod, ignore.case = T) # grep sample name in file names
    targetsdf <- targetsdf[sapply(index, length) ==1, ] # remove targetsdf entries not (uniquely) found in x.df
    
    if(!identical(sort(unname(unlist(index))), 1:nrow(x.df))) {
      return("There seem to be ambiguous sample names in targets. Can't assign them uniquely to logfile names")
    }
    
    x.df <- data.frame(x.df[unlist(index),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    x.df$filename <- x.df$sample
    rownames(x.df) <- x.df$sample
    
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
  x.melt <- melt(x.df, measure.vars=vars2plot,
                 variable="reads")
  #everything which is not a value should be a factor
  
  #now we do a violin plot of the trimmed/too_short/etc. ones and color it
  # according to the different factors given in colorByFactor 
  # prepare palette of appropriate length
  colourCount = length(unique(x.melt[,colorByFactor]))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  create.violin <- function(x.melt, color.value){
    ylab <- "# reads"
    p <- ggplot(x.melt, aes_string(x="reads",
                                   y="value",
                                   color=color.value ))+
      geom_quasirandom() +
      scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
      ylab(ylab) +
      xlab("") +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), strip.text.y = element_text(angle=0))
    return(p)
  }
  
  # one plot each element of colorByFactor
  violin.list <- lapply(colorByFactor, create.violin, x.melt=x.melt) # "colorByFactor" submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  DT::datatable(x.df[,vars4table])
}



##
## MPShelper.umiextract: get stats for barcode extraction from UMI_tools
## 
MPShelper.umiextract <- function(colorByFactor=NULL, targetsdf=targets, ...){
  
  args <- match.call()

  if("samplePattern" %in% names(args) && any(grepl("whitelist", samplePattern, ignore.case = T))) {
    foldername <- SHINYREPS_UMIEXTRACT_LOGWL
  } else {
    foldername <- SHINYREPS_UMIEXTRACT_LOG
  }
  x <- list.files(foldername,pattern='*umibarcode.log$',full.names=TRUE) 
  
  # select subset of samples or use all samples if samplePattern=NULL
  x <- selectSampleSubset(x, ...)
  
  x <- sapply(x, function(f) { 
    input_reads <- system(paste("grep \"INFO Input Reads:\"", f, "| awk '{print $6}'"), intern=TRUE)
    
    match_Regex1 <- system(paste("grep \"INFO regex matches read1:\"", f, "| awk '{print $7}'"), intern=TRUE)
    noMatch_Regex1 <- system(paste("grep \"INFO regex does not match read1:\"", f, "| awk '{print $9}'"), intern=TRUE)
    
    match_Regex2 <- system(paste("grep \"INFO regex matches read2:\"", f, "| awk '{print $7}'"), intern=TRUE)
    noMatch_Regex2 <- system(paste("grep \"INFO regex does not match read2:\"", f, "| awk '{print $9}'"), intern=TRUE)
    
    filtered_BC <- system(paste("grep \"INFO Filtered cell barcode:\"", f, "| awk '{print $7}'"), intern=TRUE)
    
    corrected_BC <- system(paste("grep \"INFO False cell barcode. Error-corrected:\"", f, "| awk '{print $8}'"), intern=TRUE)
    
    not_correctable_BC <- system(paste("grep \"INFO Filtered cell barcode. Not correctable:\"", f, "| awk '{print $9}'"), intern=TRUE)
    
    reads_output <- system(paste("grep \"INFO Reads output:\"", f, "| awk '{print $6}'"), intern=TRUE)
    
    # stats for whitelist extraction:
    BCs_passed_threshold <- system(paste("grep \"cell barcodes passed the selected threshold\"", f, "| awk '{print $5}'"), intern=TRUE)
    
    match_pattern <- system(paste("grep \"reads matched the barcode pattern\"", f, "| awk '{print $4}'"), intern=TRUE)
    
    unique_BCs <- system(paste("grep \"unique cell barcodes\"", f, "| awk '{print $5}'"), intern=TRUE)
    
    Reads_matching_BCs <- system(paste("grep \"total reads matching the selected cell barcodes\"", f, "| awk '{print $5}'"), intern=TRUE)
    
    Reads_corrected_to_BCs <- system(paste("grep \"total reads which can be error corrected to the selected cell barcodes\"", f, "| awk '{print $5}'"), intern=TRUE)
    
    return(c(input_reads=input_reads, match_Regex1=match_Regex1, noMatch_Regex1=noMatch_Regex1,
             match_Regex2=match_Regex2, noMatch_Regex2=noMatch_Regex2,
             filtered_BC=filtered_BC,
             corrected_BC=corrected_BC, not_correctable_BC=not_correctable_BC,
             reads_output=reads_output,
             BCs_passed_threshold=BCs_passed_threshold, match_pattern=match_pattern, unique_BCs=unique_BCs,
             Reads_matching_BCs=Reads_matching_BCs, Reads_corrected_to_BCs=Reads_corrected_to_BCs))
  })
  
  # for samples with very few reads some of the logfile outputs may be missing. 
  # In that case x is returned as list of character vectors of different length.
  if(class(x)[1]=="list") { # length may be >1: "matrix" "array"
    x.df <- lapply(x, function(y) {t(as.data.frame(y))})
    #x.df <- plyr::rbind.fill.matrix(x.df)
    x.df <- dplyr::bind_rows(lapply(x.df, as.data.frame, drop=F))
    x.df <- as.data.frame(x.df)
    x.df <- as.data.frame(lapply(x.df, as.numeric))
    row.names(x.df) <- basename(names(x))
  } else {
    x.df <- as.data.frame(t(x)) 
    x.df <- as.data.frame(lapply(x.df, as.numeric))
    row.names(x.df) <- basename(colnames(x))
  }
  
  
  # define variables to output
  vars4table <- colnames(x.df)
  vars2plot <- vars4table[!vars4table %in% "BCs_passed_threshold"]
  

  #reduce length of file names 
  x.df$filename_unmod <- factor(row.names(x.df))
  if(!is.na(SHINYREPS_PREFIX)) {
    row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))
  }
  row.names(x.df) <- gsub("(?:\\.R1)?\\..*$", "", row.names(x.df))
  if(nrow(x.df)>1){
    if(is.na(SHINYREPS_PREFIX)) {row.names(x.df)  <- gsub(lcPrefix(row.names(x.df) ), "", row.names(x.df) )}
    row.names(x.df)  <- gsub(lcSuffix(row.names(x.df) ), "", row.names(x.df) )
  }
  

  # passing the different factors given in targetsdf to x.df which was created from cutadapt logfile names (if 1 cell per file)
  if(!is.null(targetsdf) && !is.null(colorByFactor) && nrow(x.df) == nrow(targetsdf)) { # if targets object fits in length, add information to x.df
    
    index <- sapply(targetsdf$file, grep, x.df$filename_unmod, ignore.case = T) # grep sample name in file names
    targetsdf <- targetsdf[sapply(index, length) ==1, ] # remove targetsdf entries not (uniquely) found in x.df
    
    if(!identical(sort(unname(unlist(index))), 1:nrow(x.df))) {
      return("There seem to be ambiguous sample names in targets. Can't assign them uniquely to logfile names")
    }
    
    x.df <- data.frame(x.df[unlist(index),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    x.df$filename <- x.df$sample
    rownames(x.df) <- x.df$sample
    
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
  x.melt <- melt(x.df, measure.vars=vars2plot,
                 variable="reads")
  #everything which is not a value should be a factor
  
  #now we do a violin plot of the trimmed/too_short/etc. ones and color it
  # according to the different factors given in colorByFactor 
  # prepare palette of appropriate length
  colourCount = length(unique(x.melt[,colorByFactor]))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  create.violin <- function(x.melt, color.value){
    ylab <- ""
    p <- ggplot(x.melt, aes_string(x="reads",
                                   y="value",
                                   color=color.value ))+
      geom_quasirandom() +
      scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
      ylab(ylab) +
      xlab("") +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), strip.text.y = element_text(angle=0))
    return(p)
  }
  
  # one plot each element of colorByFactor
  violin.list <- lapply(colorByFactor, create.violin, x.melt=x.melt) # "colorByFactor" submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  DT::datatable(x.df[,vars4table])
}



##
## MPShelper.whitelistextraction: plot knee plots from UMI-tools whitelist extraction
##
MPShelper.whitelistextraction <- function(web=F, ...) {
  
  # logs folder
  if(!all(sapply(SHINYREPS_UMIEXTRACT_LOGWL, file.exists))) {
    return(paste("Whitelist extraction statistics not available for", names(which(!sapply(SHINYREPS_UMIEXTRACT_LOGWL, file.exists)))))
  }
  
  SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN[1]),error=function(e){3})
  if(SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 3L    # default to 3 columns
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(SHINYREPS_UMIEXTRACT_LOGWL, pattern="barcode_counts.png$", full.names = T)
  
  # select subset of samples or use all samples if samplePattern=NULL
  samples <- selectSampleSubset(samples, ...)
  
  df <- sapply(samples, function(f) {
    paste0("![whitelistextract img](", f, ")")
  })
  names(df) <- basename(names(df))
  #reduce length of file names 
  if(!is.na(SHINYREPS_PREFIX)) {
    names(df) <- gsub(SHINYREPS_PREFIX, "", names(df))
  }
  if(nrow(df)>1){
    if(is.na(SHINYREPS_PREFIX)) {names(df)  <- gsub(lcPrefix(names(df) ), "", names(df) )}
    names(df)  <- gsub(lcSuffix(names(df) ), "", names(df) )
  }
  
  # put sample names and output in md table of SHINYREPS_PLOTS_COLUMN columns
  while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
  samples <- names(df)
  
  df      <- matrix(df     , ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
  samples <- matrix(samples, ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
  colnames(df.names) <- rep(" ", SHINYREPS_PLOTS_COLUMN)
  
  kable(as.data.frame(df.names), align="c", output=F, format="markdown")
}





##
## MPShelper.Bustard: call the perl XML interpreter and get the MD output
##
MPShelper.Bustard <- function() {
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
## MPShelper.Trackhub: display the UCSC trackhub URL
##
MPShelper.Trackhub <- function() {
  
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
#' isoutlier function using median absolute deviation (MAD) as defined in scran
#' 
MAD <- function(x, n) {
  constant <- 1.4826 # keep consistency with sd in normally distributed data (see ?stats::mad)
  list(low =median(x) - n * constant * median(abs(x - median(x))),
       high=median(x) + n * constant * median(abs(x - median(x))))
}


##
#' plotPCAfromQCmetrics: plot PCA from qc metrics of sce object
#'
#' @param sce SingleCellExperiment object
#' @param metrics character vector with QC metric to be used from colData(sce)
#' @param anno character vector defining up to 2 cell grouping factors to be indicated in the PCA plot by color and shape
#' @param qc.drop dataframe containing logicals for passing the applied qc criteria. May contain a summary column "pass".
#' 
plotPCAfromQCmetrics <- function(sce, metrics, anno, qc.drop=NULL){
  
  pca.sce <- scater::runColDataPCA(sce, variables = metrics, name = "PCA_coldata", ncomponents = 2)  
  
  if(!is.null(qc.drop)) { # use dataframe with QC failing data 
    
    if(! "pass" %in% colnames(qc.drop)) {qc.drop$pass <- !apply(qc.drop, 1, any)} # create "pass" column if not present
    qcfiltercriteria <- colnames(qc.drop)[colnames(qc.drop) != "pass"]
    
    colData(pca.sce) <- cbind(colData(pca.sce), qc.drop[match(colnames(pca.sce), rownames(qc.drop)), ]) # add passing qc metrics to PCs
    
    for (i in qcfiltercriteria) { # create column with combined qc criteria output
      colData(pca.sce)[,paste0(i,".word")] <- ifelse(colData(pca.sce)[,i]==TRUE,i,"")
    }
    pca.sce$QC_drop <- ifelse(pca.sce$pass, "kept", apply(colData(pca.sce)[,paste0(qcfiltercriteria,".word")], 1, paste, collapse="+"))
    pca.sce$QC_drop <- sub("^\\+","",sub("\\+$","",gsub("\\+\\+","\\+",pca.sce$QC_drop)))
    
    anno <- c("QC_drop", anno) # color by failed qc criteria
  }
  
  if(length(anno)>=1) {
    dotcol <- anno[1]
    if(length(anno)>=2) { 
      dotshape <- anno[2]} else {
        dotshape <- NULL}
  } else {
    dotcol <- dotshape <- NULL
  }
  
  p_allGroups <-  plotReducedDim(pca.sce, dimred="PCA_coldata", colour_by=dotcol, shape_by=dotshape, point_size=3) +    
    #scale_fill_discrete(guide=F) +
    scale_color_brewer(palette = "Dark2", name=dotcol) +
    theme_bw()
  
  if("cells" %in% colnames(colData(pca.sce))) { # add labels for 0c and 10c controls
    p_allGroups <- p_allGroups + geom_text_repel(aes(x=PC1, y=PC2, label=cells), 
                                                 subset(data.frame(reducedDim(pca.sce, "PCA_coldata"), cells=pca.sce$cells), cells %in% c("0c", "10c"))) 
  }
  
  plot(p_allGroups)
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
  
  maxno <- as.numeric(maxno) # in case number is character
  if(!gtools::invalid(maxno) && is.numeric(maxno)) {
    samples <- samples[1:min(length(samples), maxno)]
    if(maxno > length(samples)) {cat("\nSample number restricted to", maxno)}
  }
  return(samples)
}
