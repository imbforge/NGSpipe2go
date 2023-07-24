##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("knitr")		# for markdown output
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("ngsReports")

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
VARhelper.Fastqc <- function(web=FALSE, subdir="") {
  
  # logs folder
  if(!all(sapply(file.path(SHINYREPS_FASTQC_OUT, subdir), file.exists))) {
    return(paste("Fastqc statistics not available for", names(which(!sapply(file.path(SHINYREPS_FASTQC_OUT, subdir), file.exists)))))
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.dirs(QC, recursive=F, full.names = T)
  samples <- samples[sapply(samples, function(x) {file.exists(file.path(x, "fastqc_data.txt"))})] # exclude potential subdir which is also listed by list.dirs
  
  df <- sapply(samples, function(f) {
    c(paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
      paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"),
      paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_sequence_gc_content.png)"))
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
## VARhelper.ngsReports.Fastqc: joint FastQC report of all samples in the experiment
##
VARhelper.ngsReports.Fastqc <- function(subdir="") {
  
  # output folder
  if(!file.exists(SHINYREPS_FASTQC_OUT)) {
    return("Fastqc statistics not available")
  }
  
  # Loading FastQC Data 
  f <- list.files(file.path(SHINYREPS_FASTQC_OUT, subdir), pattern="fastqc.zip$", full.names=TRUE)
  x <- ngsReports::FastqcDataList(f)
  lbls <- gsub(paste0("(^", SHINYREPS_PREFIX, "|.fastqc.zip$)"), "", names(x))
  names(lbls) <- gsub(".fastqc.zip", ".fastq.gz", names(x))
  
  print(ngsReports::plotBaseQuals(x, labels=lbls))
  print(ngsReports::plotSeqContent(x, labels=lbls) +
          theme(legend.position="right") +
          guides(fill="none", color="legend") +
          geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                     data=data.frame(base=c("T", "A", "C", "G")),
                     inherit.aes=FALSE, show.legend=TRUE) +
          scale_color_manual("", 
                             values=c("red", "green", "blue", "black"),
                             breaks=c("T", "A", "C", "G"))
  )
  print(ngsReports::plotGcContent(x, plotType="line", gcType="Genome", labels=lbls))  
}

##
## VARhelper.Fastqc.custom: prepare Fastqc summary plots
##
VARhelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, subdir="") {
  
  # logs folder
  if(!file.exists(SHINYREPS_FASTQC_OUT)) {
    return("Fastqc statistics not available")
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) paste0("/fastqc/", subdir) else file.path(SHINYREPS_FASTQC_OUT, subdir)
  
  # read fastqc results in the appropriate format
  f <- list.files(QC, pattern="\\.zip$",full.names=T)
  fastqc.stats <- ngsReports::FastqcDataList(f)
  
  # create proper name vectoir as labels
  lbls <- gsub("_fastqc.zip$", "", basename(names(fastqc.stats)))
  names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", basename(names(fastqc.stats)))
  
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
      guides(color="none",
             fill=guide_legend(title="",ncol=3)) +
      theme(axis.text.x = element_text(size=6,angle=90,hjust=0.5,vjust=0.5),
            axis.text.y = element_text(size=8),
            axis.title  = element_text(size=10),
            plot.title  = element_text(size=12),
            legend.text = element_text(size=7),
            legend.position = "top") +
      scale_x_continuous(breaks=unique(as.numeric(df$position)),
                         labels=unique(df$Base))
    
    p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=rev(lbls)) +
      labs(y = "") +
      theme(legend.position="right") +
      guides(fill="none", color="legend") +
      geom_point(mapping=aes(x=Inf, y=Inf, color=base),
                 data=data.frame(base=c("T", "A", "C", "G")),
                 inherit.aes=FALSE, show.legend=TRUE) +
      scale_color_manual("", 
                         values=c("red", "green", "blue", "black"),
                         breaks=c("T", "A", "C", "G")) 
    
  } else {
    
    p.qual <- ngsReports::plotBaseQuals(fastqc.stats, labels=lbls, plotType="boxplot") +
      theme(axis.text.x = element_text(size=5))
    p.content <- ngsReports::plotSeqContent(fastqc.stats, labels=lbls, plotType="line") +
      theme(axis.text.x = element_text(size=5), legend.position = "top")
  }
  
  # GC content line plot 
  # in case you want to add a theoretical distribution to the plot, use function plotGcContent with 
  # the following settings:
  # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", theoreticalGC=TRUE, species=SPECIES)
  # the default value for SPECIES is "Hsapiens", thus, if you don't specify it, human will be used as a default
  p.gc <- ngsReports::plotGcContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=FALSE) 
  if(!summarizedPlots) {
    p.gc <- p.gc + guides(color=guide_legend(title="",ncol=4)) + 
      theme(legend.position = "top", legend.text = element_text(size=8)) 
  }
  
  return(list(no.of.samples=length(f), p.qual=p.qual, p.content=p.content, p.gc=p.gc))
}


##
## VARhelper.fastqscreen: add FastqScreen data to and plot it as a barplot
##
#' VARhelper.fastqscreen: summarizes FastQScreen results, creates summarized barplots, only relevant contanimants shown
#'
#' @param perc.to.plot - a numeric vector of length 1 setting the percent cutoff of relevant contaminants, if any sample
#'                       shows more than perc.to.plot, contaminant will be shown in plot
#'
#' @return a list including a plot, the number of samples, and the number of plotted contaminants
#'
VARhelper.fastqscreen <- function(perc.to.plot = 1) {
  
  # logs folder
  if(!file.exists(SHINYREPS_FASTQSCREEN_OUT)) {
    return(list(errortext="FastQScreen statistics not available",
                no.of.genomes=1,
                no.of.samples=1))
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- SHINYREPS_FASTQSCREEN_OUT
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(SHINYREPS_FASTQSCREEN_OUT, pattern="_screen.txt$", recursive=T, full.names=T)
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
    targets <- read.delim(SHINYREPS_TARGET)
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    
    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    samples <- sapply(samples, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                    ifelse(sapply("R1", grepl, i), 
                                                           paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R1"),
                                                           ifelse(sapply("R2", grepl, i), 
                                                                  paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R2"),
                                                                  targets[sapply(targets$sample_ext, grepl, i),"sample"])),
                                                    gsub(paste0("^",SHINYREPS_PREFIX),"",i))})                                                    
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
    scale_y_continuous(breaks=seq(0,100,by=10),limits=c(0,100)) +
    theme_bw(base_size=10) +
    labs(x = "",
         y = "% mapped") +
    theme(axis.text.x = element_text(vjust=0.5, angle=90),
          legend.position = "top") +
    guides(fill=guide_legend(title="mapped to", ncol=3)) +
    facet_wrap(~genome,ncol=2) +
    coord_flip()
  
  return(list(p.category.wrap=p.category.wrap,
              no.of.genomes=length(unique(df$genome)),
              no.of.samples=length(unique(df$sample))))      
}


##
## VARhelper.cutadapt: get trimming statistics from the Cutadapt folder and display them
## 
#' @param targetsdf targets data.frame or character with file path to targets object
#' @param colorByFactor character with column name of sample table to be used for coloring the plot. Coloring by filename if NULL. 
#' @param sampleColumnName character with column name(s) of targets table containing file names
#' @param plotfun define function to be used for plotting
#' @param labelOutliers logical, shall outlier samples be labeled
#' @param outlierIQRfactor numeric, factor is multiplied by IQR to determine outlier
#'
#' @return plot cutadapt statistics as side effect
VARhelper.cutadapt <- function(targetsdf=SHINYREPS_TARGET, colorByFactor="group", sampleColumnName =c("file"), 
                               plotfun=VARhelper.cutadapt.plot, labelOutliers=T, outlierIQRfactor=1.5
){
  
  # logs folder
  if(!all(sapply(SHINYREPS_CUTADAPT_STATS, file.exists))) {
    return(paste("Cutadapt statistics not available for", names(which(!sapply(SHINYREPS_CUTADAPT_STATS, file.exists)))))
  }
  
  x <- list.files(SHINYREPS_CUTADAPT_STATS,pattern='*cutadapt.log$',full.names=TRUE) 
  
  # get Command line parameters of first file
  cutadaptpars <- system(paste("grep \"Command line parameters\"", x[1]), intern=T)
  
  paired <- grepl("(-p )|(--paired-output )", cutadaptpars) # output for R2 available?
  
  x <- sapply(x, function(f) { 
    
    trimmed.R1.perc <- trimmed.R2.perc <- trimmed.reads.perc <- trimmed.qual.perc <- tooshort.reads.perc <- toolong.reads.perc <- NULL 
    
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
    #namespart2 <- gsub("^.*Adapter", "Adapter", names(adapters.perc)) # note: this failed when adapter is named and the name contains the word 'Adapter'
    namespart2 <- gsub("First read: Adapter|Second read: Adapter|^Adapter", "Adapter", names(adapters.perc))
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

  if(paired) {

      indexAdapterR1 <- grep("(^-a$)|(--adapter)|(^-g$)|(--front)|(^-b$)|(--anywhere)", cutadaptpars) # index of all R1 adapters
      indexAdapterR2 <- grep("(^-A$)|(^-G$)|(^-B$)", cutadaptpars) # index of all R2 adapters
      indexAdapterSelectedR1 <- indexAdapterR1[grep("[ACGT].[[:digit:]]*}", cutadaptpars[indexAdapterR1+1])] # select e.g. polyA, polyT
      indexAdapterSelectedR2 <- indexAdapterR2[grep("[ACGT].[[:digit:]]*}", cutadaptpars[indexAdapterR2+1])] # select e.g. polyA, polyT

      # rename columns of trimmed polyX 
      if (length(indexAdapterSelectedR1)>0) {
          colnames(x.df)[grepl("Adapter", colnames(x.df)) & grepl("R1", colnames(x.df))][match(indexAdapterSelectedR1, indexAdapterR1)] <- paste0(gsub("Adapter.*$", "", colnames(x.df)[grepl("Adapter", colnames(x.df)) & grepl("R1", colnames(x.df))][match(indexAdapterSelectedR1, indexAdapterR1)]), cutadaptpars[indexAdapterSelectedR1+1])
      }
      if (length(indexAdapterSelectedR2)>0) {
          colnames(x.df)[grepl("Adapter", colnames(x.df)) & grepl("R2", colnames(x.df))][match(indexAdapterSelectedR2, indexAdapterR2)] <- paste0(gsub("Adapter.*$", "", colnames(x.df)[grepl("Adapter", colnames(x.df)) & grepl("R2", colnames(x.df))][match(indexAdapterSelectedR2, indexAdapterR2)]), cutadaptpars[indexAdapterSelectedR2+1])
      }

  } else {

      indexAdapter <- grep("(^-a$)|(--adapter)|(^-g$)|(--front)|(^-b$)|(--anywhere)", cutadaptpars) # index of all adapters applied
      indexAdapterSelected <- indexAdapter[grep("[ACGT].[[:digit:]]*}", cutadaptpars[indexAdapter+1])] # select e.g. polyA, polyT

      # rename columns of trimmed polyX
      if (length(indexAdapterSelected)>0) {
          colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)] <- paste0(gsub("Adapter.*$", "", colnames(x.df)[grepl("Adapter", colnames(x.df))][match(indexAdapterSelected, indexAdapter)]), cutadaptpars[indexAdapterSelected+1])
      }

  }

  # reduce length of file names 
  row.names(x.df) <- basename(colnames(x))
  x.df$filename_unmod <- factor(row.names(x.df))
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
    
    if(length(sampleColumnName)>1) { # melt in case of multiple file name columns (as for ChIP-Seq)
      targetsdf <- targetsdf[,unique(c(colorByFactor, sampleColumnName, "sample"))]
      targetsdf <- reshape2::melt(targetsdf, measure.vars=sampleColumnName, value.name = "filename") 
      for (i in colorByFactor) {targetsdf[, i] <- paste0(targetsdf[, i], " (", targetsdf$variable, ")")}
      targetsdf[,c(colorByFactor, "filename")] <- lapply(targetsdf[,c(colorByFactor, "filename")], factor)
      
    } else {
      targetsdf$filename <- targetsdf[,sampleColumnName]
    }
    
    targetsdf$filename <- gsub("\\..*$", "", targetsdf$filename ) # shorten filename suffix
    index <- sapply(targetsdf$filename, grep, x.df$filename_unmod, ignore.case = T) # grep sample name in file names
    if(is.list(index)) {
      #targetsdf <- targetsdf[sapply(index, length) ==1, ] # remove targetsdf entries not (uniquely) found in x.df
      targetsdf <- targetsdf[sapply(index, length)!=0,] # remove targetsdf entries not found in x.df
      index <- sapply(targetsdf$filename, grep, x.df$filename_unmod, ignore.case = T) # redo grep sample name in file names
    }
    
    if(!identical(sort(unname(unlist(index))), 1:nrow(x.df))) {
      stop("There seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
    }
    
    # x.df <- data.frame(x.df[unlist(index),], targetsdf, check.names =F) ##temp
    x.df <- data.frame(x.df[unlist(t(index)),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    if("sample" %in% colnames(x.df) && !any(duplicated(x.df$sample))) { # use sample column as identifier if present and unique
      x.df$filename <- x.df$sample
      row.names(x.df) <- x.df$sample } 
    
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


# plotting function for VARhelper.cutadapt 
VARhelper.cutadapt.plot <- function(data, color.value, labelOutliers=T, outlierIQRfactor=1.5){
  
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
 
  names(data)[names(data)==color.value] <- "color.value"
 
  p <- ggplot(data, aes(x=reads,
                        y=value,
                        color=color.value))+
    geom_quasirandom(groupOnX=TRUE) +
    geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA)  

  if(labelOutliers) {p <- p + ggrepel::geom_text_repel(data=. %>% filter(!is.na(outlier)), aes(label=filename), show.legend=F)}

  p <- p + scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
    ylab(ylab) +
    xlab("") +
    theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1)) +
    guides(color=guide_legend(title=color.value))
  
  return(p)
}




##
## VARhelper.Trackhub: display the UCSC trackhub URL
##
VARhelper.Trackhub <- function() {
  
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
					 "unmapped"             = paste0( format( x[1,] - x[3,], big.mark=",")," (", round((x[1,] - x[3,]) / x[1,] * 100, digits=2), "%)" ),
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
		
		a <- sapply(c(# "reads were filtered out during the traversal out",  #1 # this probably has to be done separately
		     "MappingQualityReadFilter",                       #1#
				 "MappingQualityAvailableReadFilter",              #2#
				 "MappedReadFilter",                               #3#
				 "NotSecondaryAlignmentReadFilter",                #4#
				 "NotDuplicateReadFilter",                         #5#
				 "PassesVendorQualityCheckReadFilter",             #6#
				 "NonZeroReferenceLengthAlignmentReadFilter",      #7
				 "GoodCigarReadFilter",                            #8#
				 "WellformedReadFilter"),function(y) {             #9#
				 	as.numeric(gsub("(.* - )?(\\d+) read.+","\\2",l[grep(y,l)])) # grep returns line number, then get the respective line ([]) and extract the first number out of it (gsub and replace the whole line with it)
				 })
		
		l.tmp <- l[grep("total reads filtered",l)]
		b <- as.numeric(gsub("(.* - )?(\\d+) total.+","\\2", l.tmp)) #10
		names(b) <- "totalReadsFilt"
		return( c(a, b) )	
		
	})
	
	# transform x from list to matrix (in extreme cases also with only one column)
	x <- do.call(cbind, x)
	# set row and column names, and output the md table
	colnames(x) <- list.files(LOG, pattern=SUFFIX)
	colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX),"",colnames(x))
	colnames(x) <- gsub(paste0(SUFFIX,"$"),"",colnames(x))
	
	df <- data.frame(#"sample name"    = colnames(x),
					 # "total reads"    = format( x[10,], big.mark=","), # no total reads in log anymore
					 "total filtered"   = paste0( format( x[10,], big.mark=",") ),
					 "CIGAR"            = paste0( format( x[8,], big.mark=",") ),
					 "Duplicate"        = paste0( format( x[5,], big.mark=",") ),
					 "MappingQuality"   = paste0( format( x[1,], big.mark=",") ),
					 "Malformed read"   = paste0( format( x[9,], big.mark=",") ),
					 "no MappingQuality"= paste0( format( x[2,], big.mark=",") ),
					 "VendorQuality"    = paste0( format( x[6,], big.mark=",") ),
					 "no RefAlignment"  = paste0( format( x[7,], big.mark=",") ),
					 "not Primary"      = paste0( format( x[4,], big.mark=",") ),
					 "unmapped"         = paste0( format( x[3,], big.mark=",") ) 
					 )
	rownames(df) <- colnames(x)
	kable(df,align=c("r","r","r","r","r","r","r","r","r","r"),output=F)
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
## report version of used tools       
Toolhelper.ToolVersions <- function() {                                                 
    tryCatch({
        ver <- read.delim(file=SHINYREPS_TOOL_VERSIONS)
        colnames(ver) <- c("Tool name","Environment", "Version")
        kable(as.data.frame(ver),output=F)
    }, error=function(e) cat("tool versions not available.\n", fill=TRUE))
}
