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
          scale_color_manual("", values=c("red", "green", "blue", "black"))
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
      scale_color_manual("", values=c("red", "green", "blue", "black")) 
    
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
    scale_y_continuous(breaks=seq(0,100,by=10)) +
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
