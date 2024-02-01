##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("reshape2")
library("ggrepel")
library("grid")
library("knitr")        # for markdown output
library("DESeq2")
library("pheatmap")
library("viridis")
library("gridExtra")
library("dplyr")
library("tidyr")
library("forcats")
library("ngsReports")
library("ggbeeswarm")

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
## smallRNAhelper.subchunkify: small function stolen from here 
## http://michaeljw.com/blog/post/subchunkify/
## to dynamically create chunks and adjust their size accordingly.
##
smallRNAhelper.subchunkify <- function(g, fig_height=7, fig_width=5) {
  g_deparsed <- paste0(deparse(
    function() {grid.draw(g)}
  ), collapse = '')
  
  sub_chunk <- paste0("`","``{r sub_chunk_", floor(runif(1) * 10000),
                      ", fig.height=", fig_height,
                      ", fig.width=", fig_width,
                      ", fig.align='center', echo=FALSE}","\n(", 
                      g_deparsed, ")()",
                      "\n`","``")
  
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}


##		      
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    tryCatch({
        ver <- read.delim(file=SHINYREPS_TOOL_VERSIONS)
        colnames(ver) <- c("Tool name","Environment", "Version")
        kable(as.data.frame(ver),output=F, format="markdown")
    }, error=function(e) cat("tool versions not available.\n", fill=TRUE))
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

##
## smallRNAhelper.Fastqc.custom: prepare Fastqc summary plots
##
smallRNAhelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, subdir="") {
  
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
                                                  gsub(paste0("^",SHINYREPS_PREFIX),"", gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.cutadapt|\\.highQ|\\.trimmed","",lbls))))})
        
        if(SHINYREPS_PAIRED == "yes") {
            x <- names(lbls)
            lbls <- paste0(lbls, ifelse(grepl("R1", names(lbls)), ".R1", ".R2"))
            names(lbls) <- x
        }
    } else {
        
        if(!is.na(SHINYREPS_PREFIX)) {
            lbls <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.cutadapt|\\.highQ|\\.trimmed","",lbls))
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

        ## use rev(lbls) to plot in reverse order
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
    # ngsReports::plotGcContent(fastqc.stats, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=TRUE, species=SPECIES)
    # the default value for SPECIES is "Hsapiens", thus, if you don't specify it, human will be used as a default
    p.gc <- ngsReports::plotGcContent(fastqc.stats, usePlotly=summarizedPlots, plotType="line", gcType="Genome", labels=lbls, theoreticalGC=FALSE) 
    if(!summarizedPlots) {
      p.gc <- p.gc + guides(color=guide_legend(title="",ncol=4)) + 
        theme(legend.position = "top", legend.text = element_text(size=8)) 
    }

    # read length distribution line plot
    p.len <- ngsReports::plotSeqLengthDistn(fastqc.stats, usePlotly=summarizedPlots, plotType="line", labels=lbls)
    if(!summarizedPlots) {
      p.len <- p.len + guides(color=guide_legend(title="",ncol=4)) + 
        theme(legend.position = "top", legend.text = element_text(size=8)) 
    }
        
    return(list(no.of.samples=length(f), 
                p.qual=p.qual, 
                p.content=p.content, 
                p.gc=p.gc,
                p.len=p.len))
}


##
## smallRNAhelper.fastqscreen: add FastqScreen data to and plot it as a barplot
##
#' smallRNAhelper.fastqscreen: summarizes FastQScreen results, creates summarized barplots, only relevant contanimants shown
#'
#' @param perc.to.plot - a numeric vector of length 1 setting the percent cutoff of relevant contaminants, if any sample
#'                       shows more than perc.to.plot, contaminant will be shown in plot
#'
#' @return a list including a plot, the number of samples, and the number of plotted contaminants
#'
smallRNAhelper.fastqscreen <- function(perc.to.plot = 1) {
  
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
                                                    gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.cutadapt|\\.highQ|\\.trimmed","",i)))})                                                    
  } else {
    if(!is.na(SHINYREPS_PREFIX)) {
      samples <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.cutadapt|\\.highQ|\\.trimmed","",samples))
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
## smallRNAhelper.RNAtypes: parse Subread count results for RNAtypes usage
##
smallRNAhelper.RNAtypes <- function() {
    
    FOLDER <- SHINYREPS_RNATYPES
    SUFFIX <- paste0(SHINYREPS_RNATYPES_SUFFIX, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return("Subread statistics not available")
    }

    if(file.exists(SHINYREPS_TARGET)){
        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file)
        add_factors <- colnames(targets)[!colnames(targets) %in% c("group", "sample", "file")]
    }

    # create a matrix using feature names as rownames, sample names as colnames
    l.infiles <- list.files(FOLDER, pattern=SUFFIX, full.names = T)
    l.counts <- lapply(l.infiles, function(f) {
       
        # extract sample name and replace it with nicer sample name given in targets file (if available)
        samplename <- gsub(paste0(SUFFIX,"$"), "", basename(f))
      
        if(file.exists(SHINYREPS_TARGET)){
            # replace files names with nicer sample names given in targets file 
            # if sample is missing in targets file, use reduced file name
            samplename <- ifelse(sum(sapply(targets$sample_ext, grepl, samplename))==1,
                                 targets[sapply(targets$sample_ext, grepl, samplename),"sample"],
                                 gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",samplename)))
        } else {
           if(!is.na(SHINYREPS_PREFIX)) { samplename <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",samplename)) }
        }
    
        df.readcount <- read.table(f, header=T, sep='\t', comment.char = '#', col.names=c("Geneid","Length",samplename))
        df.readcount$Length <- NULL
        return(df.readcount)
    })
    
    # merge counts to data frame in wide format & remove "length..." columns
    f.merge <- function(x,y) {merge(x,y,by="Geneid")}
    df.counts <- Reduce(f.merge, l.counts)
    rm(l.counts)
    
    # eradicate biotype classes that are not present (or could be summarized as "other" being <1%)
    v.countsums <- rowSums(df.counts[-1])
    v.relsums <- v.countsums / sum(v.countsums)
    df.counts$Geneid[v.relsums < 0.01] <- "other"
    
    df.counts <- aggregate(df.counts[-1],
                           by=list(Geneid =df.counts$Geneid),
                           FUN=sum)
    
    # plot
    df.counts.melt <- melt(df.counts, id.var="Geneid")
    colnames(df.counts.melt) <- c("type","sample","count")
    
    # remove possible starting "X" (in case sample name starts with a number)
    df.counts.melt$sample <- gsub("^X","",df.counts.melt$sample)

    # get reverse order for plotting
    samplenames <- gsub("^X","",colnames(df.counts)[-1])
    df.counts.melt$sample <- factor(df.counts.melt$sample, levels=samplenames[order(samplenames, decreasing=TRUE)])

    p.rnatypes <- ggplot() + 
        geom_bar(data=df.counts.melt, aes(x=sample, y=count, fill=type), position="fill", stat="identity", width = 0.8) + 
        labs(x="", y="", fill="") +
	    theme(axis.text.y = element_text(size=8)) +
	    guides(fill = guide_legend(reverse=TRUE)) + 
	    coord_flip() 
    
    return(list(p.rnatypes = p.rnatypes,
                no.of.samples = nrow(df.counts)))
}


##
## smallRNAhelper.subread: parse Subread summary stats and create a md table
##
smallRNAhelper.subread <- function(subdir="",show.ambig=FALSE) {
    
    FOLDER <- paste0(SHINYREPS_SUBREAD,"/",subdir)
    SUFFIX <- paste0(SHINYREPS_SUBREAD_SUFFIX, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return("Subread statistics not available")
    }
    
    # create a matrix using feature names as rownames, sample names as colnames
    x <- sapply(list.files(FOLDER, pattern=SUFFIX), function(f) {
        
        filename <- file(paste0(FOLDER, '/', f))
        l <- readLines(filename)
        close(filename)
        
        # get all interesting rows, extract count and category names
        assigned.and.unassigned <- l[grep("Assigned|Unassigned",l)]
        category.counts <- as.numeric(sapply(assigned.and.unassigned, function(i) {strsplit(i,"\t")[[1]][2]}))
        names(category.counts) <- sapply(assigned.and.unassigned, function(i) {strsplit(i,"\t")[[1]][1]})
        return(category.counts)
    })
    
    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
        # replace files names with nicer sample names given in targets file
	    # if sample is missing in targets file, use reduced file name
        colnames(x) <- sapply(colnames(x), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                                targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                                                gsub(paste0("^",SHINYREPS_PREFIX),"",gsub(paste0("\\",SHINYREPS_SUBREAD_SUFFIX,"|\\.miRNAmature|\\.R1|\\.cutadapt|\\.highQ|\\.trimmed"),"",i)))})
    } else {

        if(!is.na(SHINYREPS_PREFIX)) {
            colnames(x) <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub(paste0("\\",SHINYREPS_SUBREAD_SUFFIX,"|\\.miRNAmature|\\.R1|\\.cutadapt|\\.highQ|\\.trimmed"),"",colnames(x)))
        }
    }

    # remove unmapped to calculate percentages only for mapped ones
    x <- x[rownames(x) != "Unassigned_Unmapped", ]
    x <- rbind(total=x, colSums(x))
    rownames(x)[nrow(x)] <- "total"

    # create data frame for plotting (percentage)
    df.subread.perc <- data.frame(samplename = colnames(x),
                                  assigned   = 100*x["Assigned",]/x["total",],
                                  unassigned_nofeat = 100*x["Unassigned_NoFeatures",]/x["total",],
                                  row.names  = NULL)
    if (show.ambig) {
        df.subread.perc$unassigned_ambig = 100*x["Unassigned_Ambiguity",]/x["total",]
    }

    df.m <- melt(df.subread.perc)
    df.m$samplename <- factor(df.m$samplename, levels=sort(unique(df.m$samplename),decreasing=T))
    if (show.ambig) {
        df.m$variable <- factor(df.m$variable, levels=c("unassigned_nofeat","unassigned_ambig","assigned"))
        my.col <- c("#fb8500","#8ecae6","#219ebc")
        my.breaks <- c("assigned","unassigned_ambig","unassigned_nofeat")
        my.labels <- c("assigned","unassigned ambiguous","unassigned no feature")
    } else {
        df.m$variable <- factor(df.m$variable, levels=c("unassigned_nofeat","assigned"))
        my.col <- c("#fb8500","#219ebc")
        my.breaks <- c("assigned","unassigned_nofeat")
        my.labels <- c("assigned","unassigned no feature")
    }

    # create plot
    p.subread <- ggplot(df.m, aes(samplename,value,fill=variable)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) + 
        labs(x="", y="% of reads") +
        scale_fill_manual(values=my.col,
                          breaks=my.breaks,
                          labels=my.labels) +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = "top") +
        guides(fill=guide_legend(nrow=1,title="",reverse=FALSE)) +
        coord_flip()

    # create md table (omitting various values that are 0 for now)
    df <- data.frame(samplename=colnames(x),
                     total=format(x["total", ], big.mark=","),
                     assigned=paste0(format(x["Assigned", ], big.mark=","), " (", format((x["Assigned", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)"), 
		             unassigned_nofeat=paste0(format(x["Unassigned_NoFeatures", ], big.mark=","), " (", format((x["Unassigned_NoFeatures", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)"))

    if (show.ambig) {
        df$unassigned_ambig=paste0(format(x["Unassigned_Ambiguity", ], big.mark=","), " (", format((x["Unassigned_Ambiguity", ]/x["total", ])*100, digits=2, nsmall=2, trim=T), "%)")
        my.cols <- c("samplename","total","assigned","unassigned_ambig","unassigned_nofeat")
        my.colnames <- c("sample names", "all reads", "assigned (%)", "unassigned ambiguous (%)", "unassigned no feature (%)")
        my.align <- c("l", "r", "r", "r", "r")
    } else {
        my.cols <- c("samplename","total","assigned","unassigned_nofeat")
        my.colnames <- c("sample names", "all reads", "assigned (%)", "unassigned no feature (%)")
        my.align <- c("l", "r", "r", "r")
    }

    # sort alphabetically for plotting
    df <- df[order(df$samplename),]

    return(list(p.subread = p.subread,
                no.of.samples = nrow(df.subread.perc),
                subread.table = kable(df[,my.cols], align=my.align, output=F, format="markdown",row.names=FALSE,
                                      col.names=my.colnames)))
}



##
## smallRNAhelper.subread.type: parse Subread summary stats of one gene_type and create a md table
##
smallRNAhelper.subread.type <- function(subdir="",maindir="all") {
    
    FOLDER    <- paste0(SHINYREPS_SUBREAD,"/",subdir)
    ALLFOLDER <- paste0(SHINYREPS_SUBREAD,"/",maindir)
    SUFFIX    <- paste0("\\.",subdir,SHINYREPS_SUBREAD_SUFFIX_TYPE, '$')
    ALLSUFFIX <- paste0("\\",SHINYREPS_SUBREAD_SUFFIX_TYPE, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return(paste0(subdir," subread statistics not available"))
    }
    if(!file.exists(ALLFOLDER)) {
        return("Subread statistics not available")
    }
    
    # create vector of total counts
    x.type <- sapply(list.files(FOLDER, pattern=SUFFIX), function(f) {
        counts <- sum(read.table(paste0(FOLDER, '/', f))[,2])
        return(counts)
    })
    
    x.all <- sapply(list.files(ALLFOLDER, pattern=ALLSUFFIX), function(f) {
        counts <- sum(read.table(paste0(ALLFOLDER, '/', f))[,2])
        return(counts)
    })
    
    # beautify sample names
    if(file.exists(SHINYREPS_TARGET)){

        # get target names
        targets <- read.delim(SHINYREPS_TARGET)
        targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
        # replace files names with nicer sample names given in targets file
	    # if sample is missing in targets file, use reduced file name
        names(x.type) <- sapply(names(x.type), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                                    targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                                                    gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",gsub(SUFFIX,"",i))))})

        names(x.all) <- sapply(names(x.all), function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1,
                                                                  targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                                                  gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",gsub(ALLSUFFIX,"",i))))})

    } else {

        if(!is.na(SHINYREPS_PREFIX)) {
            names(x.type) <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",gsub(SUFFIX,"",names(x.type))))
            names(x.all)  <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",gsub(ALLSUFFIX,"",names(x.all))))
        }
    }

    # prepare data frames to be merged
    df.type <- data.frame(sample = names(x.type),
                          count.type  = x.type,
                          row.names=NULL)

    df.all  <- data.frame(sample = names(x.all),
                          count.all  = x.all,
                          row.names=NULL)

    # merge + get percentage
    df <- merge(df.type,df.all,by="sample")
    df$count.other <- df$count.all - df$count.type
    
    df.perc <- data.frame(sample = df$sample,
                          perc.type = 100*df$count.type/df$count.all,
                          perc.other = 100*df$count.other/df$count.all)

    df.m <- melt(df.perc)
    df.m$variable <- factor(df.m$variable, levels=c("perc.other","perc.type"))
    df.m$sample   <- factor(df.m$sample, levels=sort(unique(df.m$sample),decreasing=T))

    # create plot
    p.subread.type <- ggplot(df.m, aes(sample,value,fill=variable)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) + 
        labs(x="", y="% of reads") +
        scale_fill_manual(values=c("#fb8500","#219ebc"),
                          breaks=c("perc.type","perc.other"),
                          labels=c(paste0("assigned to ",subdir),"assigned to other RNA types")) +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = "top") +
        guides(fill=guide_legend(nrow=1,title="",reverse=FALSE)) +
        coord_flip()

    # create md table (omitting various values that are 0 for now)
    # DO NOT SHOW THIS TABLE SINCE AMBIGUOUS READS ARE COUNTED DOUBLE AND
    # THUS TOTAL NUMBERS DO NOT MATCH TOTAL SAMPLE ASSIGNED
    df.table <- data.frame(sample=df$sample,
                           total=format(df$count.all, big.mark=","),
                           assigned.type=paste0(format(df$count.type, big.mark=","), " (", format((df$count.type/df$count.all)*100, digits=2, nsmall=2, trim=T), "%)"), 
		                   assigned.other=paste0(format(df$count.other, big.mark=","), " (", format((df$count.other/df$count.all)*100, digits=2, nsmall=2, trim=T), "%)"))

    # sort alphabetically for plotting
    df.table <- df.table[order(df.table$sample),]

    return(list(p.subread.type = p.subread.type,
                no.of.samples = nrow(df.perc),
                subread.type.table = kable(df.table, align=c("l", "r", "r", "r"), output=F, format="markdown",row.names=FALSE,
                                           col.names=c("sample names", "all reads", paste0("assigned to ",subdir,"(%)"), "assigned to other (%)"))))
}



##
## smallRNAhelper.cutadapt.summary: get trimming statistics from the Cutadapt folder and display them
## 
smallRNAhelper.cutadapt.summary <- function() {

    # log files
    LOG <- SHINYREPS_CUTADAPT_STATS
    if(!file.exists(LOG)) {
        return("Trimming logs not available")
    }

    cutadapt <- sapply(list.files(LOG,pattern="*cutadapt.log$"), function(f) {

        # read file
        filename <- file(paste0(LOG, "/", f))
        l <- readLines(filename)
        close(filename)

        # get sample name
        if(file.exists(SHINYREPS_TARGET)){

            # get target names
            targets <- read.delim(SHINYREPS_TARGET)
            targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
            # replace files names with nicer sample names given in targets file
     	    # if sample is missing in targets file, use reduced file name
            samplename <- ifelse(sum(sapply(targets$sample_ext, grepl, f))==1,
                                 targets[sapply(targets$sample_ext, grepl, f),"sample"],
                                 gsub(paste0("^",SHINYREPS_PREFIX),"", gsub("\\.cutadapt\\.log$","",f)))
        } else {

            if(!is.na(SHINYREPS_PREFIX)) {
               samplename <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.cutadapt\\.log$","",f))
            }
        }

        # extract information
        removed.short <- unlist(strsplit(l[grep("Reads that were too short", l)]," "))
        removed.short.no <- as.integer(gsub(",","",removed.short[(length(removed.short)-1)]))
        removed.long <- unlist(strsplit(l[grep("Reads that were too long", l)]," "))
        removed.long.no <- as.integer(gsub(",","",removed.long[(length(removed.long)-1)]))
        kept    <- unlist(strsplit(l[grep("Reads written \\(passing filters\\)", l)]," "))
        kept.no <- as.integer(gsub(",","",kept[(length(kept)-1)]))

        # return
        c(samplename,removed.long.no,removed.short.no,kept.no)
    })

    # extract numbers
    df.cutadapt <- data.frame(samplename=cutadapt[1,],
                              allreads=as.integer(cutadapt[2,])+as.integer(cutadapt[3,])+as.integer(cutadapt[4,]),
                              toolong =as.integer(cutadapt[2,]),
                              tooshort=as.integer(cutadapt[3,]),
                              kept    =as.integer(cutadapt[4,]),
                              row.names=NULL)

    # get fractions of total counts
    df.cutadapt.perc <- data.frame(samplename = df.cutadapt$samplename,
                                   kept    =100*df.cutadapt$kept/df.cutadapt$allreads,
                                   tooshort=100*df.cutadapt$tooshort/df.cutadapt$allreads,
                                   toolong =100*df.cutadapt$toolong/df.cutadapt$allreads)

    df.m <- melt(df.cutadapt.perc)
    df.m$variable <- factor(df.m$variable, levels=c("toolong","tooshort","kept"))
    df.m$samplename <- factor(df.m$samplename, levels=sort(unique(df.m$samplename),decreasing=T))

    # create plot
    p.trim <- ggplot(df.m, aes(samplename,value,fill=variable)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) + 
        labs(x="", y="% of reads") +
        scale_fill_manual(values=c("#fb8500","#8ecae6","#219ebc"),
                          breaks=c("kept","tooshort","toolong"),
                          labels=c("properly trimmed","discarded too short","discarded too long")) +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = "top") +
        guides(fill=guide_legend(nrow=1,title="",reverse=FALSE)) +
        coord_flip()

    # prepare output table
    df <- data.frame(samplename=df.cutadapt$samplename,
                     allreads=format(df.cutadapt$allreads,big.mark=","),
                     kept=paste0(format(df.cutadapt$kep,big.mark=","), " (",format(df.cutadapt.perc$kept,nsmall=2,digits=2),"%)"),
                     tooshort=paste0(format(df.cutadapt$tooshort,big.mark=","), " (",format(df.cutadapt.perc$tooshort,nsmall=2,digits=2),"%)"),
                     toolong=paste0(format(df.cutadapt$toolong,big.mark=","), " (",format(df.cutadapt.perc$toolong,nsmall=2,digits=2),"%)"))

    # sort alphabetically for plotting
    df <- df[order(df$samplename),]

    return(list(p.trim = p.trim,
                df = df.cutadapt,
                no.of.samples = nrow(df.cutadapt.perc),
                trim.table = kable(df, align=c("l", "r", "r", "r", "r"), output=F, format="markdown", row.names=FALSE,
                                   col.names=c("sample names", "all reads", "properly trimmed (%)", "too short (%)","too long (%)"))))

}



##
## smallRNAhelper.qualityfilter: get trimming statistics from the Cutadapt folder and display them
## 
smallRNAhelper.qualityfilter <- function() {

    # log files
    LOG <- SHINYREPS_QUALITYFILTER_LOG
    if(!file.exists(LOG)) {
        return("Quality filtering logs not available")
    }

    qual <- sapply(list.files(LOG,pattern="*\\.log$"), function(f) {

        # read file
        filename <- file(paste0(LOG, "/", f))
        l <- readLines(filename)
        close(filename)

        # get sample name
        if(file.exists(SHINYREPS_TARGET)){

            # get target names
            targets <- read.delim(SHINYREPS_TARGET)
            targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
            # replace files names with nicer sample names given in targets file
     	    # if sample is missing in targets file, use reduced file name
            samplename <- ifelse(sum(sapply(targets$sample_ext, grepl, f))==1,
                                 targets[sapply(targets$sample_ext, grepl, f),"sample"],
                                 gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.R1|\\.R2|\\.cutadapt|\\.fastq_quality_filter|\\.log$","",f)))
        } else {

            if(!is.na(SHINYREPS_PREFIX)) {
               samplename <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.R1|\\.R2|\\.cutadapt|\\.fastq_quality_filter|\\.log$","",f))
            }
        }

        # extract information
        all <- as.integer(unlist(strsplit(l[grep("Input:", l)]," "))[2])
        kept <- as.integer(unlist(strsplit(l[grep("Output:", l)]," "))[2])
        discarded <- as.integer(unlist(strsplit(l[grep("discarded", l)]," "))[2])

        # return
        c(samplename,discarded,kept,all)
    })

    # extract numbers
    df.qual <- data.frame(samplename = qual[1,],
                          discarded  = as.integer(qual[2,]),
                          kept       = as.integer(qual[3,]),
                          all        = as.integer(qual[4,]),
                          row.names=NULL)

    # get fractions of total counts
    df.qual.perc <- data.frame(samplename = df.qual$samplename,
                               kept       = 100*df.qual$kept/df.qual$all,
                               discarded  = 100*df.qual$discarded/df.qual$all)

    df.m <- melt(df.qual.perc)
    df.m$variable <- factor(df.m$variable, levels=c("discarded","kept"))
    df.m$samplename <- factor(df.m$samplename, levels=sort(unique(df.m$samplename),decreasing=T))

    # create plot
    p.qual <- ggplot(df.m, aes(samplename,value,fill=variable)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) + 
        labs(x="", y="% of reads") +
        scale_fill_manual(values=c("#fb8500","#219ebc"),
                          breaks=c("kept","discarded"),
                          labels=c("high quality","discarded low quality")) +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = "top") +
        guides(fill=guide_legend(nrow=1,title="",reverse=FALSE)) +
        coord_flip()

    # prepare output table
    df <- data.frame(samplename=df.qual$samplename,
                     all=format(df.qual$all,big.mark=","),
                     kept=paste0(format(df.qual$kept,big.mark=","), " (",format(df.qual.perc$kept,nsmall=2,digits=2),"%)"),
                     discarded=paste0(format(df.qual$discarded,big.mark=","), " (",format(df.qual.perc$discarded,nsmall=2,digits=2),"%)"))

    # sort alphabetically for plotting
    df <- df[order(df$samplename),]

    return(list(p.qual = p.qual,
                df = df.qual,
                no.of.samples = nrow(df.qual.perc),
                qual.table = kable(df, align=c("l", "r", "r", "r"), output=F, format="markdown", row.names=FALSE,
                                   col.names=c("sample names", "all reads", "high quality (%)", "low quality (%)"))))

}


##
## smallRNAhelper.dedup: plot deduplication stats plot
##
smallRNAhelper.dedup <- function() {

    # log files
    LOG <- SHINYREPS_DEDUP_LOG
    if(!file.exists(LOG)) {
        return("Deduplication logs not available")
    }

    dedup <- sapply(list.files(LOG,pattern="*\\.log$"), function(f) {

        # read file
        filename <- file(paste0(LOG, "/", f))
        l <- readLines(filename)
        close(filename)

        # get sample name
        if(file.exists(SHINYREPS_TARGET)){

            # get target names
            targets <- read.delim(SHINYREPS_TARGET)
            targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
            # replace files names with nicer sample names given in targets file
     	    # if sample is missing in targets file, use reduced file name
            samplename <- ifelse(sum(sapply(targets$sample_ext, grepl, f))==1,
                                 targets[sapply(targets$sample_ext, grepl, f),"sample"],
                                 gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.R1|\\.R2|\\.cutadapt|\\.highQ|\\.dedup|\\.log$","",f)))
        } else {

            if(!is.na(SHINYREPS_PREFIX)) {
               samplename <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.R1|\\.R2|\\.cutadapt|\\.highQ|\\.dedup|\\.log$","",f))
            }
        }

        # extract information
        all <- as.integer(unlist(strsplit(l[grep("inclDup", l)]," "))[2])
        kept <- as.integer(unlist(strsplit(l[grep("exclDup", l)]," "))[2])
        discarded <- all - kept

        # return
        c(samplename,discarded,kept,all)
    })

    # extract numbers
    df.dedup <- data.frame(samplename = dedup[1,],
                           discarded  = as.integer(dedup[2,]),
                           kept       = as.integer(dedup[3,]),
                           all        = as.integer(dedup[4,]),
                           row.names=NULL)

    # get fractions of total counts
    df.dedup.perc <- data.frame(samplename = df.dedup$samplename,
                                kept       = 100*df.dedup$kept/df.dedup$all,
                                discarded  = 100*df.dedup$discarded/df.dedup$all)

    df.m <- melt(df.dedup.perc)
    df.m$variable <- factor(df.m$variable, levels=c("discarded","kept"))
    df.m$samplename <- factor(df.m$samplename, levels=sort(unique(df.m$samplename),decreasing=T))

    # create plot
    p.dedup <- ggplot(df.m, aes(samplename,value,fill=variable)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) + 
        labs(x="", y="% of reads") +
        scale_fill_manual(values=c("#fb8500","#219ebc"),
                          breaks=c("kept","discarded"),
                          labels=c("non-duplicates","removed duplicates")) +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = "top") +
        guides(fill=guide_legend(nrow=1,title="",reverse=FALSE)) +
        coord_flip()

    # prepare output table
    df <- data.frame(samplename=df.dedup$samplename,
                     all=format(df.dedup$all,big.mark=","),
                     kept=paste0(format(df.dedup$kept,big.mark=","), " (",format(df.dedup.perc$kept,nsmall=2,digits=2),"%)"),
                     discarded=paste0(format(df.dedup$discarded,big.mark=","), " (",format(df.dedup.perc$discarded,nsmall=2,digits=2),"%)"))
    
    # sort alphabetically for plotting
    df <- df[order(df$samplename),]

    return(list(p.dedup = p.dedup,
                df = df.dedup,
                no.of.samples = nrow(df.dedup.perc),
                dedup.table = kable(df, align=c("l", "r", "r", "r"), output=F, format="markdown", row.names=FALSE,
                                    col.names=c("sample names", "all reads", "non-duplicates (%)", "removed duplicates (%)"))))

}


##
## smallRNAhelper.rawfiltersummary: plot filtering stats summary plot
##
smallRNAhelper.rawfiltersummary <- function() {

    trim <- smallRNAhelper.cutadapt.summary()$df
    qual <- smallRNAhelper.qualityfilter()$df
    if (SHINYREPS_REMOVE_DUPLICATES) { dedup <- smallRNAhelper.dedup()$df }
 
    # create unique colnames
    names(trim)[!names(trim) %in% "samplename"] <- paste0(names(trim)[!names(trim) %in% "samplename"],".trim")
    names(qual)[!names(qual) %in% "samplename"] <- paste0(names(qual)[!names(qual) %in% "samplename"],".qual")
    if (SHINYREPS_REMOVE_DUPLICATES) { 
       names(dedup)[!names(dedup) %in% "samplename"] <- paste0(names(dedup)[!names(dedup) %in% "samplename"],".dedup") 
    }

    mix.df <- merge(trim[,names(trim)!="kept.trim"],qual[,names(qual)!="all.qual"],by="samplename")
    if (SHINYREPS_REMOVE_DUPLICATES) { 
       mix.df <- merge(mix.df,dedup[,names(dedup)!="all.dedup"],by="samplename")
       mix.df$kept.qual <- NULL
    }

    # get fractions of total counts
    mix.df.perc <- data.frame(samplename = mix.df$samplename,
                              100*mix.df[,!names(mix.df) %in% c("samplename","allreads.trim")]/mix.df$allreads.trim)

    names(mix.df.perc) <- gsub("kept\\.(.*)","kept",names(mix.df.perc))
    names(mix.df)      <- gsub("kept\\.(.*)","kept",names(mix.df))

    # set meaningful names to print
    meaningful.names <- c("too long","too short","low quality","duplicates","remaining/kept")
    names(meaningful.names) <- c("toolong.trim","tooshort.trim","discarded.qual","discarded.dedup","kept")

    df.m <- melt(mix.df.perc)
    df.m$variable <- factor(df.m$variable, levels=names(meaningful.names))
    df.m$samplename <- factor(df.m$samplename, levels=sort(unique(df.m$samplename),decreasing=T))

    # create plot
    p.filter <- ggplot(df.m, aes(samplename,value,fill=variable)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) + 
        labs(x="", y="% of reads") +
        scale_fill_manual(values=c("#fb8500","#ffb703","#8ecae6","#219ebc","#023047"),
                          breaks=rev(names(meaningful.names)),
                          labels=rev(meaningful.names)) +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 11),
              legend.position = "top") +
        guides(fill=guide_legend(nrow=1,title="",reverse=FALSE)) +
        coord_flip()

    # prepare output table
    df <- data.frame(samplename=mix.df$samplename,
                     all=format(mix.df$allreads.trim,big.mark=","),
                     kept=paste0(format(mix.df$kept,big.mark=","), " (",format(mix.df.perc$kept,nsmall=2,digits=2),"%)"),
                     toolong=paste0(format(mix.df$toolong.trim,big.mark=","), " (",format(mix.df.perc$toolong.trim,nsmall=2,digits=2),"%)"),
                     tooshort=paste0(format(mix.df$tooshort.trim,big.mark=","), " (",format(mix.df.perc$tooshort.trim,nsmall=2,digits=2),"%)"),
                     qual=paste0(format(mix.df$discarded.qual,big.mark=","), " (",format(mix.df.perc$discarded.qual,nsmall=2,digits=2),"%)"))

    if (SHINYREPS_REMOVE_DUPLICATES) { 
        df$dedup <- paste0(format(mix.df$discarded.dedup,big.mark=","), " (",format(mix.df.perc$discarded.dedup,nsmall=2,digits=2),"%)")
    }

    # sort alphabetically for plotting
    df <- df[order(df$samplename),]

    return(list(p.filter = p.filter,
                no.of.samples = nrow(mix.df.perc),
                filter.table = kable(df, 
                                     align=if(SHINYREPS_REMOVE_DUPLICATES) {c("l","r","r","r","r","r","r")} else {c("l","r","r","r","r","r")},
                                     output=F, format="markdown", row.names=FALSE,
                                     col.names=if(SHINYREPS_REMOVE_DUPLICATES) {
                                                  c("sample names", "all reads", "kept/remaining (%)", "too long (%)", "too short (%)", "low quality (%)", "duplicates (%)")
                                               } else {
                                                  c("sample names", "all reads", "kept/remaining (%)", "too long (%)", "too short (%)", "low quality (%)")
                                               })))

}



##
## smallRNAhelper.bowtie: parse bowtie log files and create plot and table
##
smallRNAhelper.bowtie <- function() {
    
    # log file
    LOG <- SHINYREPS_BOWTIE_RES
    SUFFIX <- paste0(SHINYREPS_BOWTIE_SUFFIX, '$')
    if(!file.exists(LOG)) {
        return("Bowtie statistics not available")
    }
    
    # look for the lines containing the strings
    # and get the values associated with this strings
    x <- sapply(list.files(LOG, pattern=SUFFIX), function(f) {
        filename <- file(paste0(LOG, "/", f))
        l <- readLines(filename)
        close(filename)

        # get sample name
        if(file.exists(SHINYREPS_TARGET)){

            # get target names
            targets <- read.delim(SHINYREPS_TARGET)
            targets$sample_ext <- gsub("\\..*$", "", targets$file )
    
            # replace files names with nicer sample names given in targets file
     	    # if sample is missing in targets file, use reduced file name
            samplename <- ifelse(sum(sapply(targets$sample_ext, grepl, f))==1,
                                 targets[sapply(targets$sample_ext, grepl, f),"sample"],
                                 gsub(paste0("^",SHINYREPS_PREFIX),"",gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",gsub(SUFFIX,"",f))))
        } else {

            if(!is.na(SHINYREPS_PREFIX)) {
               samplename <- gsub(paste0("^",SHINYREPS_PREFIX), "", gsub("\\.R1|\\.cutadapt|\\.highQ|\\.trimmed","",gsub(SUFFIX,"",f)))
            }
        }

        # extract information
        all <- gsub(".+: (.+)", "\\1", l[grep("reads processed", l)])
        unique <- gsub(".+: (.+) (\\(.*\\))", "\\1", l[grep("reads with at least one", l)])
        multi  <- gsub(".+: (.+) (\\(.*\\))", "\\1", l[grep("reads with alignments sampled due to -M", l)])
        unmapped <- gsub(".+: (.+) (\\(.*\\))", "\\1", l[grep("reads that failed to align", l)])

        # return
        c(samplename,unmapped,multi,unique,all)

    })
    
    # extract numbers
    df.bowtie <- data.frame(samplename = x[1,],
                            unmapped  = as.integer(x[2,]),
                            multi     = as.integer(x[3,]),
                            unique    = as.integer(x[4,]),
                            all       = as.integer(x[5,]),
                            row.names=NULL)

    # get fractions of total counts
    df.bowtie.perc <- data.frame(samplename = df.bowtie$samplename,
                                 unique     = 100*df.bowtie$unique/df.bowtie$all,
                                 multi      = 100*df.bowtie$multi/df.bowtie$all,
                                 unmapped   = 100*df.bowtie$unmapped/df.bowtie$all)

    # melt and combine for plotting
    df.m <- melt(df.bowtie)
    df.m$value_info <- "reads"

    df.m.perc <- melt(df.bowtie.perc)
    df.m.perc$value_info <- "perc"

    df.m.all <- rbind(df.m,df.m.perc)
    df.m.all$variable <- factor(df.m.all$variable, levels=c("unmapped","multi","unique","all"))
    df.m.all$samplename <- factor(df.m.all$samplename, levels=sort(unique(df.m.all$samplename),decreasing=T))

    ## define color scheme
    my.palette <- rev(c("#fb8500","#8ecae6","#219ebc"))
    my.color   <- "#cc0e25"

    # max number of reads in any sample
    max.reads <- max(df.m.all[df.m.all$variable=="all","value"])

    # plot a combination of percent mapped reads and total sequenced reads
    p.bowtie.perc.count <- ggplot() + 
      geom_bar(data=df.m.all[df.m.all$value_info == "reads" & df.m.all$variable!="all",], 
               mapping=aes(x = samplename, y = value, fill = variable), 
               stat = "identity", position = "fill", width = 0.8) + 
      geom_point(data=df.m.all[df.m.all$variable=="all",], 
                 mapping=aes(x = samplename, y = value/max.reads), 
                 size=2, fill=my.color, color=my.color, shape=18) +
      geom_line(data=df.m.all[df.m.all$variable=="all",],
                mapping=aes(x = samplename, y = value/max.reads, group=1),
                color=my.color, linetype="dashed", size=0.2) +
      scale_y_continuous(sec.axis = sec_axis(~ . *max.reads, name="# sequenced reads"),
                         labels = scales::percent_format()) + 
      labs(x = "",
           y = "% of reads") + 
      theme(axis.title=element_text(size = 11),
            axis.title.x.top=element_text(color=my.color), 
            axis.text=element_text(size = 10),
            axis.text.x.top=element_text(color=my.color), 
            axis.ticks.x.top=element_line(color=my.color),
            legend.text = element_text(size = 11),
            legend.position = "top") + 
      guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) + 
      scale_fill_manual(values=my.palette,
			breaks=c("unmapped","multi","unique"),
			labels=c("unmapped","multiple loci","unique")) + 
      coord_flip()

    # prepare output table
    df <- data.frame(samplename=df.bowtie$samplename,
                     all=format(df.bowtie$all,big.mark=","),
                     unique=paste0(format(df.bowtie$unique,big.mark=","), " (",format(df.bowtie.perc$unique,nsmall=2,digits=2),"%)"),
                     multi=paste0(format(df.bowtie$multi,big.mark=","), " (",format(df.bowtie.perc$multi,nsmall=2,digits=2),"%)"),
                     unmapped=paste0(format(df.bowtie$unmapped,big.mark=","), " (",format(df.bowtie.perc$unmapped,nsmall=2,digits=2),"%)"))

    # sort alphabetically for plotting
    df <- df[order(df$samplename),]

    return(list(p.bowtie = p.bowtie.perc.count,
                no.of.samples = nrow(df.bowtie.perc),
                bowtie.table = kable(df, align=c("l", "r", "r", "r", "r"), output=F, format="markdown", row.names=FALSE,
                                     col.names=c("sample names", "all reads", "unique (%)", "multi (%)", "unmapped (%)"))))

}



##
## smallRNAhelper.bowtie.params: get Bowtie parameters
##
smallRNAhelper.bowtie.params <- function() {

    # log file
    LOG <- SHINYREPS_BOWTIE_RES
    SUFFIX <- paste0(SHINYREPS_BOWTIE_SUFFIX, '$')
    if(!file.exists(LOG)) {
        return("Bowtie parameters not available")
    }
    
    # parameters are the same for all files, thus, use the first file for extraction
    parameter.file <- paste0(LOG,"/",list.files(LOG, pattern=SUFFIX)[1])

    # grep the information
    bowtie.flags <- system(paste("grep \"BOWTIE_FLAGS\"",parameter.file,"| awk 'BEGIN{ORS=\" \"}{for (i=2; i<=NF; i++) print $i;}'"), intern=TRUE)
    bowtie.ref   <- system(paste("grep \"BOWTIE_REF\"",parameter.file,"| awk '{print $2;}'"), intern=TRUE)

    # combine into a data.frame
    bowtie.para.df <- data.frame(setting=c("bowtie flags","bowtie reference"),set_values=c(bowtie.flags,bowtie.ref))

    # print
    kable(bowtie.para.df, align=c("l", "r"), output=F, format="markdown")
}



##
## DESeq2 DE analysis
##

#' Create MDS plot from DESeq2 object. 
#'
#' @param dds - a DESeq2 analysis S4 object
#' @param rld - a rlog transformed DESeq2 object
#'
#' @return MDS plot
#'
#' @description This function will use a DESeq2 differential analysis object to create an MDS plot,
#'              which is labelled non-overlapping gene names.
#' 
#' @return A printed plot (ggplot2) object.
#' 
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           rld <- rlog(dds)
#'           # rld and dds are taken from environment
#'           smallRNAhelper.DESeq2.MDS()
#' 
smallRNAhelper.DESeq2.MDS <- function() {

	pca.data <- plotPCA(rld, intgroup=colnames(colData(dds))[1], returnData=TRUE)
	percentVar <- round(100 * attr(pca.data, "percentVar"))
	ggplot(pca.data, aes(PC1, PC2, color=group)) +
 		geom_point(size=2,alpha=0.7) +
 		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		coord_fixed() + 
  		scale_color_manual(values=define.group.palette(length(levels(colData(dds)[,"group"])))) +
  		geom_text_repel(aes(label=colData(dds)[,"replicate"]), show.legend=FALSE) + 
  		theme_bw()

}


#' Create (pairwise) MDS plot from DESeq2 object.
#'
#' @param pairwise.dds - a list of DESeq2 analysis S4 objects
#' @param i - index of the list element to be plotted
#'
#' @return MDS plot
#'
#' @description This function will use a DESeq2 differential analysis object to create an MDS plot,
#'              which is labelled non-overlapping gene names.
#'
#' @return A printed plot (ggplot2) object.
#'
#' @examples Taken from DE_DESeq2.R
#'           pairwise.dds <- lapply(conts[,1],function(cont) {
#'                       ...
#'           })
#'           # pairwise.dds are taken from environment
#'           smallRNAhelper.DESeq2.pairwisePCA(i)
#'
smallRNAhelper.DESeq2.pairwisePCA <- function(i=1) {

	pca.data <- plotPCA(rlog(pairwise.dds[[i]]), intgroup=colnames(colData(pairwise.dds[[i]]))[1], returnData=TRUE)
     	percentVar <- round(100 * attr(pca.data, "percentVar"))
     	print(ggplot(pca.data, aes(PC1, PC2, color=group)) +
		   geom_point(size=2,alpha=0.7) +
		   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
		   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		   coord_fixed() + 
	 	   scale_color_manual(values=brewer.pal(9,"Set1")[1:2]) + 
		   geom_text_repel(aes(label=colData(pairwise.dds[[i]])[,"replicate"]), show.legend=FALSE) +
	 	   theme_bw())
}


## smallRNAhelper.DESeq2.heatmap 
#' Heatmap of sample to sample distances, of top variant 'n' genes of the counts-per-million table or 
#' of top 'n' highest mean genes of the regularized log transformed counts-per-million table 
#' with or without rlog normalization.
#'
#' @param i numeric index, only needed if dds is a list of objects to select a list element.
#' @param dds - rlog transformed DESeq2 object.
#' @param logTransform logical, if TRUE, dds will be transformed using DESeq2::rlog().
#' @param anno_factors character vector with factors given in dds to be used for heatmap column annotation.
#' @param type character with type of heatmap to be drawn. One of "distance", "cluster_sd" or "cluster_mean".
#' @param n amount of transcripts/rows to be plotted (default = 40). Applies only for type "cluster_sd" and "cluster_mean".
#' 
#' @return A plot (base plotting) object.
#'
#' @examples Taken from DE_DeSeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           rld <- rlog(dds)
#'           # rld is taken from environment
#'           smallRNAhelper.DESeq2.heatmap()
#'           
smallRNAhelper.DESeq2.heatmap <- function(i=NULL, dds=rld, logTransform=F, anno_factors = c("group","replicate"), type="distance", n=40, mature=FALSE) {
  
  # for pairwise object select ith element
  if(!is.null(i)) {dds <- dds[[i]]}
  
  # extract assay, optionally logtransform (as for the pairwise.dds objects) and over-write gene_id in rownames with gene_name
  if(logTransform) {
    if(class(dds) != "DESeqTransform") {
      dds <- DESeq2::rlog(dds)
    } else {
      warning("The dds object is of class DESeqTransform, which means it is already transformed. Therefore, the log transformation selected here is omitted.")
    }
  } 
  assay.rld <- assay(dds)

  if (mature==FALSE) {
      gene_ids <- gtf$gene_name[match(rownames(assay(dds)), gtf$gene_id)]
      rownames(assay.rld) <- ifelse(is.na(gene_ids), rownames(assay(dds)), gene_ids)
  }
  
  ## set parameter specific for heatmap type
  switch(type,
         distance={
           # extract sample to sample distances
           sampleDists <- dist(t(assay.rld))
           sampleDistMatrix <- as.matrix(sampleDists)
           # sample to sample distance heatmap w/ 0's on diagonal set to NA
           # since they otherwise shift the scale too much
           sampleDistMatrix.diagNA <- sampleDistMatrix
           sampleDistMatrix.diagNA[sampleDistMatrix.diagNA==0] <- NA
           data2plot <- sampleDistMatrix.diagNA
           # set color scheme
           col <- plasma(255)
           hmmain = NA
           hm_fontsize_row = 8
           hm_fontsize_col = 8
           clustering_distance_rows = sampleDists
           clustering_distance_cols = sampleDists
           cluster_rows=TRUE
         },
         cluster_sd={
           # pick top n most variable genes
           select.highestSD <- order(apply(assay.rld,1,sd),decreasing=TRUE)[1:n]
           data2plot <- assay.rld[select.highestSD,]
           # set color scheme
           col <- colorRampPalette(brewer.pal(9,"GnBu"))(255)
           hmmain=paste("Normalized expression values of",n,"most variable genes") # heatmap title
           hm_fontsize_row = 6
           hm_fontsize_col = 8
           clustering_distance_rows = "euclidean"
           clustering_distance_cols = "euclidean"
           cluster_rows=TRUE
         },
         cluster_mean={
           # pick top n most variable genes
           select.highestMean <- order(rowMeans(assay.rld),decreasing=TRUE)[1:n]
           data2plot <- assay.rld[select.highestMean,]
           # set color scheme
           col <- rev(heat.colors(255))
           hmmain=paste("Normalized expression values of",n,"genes with highest mean") # heatmap title
           hm_fontsize_row = 6
           hm_fontsize_col = 8
           clustering_distance_rows = "euclidean"
           clustering_distance_cols = "euclidean"
           cluster_rows=FALSE
         }
  )
  
  ## extract information for legend if available
  anno_factors <- anno_factors[anno_factors %in% colnames(colData(dds))]
  if(length(anno_factors)>0) {
    
    legend.df <- data.frame(colData(dds)[,c(anno_factors),drop=F],row.names=rownames(colData(dds)))
    
    # fix group colors for legend additional factors if available
    legend_colors <- list()
    
    for (f in anno_factors) {
      if (f == "group") {
        mypalette <- define.group.palette(length(levels(colData(dds)[,f])))
        legend_colors[[f]] <- setNames(mypalette, levels(colData(dds)[,f]))
      } else {
        mypalette <- define.replicate.palette(length(unique(colData(dds)[,f])))
        legend_colors[[f]] <- setNames(mypalette, sort(unique(colData(dds)[,f])))
      }
    }
  } else {
    legend.df <- NA
    legend_colors <- NA
  }  
  
  
  # plot heatmap
  pheatmap(data2plot,
           cluster_rows=cluster_rows,
           cluster_cols=TRUE,
           clustering_distance_rows=clustering_distance_rows,
           clustering_distance_cols=clustering_distance_cols,
           annotation_col=legend.df,
           annotation_colors=legend_colors,
           col=col,
           border_color=NA,
           main=hmmain, 
           fontsize_row=hm_fontsize_row,
           fontsize_col=hm_fontsize_col,
           treeheight_row=20,
           treeheight_col=20,
           annotation_names_col=T)
}


## smallRNAhelper.MAplot: MA plots
#' MA plots of DESeq2 results
#'
#' @param i - iterator to select contrast from conts [default = 1]
#' @param fdr - FDR cutoff [default = 0.01]
#' @param conts - list of contrasts
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#' @return MA plot with genes padj < fdr hichlighted
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           conts <- c("shMed12vsshNMC=(shMed12-shNMC)")
#'           # res & conts are taken from environment
#'           smallRNAhelper.DESeq2.MAplot()
#'           
smallRNAhelper.DESeq2.MAplot <- function(i=1, fdr=.01) {
	cont.name <- gsub("(.+)=(.+)","\\1",conts[i,1])
	plotMA(res[[i]], main=cont.name, alpha=fdr)
}

## smallRNAhelper.DEgenes: show the DE results
#' Show the DE results ordered by adjusted pValue.
#'
#' @param i - iterator to select precalculated DESeq2 result of contrast from result list
#'            [default = 1]
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#' @return data frame ordered by adjusted pvalue
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           # res is taken from environment
#'           smallRNAhelper.DESeq2.DEgenes(i=2)
#'           
smallRNAhelper.DESeq2.DEgenes <- function(i=1) {
    ord  <- order(-log(res[[i]]$padj), 
                   abs(res[[i]]$log2FoldChange), 
                  decreasing=TRUE)
    res[[i]][ord, ]
}


## smallRNAhelper.DESeq2.VolcanoPlot: Volcano plots from DEseq2 results
#' Produce volcano plots from DEseq2 results.
#'
#' @param i - integer, iterator to determine which contrast's DESeq2 result to plot
#' @param fdr - float, FDR cut off of genes to be highlighted
#' @param top - integer, count of genes to be labelled in Volcano plot
#' @param web - boolean
#' @param res - DE Seq2 result data frame containing log2FC, pvalue, padj, gene_name, ...;
#'              rownames(res) <- ENSEMBL gene IDs
#'
#'
#' @return ggplotly object, if web = TRUE
#'         printed ggplot2 object, if web = FALSE
#'
#' @examples Taken from DE_DeSeq2.R
#'           res <- lapply(conts, function(cont){
#'              ...
#'           })
#'           
#'           # res is taken from environment
#'           smallRNAhelper.DESeq2.VolcanoPlot(i = 1, fdr = 0.01, top = 20, web = FALSE)
#'           
smallRNAhelper.DESeq2.VolcanoPlot <- function(i=1, fdr=.01, top=25, web=TRUE) {
    # gather results
    d <- as.data.frame(smallRNAhelper.DESeq2.DEgenes(i))
    x.limit <- max(abs(d$log2FoldChange), na.rm=T) # find the maximum spread of the x-axis
    
    # plotting
    p <- ggplot(d) +
            geom_point(mapping=aes(log2FoldChange, -log10(padj), size=log10(baseMean+1), color=padj < fdr), alpha=0.1) +
            theme_bw() +
            xlim(-x.limit, x.limit) +
            ylab("-log10 adj. p-value") +
            xlab("log2 fold change") + 
            scale_color_manual(values=c("black", "red")) +
            scale_size_continuous("mean norm. counts (log10)") +
	    guides(color="none", 
		   size = guide_legend(nrow=1)) +
       	    theme(legend.position = "top")

    # add name of top genes
    if(top > 0) p <- p + geom_text_repel(data=d[1:min(top, nrow(d)),],
                                         mapping=aes(log2FoldChange, -log10(padj), label=gene_name),
                                         color="black")

    # return plot
    print(p)
    #return(list(p = p))
}


#' Barplot of sign. vs expressed genes 
#'
#' @param i - iterator to select contrast from conts [default = 1]
#' @param rld - rlog transformed DESeq2 object
#' @return three barplots.
#'
#' @examples Taken from DE_DESeq2.R
#'           dds <- DESeqDataSetFromHTSeqCount(sampleTable=targets,
#'                                             directory=cwd, 
#'                                             design=as.formula(mmatrix))
#'           dds <- DESeq(dds)
#'           rld <- rlog(dds)
#'           # rld is taken from environment
#'           smallRNAhelper.fractionSign.pairwise()
#'           
smallRNAhelper.fractionSign.pairwise <- function(i=1,fdr=.01) {

    res.df <- as.data.frame(res[[i]])

    # extract expressed/tested entries
    res.tested.df <- res.df[!is.na(res.df$padj),]

    # extract type of RNA
    res.tested.df$gene_type <- data.frame(gtf[match(rownames(res.tested.df), gtf$gene_id),])[,SHINYREPS_GENETYPE]

    # get counts per RNA type
    rnatypes.total.df <- as.data.frame(table(res.tested.df$gene_type))
    rnatypes.sign.df  <- as.data.frame(table(res.tested.df[res.tested.df$padj<fdr,"gene_type"]))
    names(rnatypes.total.df) <- c("rnatype","count.total")

    if (nrow(rnatypes.sign.df) > 0) {

        names(rnatypes.sign.df)  <- c("rnatype","count.sign")
        
        # merge and replace NA by 0
        rnatypes.all.df <- merge(rnatypes.total.df,rnatypes.sign.df,by="rnatype",all=TRUE)
        rnatypes.all.df[is.na(rnatypes.all.df)] <- 0

        # get percentages
        rnatypes.all.perc.df <- data.frame(rnatype=rnatypes.all.df$rnatype,        
                                           perc.total=100*(rnatypes.all.df$count.total/sum(rnatypes.all.df$count.total)),
                                           perc.sign=100*(rnatypes.all.df$count.sign/sum(rnatypes.all.df$count.sign)))

        # get percentages relative to it's own category
        rnatypes.all.perc.own.df <- data.frame(rnatype=rnatypes.all.df$rnatype,
                                               perc.own=100*(rnatypes.all.df$count.sign/rnatypes.all.df$count.total))

        # sort levels
        rnatypes.all.df$rnatype      <- factor(rnatypes.all.df$rnatype, 
                                               levels=rnatypes.all.df$rnatype[order(rnatypes.all.df$count.total,decreasing=FALSE)])
        rnatypes.all.perc.df$rnatype <- factor(rnatypes.all.perc.df$rnatype,
                                               levels=rnatypes.all.perc.df$rnatype[order(rnatypes.all.df$count.total,decreasing=TRUE)]) 
        rnatypes.all.perc.own.df$rnatype <- factor(rnatypes.all.perc.own.df$rnatype,
                                                   levels=rnatypes.all.df$rnatype[order(rnatypes.all.df$count.total,decreasing=FALSE)])

        # melt
        rnatypes.all.df.m      <- melt(rnatypes.all.df)
        rnatypes.all.perc.df.m <- melt(rnatypes.all.perc.df)

        rnatypes.all.df.m$variable <- factor(rnatypes.all.df.m$variable, levels=rev(c("count.total","count.sign")))

        # get how many colors are needed
        colorCount <- nrow(rnatypes.all.perc.df)

        # basic color palette
        getPalette = colorRampPalette(brewer.pal(12, "Paired"))

        # plot w/ extended color palette (because rainbow palette is hard to read)
        p.perc100 <- ggplot(rnatypes.all.perc.df.m,aes(variable,value)) + 
            geom_col(aes(fill=rnatype),position="stack") + 
            scale_fill_manual(values = getPalette(colorCount)) +
            scale_x_discrete(breaks=c("perc.total", "perc.sign"),
                             labels=c("detected", "sign. changed")) +
            theme(axis.text.x=element_text(colour="grey30",size=12,angle=30,vjust=1,hjust=1),
                  axis.text.y=element_text(colour="grey30",size=12),
                  axis.title=element_text(colour="grey30",size=14),
                  title=element_text(size=14,hjust=0.5),
                  legend.text=element_text(size=8,colour="grey30")) +
            labs(x="",
                 y="%",
                 title="") + # title="detected and significantly changed genes (percentage)") +
            guides(fill=guide_legend(title="",ncol=2))

        p.count <- ggplot(rnatypes.all.df.m,aes(rnatype,value)) + 
            geom_col(aes(fill=variable),position="dodge") + 
            scale_fill_manual(values = c("#219ebc","#fb8500"),
                              breaks=c("count.total", "count.sign"),
                              labels=c("detected", "significantly changed")) +
            theme(axis.text.x=element_text(colour="grey30",size=12),
                  axis.text.y=element_text(colour="grey30",size=12),
                  axis.title=element_text(colour="grey30",size=14),
                  title=element_text(size=14,hjust=0.5),
            legend.text=element_text(size=12,colour="grey30")) +
            labs(x="",
                 y="gene count") + # title="detected and significantly changed genes (absolute count)") +
            theme(legend.position="top")+
            guides(fill=guide_legend(title="")) +
            coord_flip()

        p.perc <- ggplot(rnatypes.all.perc.df.m,aes(rnatype,value)) + 
            geom_col(aes(fill=variable),position="dodge") + 
            scale_fill_manual(values = c("#219ebc","#fb8500"),
                              breaks=c("perc.total", "perc.sign"),
                              labels=c("detected", "significantly changed")) +
            theme(axis.text.x=element_text(colour="grey30",size=12,angle=45,vjust=1,hjust=1),
                  axis.text.y=element_text(colour="grey30",size=12),
                  axis.title=element_text(colour="grey30",size=14),
                  title=element_text(size=14,hjust=0.5),
                  legend.text=element_text(size=12,colour="grey30")) +
            labs(x="",
                 y="% of all") + # title="detected and significantly changed genes (percentage)") +
            theme(legend.position="top")+
            guides(fill=guide_legend(title=""))

        p.perc.own <- ggplot(rnatypes.all.perc.own.df,aes(rnatype,perc.own)) +
            geom_col(fill="grey50",position="dodge",width=0.8) +
            theme(axis.text.x=element_text(colour="grey30",size=12),
                  axis.text.y=element_text(colour="grey30",size=12),
                  axis.title=element_text(colour="grey30",size=14),
                  title=element_text(size=14,hjust=0.5),
                  legend.text=element_text(size=12,colour="grey30")) +
            labs(x="",
                 y="% significant",
                 title="") + # title="Percentage of sign. changed genes in each type") +
            theme(legend.position="top")+
            guides(fill=guide_legend(title="")) +
            coord_flip()

        return(list(p.perc100 = p.perc100,
                    p.count   = p.count,
                    p.perc    = p.perc,
                    p.perc.own= p.perc.own,
                    rnatypes.total = nrow(rnatypes.all.df),
                    rnatypes.sign  = nrow(rnatypes.sign.df)))
    } else {
        #cat("NO SIGNIFICANT GENES FOUND.\n")
        return(rnatypes.sign = nrow(rnatypes.sign.df))
    }

}

























