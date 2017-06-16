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
## smallRNAhelper.bowtie.table: parse bowtie log files and create a md table
##
smallRNAhelper.bowtie.table <- function() {
    
    # log file
    LOG <- SHINYREPS_BOWTIE_RES
    SUFFIX <- paste0(SHINYREPS_BOWTIE_SUFFIX, '$')
    if(!file.exists(LOG)) {
        return("Bowtie statistics not available")
    }
    
    # look for the lines containing the strings
    # and get the values associated with this strings
    x <- sapply(list.files(LOG, pattern=SUFFIX), function(f) {
        f <- file(paste0(LOG, "/", f))
        l <- readLines(f)
        close(f)
        
        sapply(c("reads processed",
                 "reads with at least one reported alignment",
                 "reads that failed to align",
                 "reads with alignments sampled due to -M"), function(x) {
                     gsub(".+: (.+)", "\\1", l[grep(x, l)])
                 })    
    })
    
    # set row and column names, and output the md table
    colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
    colnames(x) <- gsub(".cutadapt|.highQ|.deduped|.trimmed","", colnames(x))

    df <- data.frame(input_reads=format(as.numeric(x[1, ]), big.mark=","),
                     uniquely_mapped=paste(format(as.numeric(gsub("([[:digit:]]+) .+","\\1",x[2, ])),big.mark=","),gsub("([[:digit:]]+) (.+)","\\2",x[2, ])),
                     multi_mapped=paste(format(as.numeric(gsub("([[:digit:]]+) .+","\\1",x[4, ])),big.mark=","),gsub("([[:digit:]]+) (.+)","\\2",x[4, ])),
                     unmapped=paste(format(as.numeric(gsub("([[:digit:]]+) .+","\\1",x[3, ])),big.mark=","),gsub("([[:digit:]]+) (.+)","\\2",x[3, ])),
                     row.names=rownames(df))

    kable(df, align=c("r", "r", "r", "r"), output=F)
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
    bowtie.ref   <- shorten(system(paste("grep \"BOWTIE_REF\"",parameter.file,"| awk '{print $2;}'"), intern=TRUE))

    # combine into a data.frame
    bowtie.para.df <- data.frame(setting=c("bowtie flags","bowtie reference"),set_values=c(bowtie.flags,bowtie.ref))

    # print
    kable(bowtie.para.df, align=c("l", "r"), output=F)
}



##
## smallRNAhelper.bowtie.plot: plot bowtie mapping stats
##
smallRNAhelper.bowtie.plot <- function() {

    # plot folder
    if(!file.exists(SHINYREPS_BOWTIE_PLOT_DIR)) {
        return("Mapping statistics plot not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_BOWTIE_PLOT_DIR

    # construct the image url from the folder contents (skip current dir .)
    f <- list.files(SHINYREPS_BOWTIE_PLOT_DIR, pattern="^totalReads.png$")
    paste0("![](", QC, "/", basename(f), ")")

}



##
## smallRNAhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
smallRNAhelper.Fastqc <- function() {

    # logs folder
    if(!file.exists(SHINYREPS_FASTQC_OUT)) {
        return("Fastqc statistics not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_FASTQC_OUT

    # construct the image url from the folder contents (skip current dir .)
    samples <- list.dirs(SHINYREPS_FASTQC_OUT,recursive=F)[grep("cutadapt",list.dirs(SHINYREPS_FASTQC_OUT,recursive=F),invert=T)]
    df <- sapply(samples, function(f) {
        c(paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_quality.png)"),
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"))
    })

    # set row and column names, and output the md table
    df <- as.data.frame(t(df))
    rownames(df) <- gsub(paste0("^", SHINYREPS_PREFIX), "", basename(samples))
    rownames(df) <- gsub(paste0("_fastqc$"), "", rownames(df))
    colnames(df) <- c("Read qualities", "Sequence bias")
    kable(df, output=F, align="c")
}


##
## smallRNAhelper.filteredFastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
smallRNAhelper.filteredFastqc <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_FASTQC_OUT)) {
        return("Fastqc statistics not available")
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_FASTQC_OUT
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.dirs(SHINYREPS_FASTQC_OUT,recursive=F)[grep("cutadapt",list.dirs(SHINYREPS_FASTQC_OUT,recursive=F))]
    df <- sapply(samples, function(f) {
        c(paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_quality.png)"), 
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"),
          paste0("![fastqc img](", QC, "/", basename(f), "/Images/sequence_length_distribution.png)"))
    })

    # set row and column names, and output the md table
    df <- as.data.frame(t(df))
    rownames(df) <- gsub(paste0("^", SHINYREPS_PREFIX), "", basename(samples))
    rownames(df) <- gsub(paste0("_fastqc$"), "", rownames(df))
    rownames(df) <- gsub(".cutadapt|.highQ|.deduped|.trimmed","", rownames(df))
    colnames(df) <- c("Read qualities", "Sequence bias", "Length distribution")
    kable(df, output=F, align="c")
}


##
## smallRNAhelper.cutadapt: plot adapter trimming stats plot
##
smallRNAhelper.cutadapt <- function() {

    # logs folder
    if(!file.exists(SHINYREPS_CUTADAPT_PLOT_DIR)) {
        return("cutadapt statistics not available")
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_CUTADAPT_PLOT_DIR

    # construct the image url from the folder contents (skip current dir .)
    f <- list.files(SHINYREPS_CUTADAPT_PLOT_DIR, pattern="^trimmedReads.png$")
    paste0("![](", QC, "/", basename(f), ")")
}


##
## smallRNAhelper.qualityfilter: plot quatliy filter stats plot
##
smallRNAhelper.qualityfilter <- function() {

    # logs folder
    if(!file.exists(SHINYREPS_QUALITYFILTER_PLOT_DIR)) {
        return("Quality filter statistics not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_QUALITYFILTER_PLOT_DIR

    # construct the image url from the folder contents (skip current dir .)
    f <- list.files(SHINYREPS_QUALITYFILTER_PLOT_DIR, pattern="^qualityFilteredReads.png$")
    paste0("![](", QC, "/", basename(f), ")")
}


##
## smallRNAhelper.dedup: plot deduplication stats plot
##
smallRNAhelper.dedup <- function() {

    # logs folder
    if(!file.exists(SHINYREPS_DEDUP_PLOT_DIR)) {
        return("De-duplication statistics not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_DEDUP_PLOT_DIR

    # construct the image url from the folder contents (skip current dir .)
    f <- list.files(SHINYREPS_DEDUP_PLOT_DIR, pattern="^dedupReads.png$")
    paste0("![](", QC, "/", basename(f), ")")
}


##
## smallRNAhelper.rawfiltersummary: plot filtering stats summary plot
##
smallRNAhelper.rawfiltersummary <- function() {

    # logs folder
    if(!file.exists(SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR)) {
        return("Raw data filter summary statistics not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR

    # construct the image url from the folder contents (skip current dir .)
    f <- list.files(SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR, pattern="^allTrimmingStats.png$")
    paste0("![](", QC, "/", basename(f), ")")
}


##
## smallRNAhelper.fastqscreen: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
smallRNAhelper.fastqscreen <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_FASTQSCREEN_OUT)) {
        return("FastQScreen statistics not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_FASTQSCREEN_OUT

    SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){4})
    if(SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
    }
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_FASTQSCREEN_OUT, pattern=".png$", recursive=T, full.names=T)
    df <- sapply(samples, function(f) {
        paste0("![fastqscreen img](", f, ")")
    })
    
    # put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        gsub(".cutadapt|.highQ|.deduped|.trimmed|_screen.png)","", x)
    })

    df      <- matrix(df     , ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    samples <- matrix(samples, ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    colnames(df.names) <- rep(" ", SHINYREPS_PLOTS_COLUMN)
    
    kable(as.data.frame(df.names), align="c", output=F, format="markdown")
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
    
    # create a matrix using feature names as rownames, sample names as colnames
    l.infiles <- list.files(FOLDER, pattern=SUFFIX, full.names = T)
    l.counts <- lapply(l.infiles, function(f) {
        # f <- list.files(FOLDER, pattern=SUFFIX, full.names = T)[1]
        samplename <- basename(f)
        samplename <- gsub(SUFFIX, "", samplename)
	samplename <- gsub(".cutadapt|.highQ|.deduped|.trimmed", "", samplename)
        df.readcount <- read.table(f, header=T, sep='\t', comment.char = '#', col.names=c("Geneid","Length",samplename))
        df.readcount$Length <- NULL
	return(df.readcount)
   })
    
    # merge counts to data frame in wide format & remove "length..." columns
    f.merge <- function(x,y) {merge(x,y,by="Geneid")}
    df.counts <- Reduce(f.merge, l.counts)
    colnames(df.counts) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(df.counts))
    rm(l.counts)
    
    # eradicate biotype classes that are not present (or could be summarized as "other" being <1%)
    v.countsums <- rowSums(df.counts[-1])
    v.relsums <- v.countsums / sum(v.countsums)
    df.counts$Geneid[v.relsums < as.numeric(SHINYREPS_RNATYPES_CUTOFF)] <- "other"
    
    df.counts <- aggregate(df.counts[-1],
                           by=list(Geneid =df.counts$Geneid),
                           FUN=sum)
    
    # plot
    df.counts.melt <- melt(df.counts, id.var="Geneid")
    colnames(df.counts.melt) <- c("type","sample","count")
    
    plot <- ggplot() + 
        geom_bar(data=df.counts.melt, aes(x=sample, y=count, fill=type), position="fill", stat="identity") + 
        labs(x="", y="", fill="") + 
	theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) 
	#scale_fill_brewer(palette = "Set1")
    
    return(plot)
}



##
## smallRNAhelper.subread: parse Subread summary stats and create a md table
##
smallRNAhelper.subread <- function() {
    
    FOLDER <- SHINYREPS_SUBREAD
    SUFFIX <- paste0(SHINYREPS_SUBREAD_SUFFIX, '$')
    
    # check if folder exists
    if(!file.exists(FOLDER)) {
        return("Subread statistics not available")
    }
    
    # create a matrix using feature names as rownames, sample names as colnames
    x <- sapply(list.files(FOLDER, pattern=SUFFIX), function(f) {
        
        f <- file(paste0(FOLDER, '/', f))
        l <- readLines(f)
        close(f)
        
        
        sapply(c("Assigned",
                 "Unassigned_Ambiguity",
                 "Unassigned_MultiMapping",
                 "Unassigned_NoFeatures",
                 "Unassigned_Unmapped",
                 "Unassigned_MappingQuality",
                 "Unassigned_FragmentLength",
                 "Unassigned_Chimera",
                 "Unassigned_Secondary",
                 "Unassigned_Nonjunction",
                 "Unassigned_Duplicate"), function(y) {
                    as.numeric(  gsub( ".+\t(.+)", "\\1", l[grep(y, l)] )  )
                 })    
        
    })
    
    # correct column names
    colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
    colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
    colnames(x) <- gsub(".cutadapt|.highQ|.deduped|.trimmed","", colnames(x))
    
    # create md table (omitting various values that are 0 for now)
    #from x we romeove the ones which are unmapped to calculate percentages
    #only for the mapped ones
    x <- x[rownames(x) != "Unassigned_Unmapped", ]
    x <- rbind(total=x, colSums(x))
    rownames(x)[nrow(x)] <- "total"
    df <- data.frame(assigned=paste0(format(x[1, ], big.mark=","), " (", format((x[1, ]/x["total", ])*100, digits=2, nsmall=2), "%)"), 
                     unassigned_ambiguous=paste0(format(x[2, ], big.mark=","), " (", format((x[2, ]/x["total", ])*100, digits=2, nsmall=2), "%)"), 
                     unassigned_nofeature=paste0(format(x[4, ], big.mark=","), " (", format((x[4, ]/x["total", ])*100, digits=2, nsmall=2), "%)"))
    rownames(df) <- colnames(x)
    kable(df, align=c("r", "r", "r"), output=F)
    
}


##
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    table.content <- data.frame(  tool=c(
                                          "FastQC",
					  "cutadapt",  
                                          "Bowtie", 
                                          "samtools", 
                                          "Subread", 
                                          "FastQScreen"
                                          ), 
                                  version=c(
                                          Toolhelper.VersionReporter("FastQC",   SHINYREPS_FASTQC_LOG ), 
					  Toolhelper.VersionReporter("cutadapt", SHINYREPS_CUTADAPT_LOG ), 
                                          Toolhelper.VersionReporter("Bowtie",   SHINYREPS_BOWTIE_LOG ), 
                                          Toolhelper.VersionReporter("samtools", SHINYREPS_BAMINDEX_LOG), 
                                          Toolhelper.VersionReporter("Subread",  SHINYREPS_SUBREAD_LOG), 
                                          Toolhelper.VersionReporter("FastQScreen", SHINYREPS_FASTQSCREEN_LOG)
                                          )
                               )
    kable(table.content)
}

##
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    ver <- read.table(file=SHINYREPS_TOOL_VERSIONS,sep="=")
         ver$V1 <- strsplit(as.character(ver$V1),"_VERSION")
         colnames(ver) <- c("Tool name","Version")
               kable(as.data.frame(ver),output=F)
}
