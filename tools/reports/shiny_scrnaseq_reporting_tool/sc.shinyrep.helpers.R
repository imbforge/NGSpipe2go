##################################
##
## helper functions to create the plots for the Shiny report
##
##################################


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
#'         SHINYREPS_ORG <- "human"
#'         SHINYREPS_DB <- "hg38"
#'              
#' @description File content should be:
#'              SHINYREPS_PROJECT=projects/example_project
#'              SHINYREPS_ORG=human
#'              SHINYREPS_DB=hg38
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
    group.palette <- RColorBrewer::brewer.pal(9,"Set1")[1:num]
  } else {
    group.palette <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(num)
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
    replicate.palette <- RColorBrewer::brewer.pal(8,"Dark2")[1:num]
  } else {
    replicate.palette <- colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(num)
  }
  return(replicate.palette)
}




##
## DEhelper.STAR: parse STAR log files and create a md table
##
DEhelper.STARparms <- function() {
    
    # log file
    LOG <- SHINYREPS_STAR_LOG
    SUFFIX <- paste0(SHINYREPS_STARparms_SUFFIX, '$')
     if(!all(sapply(LOG, file.exists))) {
      return(paste("STAR statistics not available for", names(which(!sapply(LOG, file.exists)))))
    }
    
    # look for the lines containing the strings and get the values associated with this strings
    parseLog <- function(f) {
        # read in the lines
        #f <- file(paste0(LOG, "/", f))
        f <- file(f)
        l <- readLines(f)
        close(f)
        
        # get the version number from the first line (STAR svn revision compiled=STAR_2.3.1z13_r470)
        v <- unlist(strsplit(l[1], "="))[2]
        
        # get the redifined parameters and parse them in a key-value data.frame
        redefined <- l[grep("\\s+\\~RE-DEFINED$", l)]
        redefined <- sapply(redefined, function(x) {
            x <- unlist(strsplit(gsub("\\s+\\~RE-DEFINED$", "", x), "\\s+"))
            x[2] <- if(length(x) < 2) "" else x[2]
            x[2] <- shorten(x[2])
            x[c(1, 2)]
        })

        # take the 1st time a parm appears re-defined (skip re-defined in index generation)
        redefined <- redefined[, !duplicated(redefined[1, ])] 

        # create a vector with the values
        x <- redefined[2, ]
        names(x) <- redefined[1, ]
        
        # put the STAR version number
        x <- c(version=v, x)
        return(x)
    }
    df <- sapply(list.files(LOG, pattern=SUFFIX, full.names = T), parseLog)
    colnames(df) <- basename(colnames(df))
    
    # remove variable lines (lines depending on the fastq.gz file name)
    # and check if all the columns contain the same value. Display a warning otherwise
    df <- df[!grepl("(outFileNamePrefix|outTmpDir|readFilesIn)", rownames(df)), ]
    l <- apply(df, 1, function(x) length(unique(x)))    # rows differing (l > 1)
    df <- as.data.frame(df[, 1, drop=F])    # keep only the first column
    colnames(df) <- "parms"
    df$warning[l > 1] <- "Some files aligned with a different parm. Check logs"
    
    # set row and column names, and output the md table
    if(all(is.na(df$warning))) {
        #kable(df[, 1, drop=F], align=c("r"), output=F)
      kable(df[, 1, drop=F]) |> kable_styling()
    } else {
        #kable(df, align=c("r", "r"), output=F)
      kable(df) |> kable_styling()
      }
}

##
## DEhelper.STAR.violin: parse STAR log files and create violin plot and DataTable
##
DEhelper.STAR.violin <- function(colorByFactor=NULL, targetsdf=targets, ...) {
    
    # log file
    LOG <- SHINYREPS_STAR_LOG
    SUFFIX <- paste0(SHINYREPS_STAR_SUFFIX, '$')
    if(!all(sapply(LOG, file.exists))) {
      return(paste("STAR statistics not available for", names(which(!sapply(LOG, file.exists)))))
    }
    
    # look for the lines containing the strings
    # and get the values associated with this strings
    x <- list.files(LOG, pattern=SUFFIX, full.names = T)
  
    # select subset of samples 
    x <- selectSampleSubset(x, ...)

    x <- sapply(x, function(f) {
        f <- file(f)
        l <- readLines(f)
        close(f)
        
        sapply(c("Number of input reads",
                 "Uniquely mapped reads number",
                 "Uniquely mapped reads %",
                 "Number of reads mapped to multiple loci",
                 "Number of reads mapped to too many loci",
                 "% of reads mapped to multiple loci",
                 "% of reads mapped to too many loci",
                 "% of reads unmapped: too many mismatches",
                 "% of reads unmapped: too short",
                 "% of reads unmapped: other"), function(x) {
                     as.numeric(gsub("%", "", gsub(".+\\|\t(.+)", "\\1", l[grep(x, l)])))
                 })    
    })
    
    # set row and column names, and output the md table
    colnames(x) <- basename(colnames(x))
    colnames(x) <- gsub(Biobase::lcSuffix(colnames(x)), "", colnames(x)) # remove longest common suffix
    colnames(x) <- gsub(Biobase::lcPrefix(colnames(x)), "", colnames(x)) # remove longest common prefix
  
    df.stacked <- data.frame(filename = gsub("\\.R[12]\\.*$", "", colnames(x)),
                             input = x[1, ],
                             unique = x[2, ],
                             unique_perc = 100*(x[2, ]/x[1, ]),
                             multi = x[4, ],
                             multi_perc = 100*(x[4, ]/x[1, ]),
                             too_many_loci = x[7,],
                             mapped_perc = 100*((x[4, ]+x[2, ]) / x[1, ]),
                             unmapped = x[8, ] + x[9, ] + x[10, ])
    
## prepare groupwise plots
  if(!is.null(colorByFactor) && nrow(df.stacked) == nrow(targetsdf)) {
    # we want to plot the input reads and the mapped and the multi mapped reads numbers into different plots with different axises separated by one feature
    # we want to plot the same thing as percentages of the total
    # we have to plot per feature and then rearrange
    # we add one plot for the color value where we plot the percentages and color them according to the amount of input reads
  
    targetsdf$samplemod <- gsub(paste0(Biobase::lcSuffix(targetsdf$file ), "$"), "", targetsdf$file ) # shorten filename suffix
    #if(!is.na(SHINYREPS_PREFIX)) {targetsdf$samplemod  <- gsub(SHINYREPS_PREFIX, "", targetsdf$samplemod)}
    targetsdf$samplemod <- gsub(paste0("^", Biobase::lcPrefix(targetsdf$samplemod )), "", targetsdf$samplemod ) # shorten filename prefix
    
    index <- as.numeric(sapply(targetsdf$samplemod, function(x) grep(x, df.stacked$filename, ignore.case = T))) # grep for sample name in shortened file names
    if((nrow(df.stacked) != length(index)) || any(is.na(index))) {
      stop("\nThere seem to be ambiguous sample names in targets. Can't assign them uniquely to STAR logfile names")
    }
    
    if(any(!colorByFactor %in% colnames(targetsdf))) {
      if(all(!colorByFactor %in% colnames(targetsdf))) {
        cat("\nNone of the column names given in colorByFactor is available. Perhaps sample names are not part of fastq file names? Using filename instead.\n")
        colorByFactor <- "filename"
      } else { # one plot each element of colorByFactor
        cat("\n", colorByFactor[!colorByFactor %in% colnames(targetsdf)], "not available. Using", colorByFactor[colorByFactor %in% colnames(targetsdf)], "\n")
        colorByFactor <- colorByFactor[colorByFactor %in% colnames(targetsdf)]
      }
    }

    targetsdf$filename <- df.stacked$filename[index]
    df.stacked <- merge(df.stacked, targetsdf[,unique(c("filename", "sample", colorByFactor)), drop=F], by="filename")
    rownames(df.stacked) <- df.stacked$sample
    df.stacked <- df.stacked[order(rownames(df.stacked)),, drop=F]
    df.stacked <- df.stacked[,!apply(df.stacked,2, function(x) any(is.na(x))), drop=F] # remove NA columns from unsuccessful matching
    
  } else {
    # if colorByFactor == NULL or targets does not fit to number of files
    colorByFactor <- "filename"
  } # end groupwise plots
    
    
    # melt data frame for plotting
    df.melt  <- reshape2::melt(df.stacked, id.vars=unique(c(colorByFactor, "filename")), variable.name="map_feature")
    
    map.feature.plots <- lapply(colorByFactor, function(color.value){
      p <- ggplot(df.melt[df.melt$map_feature=="input",], aes_string("map_feature", "value", color=color.value)) +
        ggbeeswarm::geom_quasirandom() +
        scale_color_brewer(type= "qual", palette=2)  + 
        facet_wrap(~map_feature, scales="free") +
        scale_y_log10() +
        xlab(NULL) + theme(axis.text.x = element_text(size = 10)) +
        ylab("# Reads")
      p.perc <- ggplot(df.melt[grepl("perc", df.melt$map_feature),],
                       aes_string("map_feature", "value", color=color.value)) +
        ggbeeswarm::geom_quasirandom() +
        scale_color_brewer(type= "qual", palette=2)  + 
        facet_wrap(~map_feature, scales="free") +
        xlab(NULL) + theme(axis.text.x = element_text(size = 10)) +
        guides(color="none") +
        ylab("% of input reads")  
        # theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
      return(list(p.perc,p))
    })
    
    df.stats <- data.frame(sample=df.stacked$filename,
                     input_reads=format(df.stacked$input, big.mark=","), 
                     uniquely_mapped=paste0(format(df.stacked$unique, big.mark=","), " (", format(df.stacked$unique_perc, nsmall=2, digits=2), "%)"), 
                     multi_mapped=paste0(format(df.stacked$multi, big.mark=","), " (", format(df.stacked$multi_perc, nsmall=2, digits=2), "%)"), 
                     too_many_loci=paste0(format(df.stacked$too_many_loci,nsmall=2, digits=2), "%"),
                     unmapped=paste0(format(df.stacked$unmapped, nsmall=2, digits=2), "%")
    )
    
    stat=DT::datatable(df.stats, rownames=F, options = list(pageLength= 20))
    
    return(list(p_mapped=map.feature.plots, stat=stat))
}



##
## DEhelper.STAR: parse STAR log files and create a md table
##
DEhelper.STAR <- function(targetsdf=targets) {
  
  # log file
  LOG <- SHINYREPS_STAR_LOG
  SUFFIX <- paste0(SHINYREPS_STAR_SUFFIX, '$')
  if(!file.exists(LOG)) {
    return("STAR statistics not available")
  }
  
  # look for the lines containing the strings
  # and get the values associated with this strings
  x <- sapply(list.files(LOG, pattern = SUFFIX), function(f) {
    f <- file(paste0(LOG, "/", f))
    l <- readLines(f)
    close(f)
    
    sapply(c("Number of input reads",
             "Uniquely mapped reads number",
             "Uniquely mapped reads %",
             "Number of reads mapped to multiple loci",
             "% of reads mapped to multiple loci",
             "Number of reads mapped to too many loci",
             "% of reads mapped to too many loci",
             "% of reads unmapped: too many mismatches",
             "% of reads unmapped: too short",
             "% of reads unmapped: other"), function(x) {
               as.numeric(gsub("%", "", gsub(".+\\|\t(.+)", "\\1", l[grep(x, l)])))
             })    
  })
  
  # set row and column names, and output the md table
  colnames(x) <- gsub(paste0(SUFFIX, "$"), "", colnames(x))
  
  ## row.names need to be set in case of only one sample (if more than one sample,
  ## rownames would be set automatically, but doesnt harm to set them)
  df_values <- as.data.frame(t(x[1:7,]),row.names=colnames(x))
  
  df_values["unmapped"] <- x[1, ] - x[2, ] - x[4, ] - x[6,]
  df_values["% unmapped"] <- x[8, ] + x[9, ] + x[10, ]
  df_values$sample <- rownames(df_values)
  # we clean up the colnames a little to make them shorter and nicer
  colnames(df_values) <- gsub("of reads mapped to ", "",
                              gsub("Number of reads mapped to ", "",
                                   gsub("Uniquely mapped reads %","% unique",
                                        gsub("Uniquely mapped reads number","unique",
                                             gsub("Number of input reads","input",colnames(df_values))))))
  
  # if we have a differential expression analysis
  # we refactor the samples depending on group/subject or alternatively on the
  # amount of unique_mapping reads
  if(file.exists(SHINYREPS_TARGET)){
    
    #targets <- read.delim(SHINYREPS_TARGET, comment.char = "#")
    targets <- targetsdf
    targets$sample_ext <- gsub("\\..*$", "",targets$file )
    add_factors <- colnames(targets)[!colnames(targets) %in% c("group", "sample", "file")]
    
    # replace files names with nicer sample names given in targets file 
    # if sample is missing in targets file, use reduced file name
    df_values$sample <- sapply(df_values$sample, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1 && !any(duplicated(targets$sample)),   
                                                                      targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                                      gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
  } else{
    
    # remove sample prefix from sample names (suffix was already removed before)
    df_values$sample <- gsub(paste0("^",SHINYREPS_PREFIX), "", df_values$sample)
  }
  
  # sort sample alphabetically or based on number of (uniquely) mapped reads
  if (SHINYREPS_SORT_ALPHA) {
    df_values$sample <- fct_reorder(df_values$sample, df_values$sample, .desc=TRUE)
  } else {
    df_values$sample <- fct_reorder(df_values$sample, df_values$`% unique`)
  }
  
  df_melt <- reshape2::melt(df_values, value.name = "reads", variable.name = "mapping_stat")
  df_melt$value_info <- ifelse(grepl("%", df_melt$mapping_stat), "perc", "reads")
  
  ## melt does not work properly in case of only one sample (sample name gets lost)
  if (nrow(df_values) == 1) {
    df_melt$sample <- rownames(df_values)
  }
  
  ## invert order for plotting
  df_melt$mapping_stat <- fct_rev(df_melt$mapping_stat)
  
  ## define color scheme
  my.palette <- rev(c("#e08a28","#dfc27d","#acd2f7","#5ea1e3"))
  my.color   <- "#cc0e25"
  
  # plot showing percent mapped reads
  p_perc <- ggplot(df_melt[df_melt$value_info == "perc",], 
                   aes(x = sample, y = reads, fill = mapping_stat )) +
    geom_bar(stat     = "identity", 
             position = "stack") +
    labs(x = "", 
         y = "% of sequenced reads", 
         title = "Mapping summary (percentage)") + 
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5,size=12),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) +
    scale_fill_manual(values=my.palette) +
    coord_flip()
  
  # max number of reads in any sample
  max.reads <- max(df_melt[df_melt$mapping_stat=="input","reads"])
  
  # plot a combination of percent mapped reads and total sequenced reads
  p_perc_count <- ggplot() + 
    geom_bar(data=df_melt[df_melt$value_info == "reads" & df_melt$mapping_stat!="input",], 
             mapping=aes(x = sample, y = reads, fill = mapping_stat), 
             stat = "identity", position = "fill", width = 0.8) + 
    geom_point(data=df_melt[df_melt$mapping_stat=="input",], 
               mapping=aes(x = sample, y = reads/max.reads), 
               size=2, fill=my.color, color=my.color, shape=18) +
    geom_line(data=df_melt[df_melt$mapping_stat=="input",],
              mapping=aes(x = sample, y = reads/max.reads, group=1),
              color=my.color, linetype="dashed", size=0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ . *max.reads, name="# sequenced reads"),
                       labels = scales::percent_format()) + 
    labs(x = "",
         y = "% of sequenced reads") + 
    theme(axis.text.y = element_text(size = 8),
          axis.title.x.top=element_text(color=my.color), 
          axis.text.x.top=element_text(color=my.color), 
          axis.ticks.x.top=element_line(color=my.color),
          plot.title = element_text(hjust = 0.5,size=12),
          legend.position = "top") + 
    guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) + 
    scale_fill_manual(values=my.palette) + 
    coord_flip()
  
  # in case of non-alphabetically sorting, re-order based on total number of reads
  # before plotting reads counts
  if (SHINYREPS_SORT_ALPHA == FALSE) {
    df_melt$sample <- factor(df_melt$sample, levels=df_values$sample[order(df_values$input)])
  }
  
  p_count <- p_perc %+%
    df_melt[df_melt$value_info == "reads" & df_melt$mapping_stat != "input",] +
    labs(x = "",
         y = "# sequenced reads",
         title = "Mapping summary (read counts)") 
  
  # prepare data for summary table  
  rownames(df_values) <- df_values$sample
  df_values <- df_values[, colnames(df_values) != "sample"]
  
  # reformat individual columns
  df_values[, grepl("%", colnames(df_values))] <- as.data.frame(
    lapply(
      df_values[, grepl("%", colnames(df_values))], function(x){
        paste0(format(x, nsmall=2), "%") 
      }))
  df_values[, !grepl("%", colnames(df_values))] <- as.data.frame(
    lapply(
      df_values[, !grepl("%", colnames(df_values))], function(x){
        format(x, big.mark=",") 
      }))
  
  # create data frame to print as table
  df_values_print <- data.frame(all.reads     = df_values$input,
                                unique        = paste0(df_values[,"unique"]," (",df_values[,"% unique"],")"),
                                multiple.loci = paste0(df_values[,"multiple loci"]," (",df_values[,"% multiple loci"],")"),
                                too.many.loci = paste0(df_values[,"too many loci"]," (",df_values[,"% too many loci"],")"),
                                unmapped      = paste0(df_values[,"unmapped"]," (",df_values[,"% unmapped"],")"), 
                                row.names     = rownames(df_values))
  
  # sort alphabetically based on sample names (=rownames) or based on % unique
  if (SHINYREPS_SORT_ALPHA == TRUE) {
    df_values_print <- df_values_print[order(rownames(df_values_print)),]
  } else {
    df_values_print <- df_values_print[rownames(df_values[order(df_values$`% unique`, decreasing=TRUE),]),]
  }
  colnames(df_values_print) <- gsub("\\."," ",colnames(df_values_print))
  
  return( list(p_perc = p_perc,
               p_perc_count = p_perc_count,
               p_count = p_count,
               stat= DT::datatable(df_values_print, options = list(pageLength= 20))
               #stat = kable(df_values_print, align=c("r", "r", "r", "r", "r"), format="markdown", output=F)
               )
  )
  
}



##
## DEhelper.Fastqc.custom: prepare Fastqc summary plots
##
DEhelper.Fastqc.custom <- function(web=FALSE, summarizedPlots=TRUE, targetsdf=targets, subdir="", ...) {
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
  names(lbls) <- gsub("_fastqc.zip", ".fastq.gz", basename(names(fastqc.stats)))
  
  if(file.exists(SHINYREPS_TARGET)){ 
    # get target names
    #targets <- read.delim(SHINYREPS_TARGET, comment.char = "#")
    targets <- targetsdf
    targets$sample_ext <- gsub("\\..*$", "",targets$file )

    # replace files names with nicer sample names given in targets file 
    # if sample is missing in targets file, use reduced file name
    #cat('In if - before gsub')
    lbls <- sapply(lbls, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1 && !any(duplicated(targets$sample)),
                                              targets[sapply(targets$sample_ext, grepl, i),"sample"],
                                              gsub("_001.fastq.gz", "", gsub(paste0("^",SHINYREPS_PREFIX),"",basename(i))))})
    
    if(any(grepl("[\\._]R1[\\._]|[\\._]R2[\\._]", names(lbls)))) {  # SHINYREPS_PAIRED == "yes" # for scRNA-Seq we may have se mapping but still an R2 with barcode
      x <- names(lbls)
      lbls <- paste0(lbls, ifelse(grepl("[_\\.]R1", names(lbls)), "_R1", "_R2"))
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
  return(list(no.of.samples=length(f), p.qual=p.qual, p.content=p.content, p.gc=p.gc))
}




#' DEhelper.fastqscreen: summarizes FastQScreen results, creates summarized barplots, only relevant contanimants shown
#'
#' @param perc.to.plot - a numeric vector of length 1 setting the percent cutoff of relevant contaminants, if any sample
#'                       shows more than perc.to.plot, contaminant will be shown in plot
#'
#' @return a list including a plot, the number of samples, and the number of plotted contaminants
#'
DEhelper.fastqscreen <- function(subdir="", targetsdf=targets, perc.to.plot = 1, ncol=2, ...) {
  
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
  samples <- samples[sapply(samples, function(x) {any(grepl(".*screen.html", list.files(x)))})] # exclude incomplete results
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
    #targets <- read.delim(SHINYREPS_TARGET, comment.char = "#")
    targets <- targetsdf
    targets$sample_ext <- gsub("\\..*$", "",targets$file )

    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    samples <- sapply(samples, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1 && !any(duplicated(targets$sample)),
                                                    ifelse(sapply("\\.R1", grepl, i), 
                                                           paste0(targets[sapply(targets$sample_ext, grepl, i),"sample"], ".R1"),
                                                           ifelse(sapply("\\.R2", grepl, i), 
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
    facet_wrap(~genome,ncol=ncol) +
    coord_flip()
  
  return(list(p.category.wrap=p.category.wrap,
              no.of.genomes=length(unique(df$genome)),
              no.of.samples=length(unique(df$sample)),
              no.of.rows = ceiling(length(unique(df$genome))/ncol)
              )
        )      
}




##
## DEhelper.RNAtypes: parse Subread count results for RNAtypes usage
##
DEhelper.RNAtypes <- function(targetsdf=SHINYREPS_TARGET, ...) {
  
  FOLDER <- SHINYREPS_RNATYPES
  SUFFIX <- paste0(SHINYREPS_RNATYPES_SUFFIX, '$')
  
  # check if folder exists
  if(!file.exists(FOLDER) || length(list.files(FOLDER))==0) {
    return("Subread statistics not available")
  }
  
  # list files to plot
  plotfiles <- list.files(FOLDER, pattern=SUFFIX, full.names=T)
  
  # select subset of plotfiles in case there are too many separate files (select all by default)
  plotfiles <- selectSampleSubset(plotfiles, grepInBasename=T, ...)
  
  names(plotfiles)  <- gsub(SUFFIX, "", basename(plotfiles))
  names(plotfiles)  <-  gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", names(plotfiles))
  
  l <- lapply(1:length(plotfiles), function(i) {
    df.readcount <- read.table(plotfiles[i], header=T, sep='\t', comment.char = '#', col.names=c("Geneid","Length",names(plotfiles)[i]))
    df.readcount$Length <- NULL
    return(df.readcount)
  })

  # merge counts to data frame in wide format & remove "length..." columns
  f.merge <- function(x,y) {merge(x,y,by="Geneid")}
  df.counts <- Reduce(f.merge, l)

  # eradicate biotype classes that are not present (or could be summarized as "other" being <1%)
  v.countsums <- rowSums(df.counts[-1])
  v.relsums <- v.countsums / sum(v.countsums)
  df.counts$Geneid[v.relsums < 0.01] <- "other"
  
  df.counts <- aggregate(df.counts[-1],
                         by=list(Geneid =df.counts$Geneid),
                         FUN=sum)
  
  df <- reshape2::melt(df.counts, id.var="Geneid")
  colnames(df) <- c("type","plotfile","count")
  
  if(class(targetsdf)=="data.frame" || file.exists(targetsdf)){
    
    # get targets
    if(class(targetsdf)=="data.frame") {
      targets <- targetsdf
    } else {
      targets <- read.delim(targetsdf, comment.char = "#")
    }
    
    targets$file <- gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", targets$file)
    
    # now left_join by file name
    df <- dplyr::left_join(df, targets, by=dplyr::join_by(plotfile==file))
    
    # use sample column as plotfile if no multiple sample names per plotfile present
    # and sample names are shorter than file names (e.g. for ParseBio the sample column 
    # may be longer if it consist of multiple concatenated sample names)
    if("sample" %in% colnames(df) && all(!is.na(df$sample)) && 
       identical(df[!duplicated(df$sample),], df[!duplicated(df$plotfile),]) && 
       max(nchar(df$sample)) < max(nchar(df$plotfile))) {
      df$plotfile <- df$sample
    }
  }   
  
  # prune prefix from plotfile if present
  if(!is.na(SHINYREPS_PREFIX)) {
    df$plotfile <- gsub(paste0("^",SHINYREPS_PREFIX), "", df$plotfile)
  } 
  
  # remove possible starting "X" (in case sample name starts with a number)
  df$plotfile <- gsub("^X","",df$plotfile)
  
  # check if targets.txt info is merged properly, factorize and sort.
  df$plotfile <- factor(df$plotfile, levels = rev(sort(unique(df$plotfile))))
  if("group" %in% colnames(df) && all(!is.na(df$group)) && SHINYREPS_SORT_ALPHA==FALSE) {
    df$group <- factor(df$group, levels = sort(unique(df$group)))
    df <- df[order(df$group, df$plotfile),]
  } else {
    df <- df[order(df$plotfile),]
  }

  plot <- ggplot() +
    geom_bar(data=df, aes(x=plotfile, y=count, fill=type), position="fill", stat="identity") +
    labs(x="", y="", fill="") +
    theme_bw() +
    theme(axis.text.y = element_text(size=8), legend.position = "top") +
    guides(fill = guide_legend(nrow=1, title="", reverse=TRUE)) +
    scale_fill_brewer(palette = "Dark2") +
    coord_flip()
  
  return(plot)
}




##
## DEhelper.geneBodyCov2: go through geneBodyCov output dir and plot into one plot
##
DEhelper.geneBodyCov2 <- function(web=F, targetsdf=SHINYREPS_TARGET, ...) {
  
  # logs folder
  if(!file.exists(SHINYREPS_GENEBODYCOV_LOG) || length(list.files(SHINYREPS_GENEBODYCOV_LOG))==0) {
    return("geneBodyCov statistics not available")
  }
  
  # construct the folder name, which is different for web and noweb
  QC <- if(web) "/geneBodyCov" else SHINYREPS_GENEBODYCOV_LOG
  
  # list files to plot
  plotfiles <- list.files(QC, pattern=".csv$", full.names=T)
  
  # select subset of plotfiles in case there are too many separate files (select all by default)
  plotfiles <- selectSampleSubset(plotfiles, grepInBasename=T, ...)
  
  names(plotfiles) <- gsub("_geneBodyCov.csv", "", basename(plotfiles)) # prune file names
  names(plotfiles) <-  gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", names(plotfiles))
  
  l <- lapply(1:length(plotfiles), function(i) {
    df <- read.csv(plotfiles[[i]]) 
    colnames(df) <- c("perc", names(plotfiles)[i])
    return(df[,2, drop=F])
  })
  df <- do.call(cbind, l) # coverage df with files as columns percentage as rows
  
  df$perc <- 1:nrow(df) #these are the 100 bins used in the original plot
  df <- reshape2::melt(df, id.vars="perc", variable.name="plotfile", value.name="cov")
  

  if(class(targetsdf)=="data.frame" || file.exists(targetsdf)){
    
    # get targets
    if(class(targetsdf)=="data.frame") {
      targets <- targetsdf
    } else {
      targets <- read.delim(targetsdf, comment.char = "#")
    }
    
    targets$file <- gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", targets$file)

    # now left_join by file name
    df <- dplyr::left_join(df, targets, by=dplyr::join_by(plotfile==file))
    
    # use sample column as plotfile if no multiple sample names per plotfile present
    # and sample names are shorter than file names (e.g. for ParseBio the sample column 
    # may be longer if it consist of multiple concatenated sample names)
    if("sample" %in% colnames(df) && all(!is.na(df$sample)) && 
       identical(df[!duplicated(df$sample),], df[!duplicated(df$plotfile),]) && 
       max(nchar(df$sample)) <= max(nchar(df$plotfile))) {
      df$plotfile <- df$sample
    }
  }   
 
  # prune prefix from plotfile if present
  if(!is.na(SHINYREPS_PREFIX)) {
         df$plotfile <- gsub(paste0("^",SHINYREPS_PREFIX), "", df$plotfile)
  } 
  
  # check if targets.txt info is merged properly and factorize. Remove column if existing but NAs present.
  df$plotfile <- factor(df$plotfile, levels = sort(unique(df$plotfile)))
  if("group" %in% colnames(df) && all(!is.na(df$group))) {
    df$group <- factor(df$group, levels = sort(unique(df$group)))
  } else {df$group <- NULL}
  if("replicate" %in% colnames(df) && all(!is.na(df$replicate))) {
    df$replicate <- factor(df$replicate, levels = sort(unique(df$replicate)))
  } else {df$replicate <- NULL}
  if("cells" %in% colnames(df) && all(!is.na(df$cells))) {
    pal_control_wells <- c("1c" = "grey70", "0c"  = "grey50", "10c" = "grey0")
    df$cells <- factor(df$cells, levels = names(pal_control_wells))
    df <- df[order(df$cells),]
   } else {df$cells <- NULL}
  
    
  # plotly plot per file
  plot_df <- plotly::highlight_key(df, ~plotfile) 
  
  if("cells" %in% colnames(df)) {
    p <- ggplot(plot_df, aes(x=perc, y=cov, group=plotfile, color= cells )) + geom_line() +
      scale_color_manual(values = pal_control_wells) + labs(color = "") 
  } else {
    p <- ggplot(plot_df, aes(x=perc, y=cov, group=plotfile)) + geom_line(color="darkgrey") 
  }
  p <- p + labs(title="",
                x = "Gene body percentile 5' -> 3'",
                y = "Averaged normalised coverage") +
    ylim(0,1) + theme_bw() 
  
  gg <- plotly::ggplotly(p) 

  plot_list <- list()
  plot_list[["plotly"]] <- gg
  
  
  # create static plots as well
  num_samples <- length(unique(df$plotfile))
  if(num_samples < 10){ #we only have 9 colors Set1
    sample.cols <- RColorBrewer::brewer.pal(num_samples, "Set1")
  }else{
    sample.cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_samples)
  }

  if(all(c("group","replicate") %in% colnames(df)) && length(unique(df$group))>1){
    num_groups <- length(unique(df$group))
    num_replicates <- length(unique(df$replicate))
    plot_list[["num_groups"]] <- num_groups
    plot_list[["num_replicates"]] <- num_replicates
    
    p_per_replicate_splitByGroup <- ggplot(df, aes(x=perc,
                                                   y=cov,
                                                   group=plotfile,
                                                   color=replicate)) +   
      geom_line() +
      labs(x = "Gene body percentile 5' -> 3'",
           y = "Averaged normalised coverage") +
      scale_color_manual(values=define.replicate.palette(num_replicates)) + 
      ylim(0,1) +
      theme_bw() + 
      theme(legend.position="top") +
      guides(color=guide_legend(ncol=6)) +
      facet_wrap(~group,ncol=2)
    
    plot_list[["p_per_replicate_splitByGroup"]] <- p_per_replicate_splitByGroup
    
    if(length(unique(df$replicate))>1) {
    
      p_per_group_splitByReplicate <- ggplot(df, aes(x=perc,
                                                     y=cov,
                                                     group=plotfile,
                                                     color=group)) +
        geom_line() +
        labs(x = "Gene body percentile 5' -> 3'",
             y = "Averaged normalised coverage") +
        scale_color_manual(values=define.group.palette(num_groups)) +
        ylim(0,1) +
        theme_bw() +
        theme(legend.position="top") +
        guides(color=guide_legend(ncol=3)) +
        facet_wrap(~replicate,ncol=2)
    
    plot_list[["p_per_group_splitByReplicate"]] <- p_per_group_splitByReplicate
    }
 
  }
  return(plot_list)
}




##
##DEhelper.strandspecificity: get the strand specificity from the qc and display them
##
DEhelper.strandspecificity <- function(targetsdf=SHINYREPS_TARGET, ...){
  
  # logs folder
  if(!file.exists(SHINYREPS_INFEREXPERIMENT_LOGS)) {
    return("Strand specificity statistics not available")
  }
  
  filelist <- list.files(path=SHINYREPS_INFEREXPERIMENT_LOGS, full.names=TRUE)
  filelist <- selectSampleSubset(filelist, ...)
  strandspecificity <- tryCatch(lapply(filelist, read.table, sep=":", skip=3, header=FALSE, row.names=1, blank.lines.skip=TRUE), 
                                error=function(e) NULL)
  if (is.null(strandspecificity)) {
    return("Strand specificity statistics not readable")  # exits the function
  }
  
  strandspecifity <- do.call(cbind, strandspecifity)
  samplenames <- gsub("_inferexperiment.txt", "", basename(filelist))
  
  if(class(targetsdf)=="data.frame" || file.exists(targetsdf)){
    
    # get target names
    if(class(targetsdf)=="data.frame") {
      targets <- targetsdf
    } else {
      targets <- read.delim(targetsdf, comment.char = "#")
    }
    targets$sample_ext <- gsub("\\..*$", "",targets$file )

    # replace files names with nicer sample names given in targets file
    # if sample is missing in targets file, use reduced file name
    samplenames <- sapply(samplenames, function(i) { ifelse(sum(sapply(targets$sample_ext, grepl, i))==1 && !any(duplicated(targets$sample)),   
                                                            targets[sapply(targets$sample_ext, grepl, i),"sample"], 
                                                            gsub(paste0("^",SHINYREPS_PREFIX),"",i))})
  } else {
    
    if(!is.na(SHINYREPS_PREFIX)) {
      samplenames <- gsub(paste0("^",SHINYREPS_PREFIX), "", samplenames)
    }
  }
  
  colnames(strandspecifity) <- samplenames 
  rownames(strandspecifity) <- c("ambiguous", "sense", "antisense")
  
  # sort alphabetically for report
  strandspecifity <- as.data.frame(t(strandspecifity)[order(rownames(t(strandspecifity))),,drop=F])
  
  # add another column specifying the strandedness based on which strandedness is expected [no|yes|reverse]
  if (SHINYREPS_STRANDEDNESS == "yes") {
    strandspecifity$`strandedness [%]` <- 100*round(strandspecifity$sense/(strandspecifity$sense+strandspecifity$antisense),digits=4)
  } else {
    if (SHINYREPS_STRANDEDNESS == "reverse") {
      strandspecifity$`strandedness [%]` <- 100*round(strandspecifity$antisense/(strandspecifity$sense+strandspecifity$antisense),digits=4)
    }
  }
  
  # kable(strandspecifity, output=F, format="markdown", align=c("c"))
  DT::datatable(strandspecifity, options = list(pageLength= 20))
}


##
## DEhelper.cutadapt: get trimming statistics from the Cutadapt folder and display them
## 
#' @param targetsdf targets data.frame or character with file path to targets object
#' @param colorByFactor character with column name of sample table to be used for coloring the plot. Coloring by filename if NULL. 
#' @param sampleColumnName character with column name(s) of targets table containing file names
#' @param plotfun define function to be used for plotting
#' @param labelOutliers logical, shall outlier samples be labeled
#' @param outlierIQRfactor numeric, factor is multiplied by IQR to determine outlier
#'
#' @return plot cutadapt statistics as side effect
DEhelper.cutadapt <- function(targetsdf=targets, colorByFactor="group", sampleColumnName =c("file"), 
                              plotfun=DEhelper.cutadapt.plot, labelOutliers=T, outlierIQRfactor=1.5, 
                              ...){
  
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
    if(is.na(SHINYREPS_PREFIX)) {row.names(x.df)  <- gsub(Biobase::lcPrefix(row.names(x.df) ), "", row.names(x.df) )}
    row.names(x.df)  <- gsub(Biobase::lcSuffix(row.names(x.df) ), "", row.names(x.df) )
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
    
    x.df <- data.frame(x.df[unlist(t(index)),], targetsdf, check.names =F)
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    if("sample" %in% colnames(x.df) && !any(duplicated(x.df$sample))) { # use sample column as identifier if present and unique
      x.df$filename <- x.df$sample
      row.names(x.df) <- x.df$sample} 
    
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


# plotting function for DEhelper.cutadapt 
DEhelper.cutadapt.plot <- function(data, color.value, labelOutliers=T, outlierIQRfactor=1.5){
  
  is_outlier <- function(x) { # function for identification of outlier
    if(IQR(x)!=0) {
    return(x < quantile(x, 0.25) - outlierIQRfactor * IQR(x) | x > quantile(x, 0.75) + outlierIQRfactor * IQR(x))
    } else {
    return(x < mean(x) - outlierIQRfactor * mean(x) | x > mean(x) + outlierIQRfactor * mean(x))
    }
  }
  
  data <- data |>
    dplyr::group_by(reads) |>
    dplyr::mutate(outlier=is_outlier(value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(outlier=ifelse(outlier,filename,as.numeric(NA))) |>
    as.data.frame()
  
  ylab <- "% reads"
  
  # prepare palette of appropriate length according to the different factors given in colorByFactor
  colourCount = length(unique(data[,color.value]))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

  names(data)[names(data)==color.value] <- "color.value" 
 
  p <- ggplot(data, aes(x=reads,
                        y=value,
                        color=color.value ))+
    geom_quasirandom(groupOnX=TRUE) +
    geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA)  

    if(labelOutliers) {p <- p + ggrepel::geom_text_repel(data=. |> filter(!is.na(outlier)), aes(label=filename), show.legend=F)}

  p <- p + scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
           ylab(ylab) +
           xlab("") +
           theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1)) +
           guides(color=guide_legend(title=color.value))
  
  return(p)
}



##
##DEhelper.umicount: get deduplication stats from UMI_tools count
## 
DEhelper.umicount <- function(colorByFactor=NULL, targetsdf=targets, ...){
  
  
  x <- list.files(SHINYREPS_UMICOUNT_LOG,pattern='*umicount.log$',full.names=TRUE) 
  # select subset of samples for fastqc figures (e.g. merged singlecell pools) or use all samples for samplePattern=NULL
  
  x <- selectSampleSubset(x, ...)
  if(length(x) == 0) {
      return("No samples matched with this pattern...")
  }
  
  x <- sapply(x, function(f) { 
    input_reads_total <- system(paste("grep \"INFO Input Reads\"", f, "| awk '{print $6}'"), intern=TRUE)
    
    skipped_reads_total <- system(paste("grep \"INFO Read skipped, no tag\"", f, "| awk '{print $8}'"), intern=TRUE)
    
    counted_reads_total <- system(paste("grep \"INFO Number of (post deduplication) reads counted\"", f, "| awk '{print $10}'"), intern=TRUE)
    
    return(c(unique(input_reads_total), unique(skipped_reads_total), unique(counted_reads_total)))
  })
  
  # set row and column names
  x.df <- as.data.frame(t(x)) 
  colnames(x.df) <- c("input_reads_total", "skipped_reads_total","counted_reads_total")
  x.df <- as.data.frame(lapply(x.df, as.numeric))
  x.df$skipped_reads <- round(100* (x.df$skipped_reads_total / x.df$input_reads_total ), 2)
  x.df$counted_reads <- round(100* (x.df$counted_reads_total / x.df$input_reads_total), 2)
  
  
  # reduce size of file names 
  row.names(x.df) <- basename(colnames(x))
  row.names(x.df)  <- gsub(Biobase::lcSuffix(row.names(x.df) ), "", row.names(x.df) )
  row.names(x.df)  <- gsub(Biobase::lcPrefix(row.names(x.df) ), "", row.names(x.df) )
  #if(!is.na(SHINYREPS_PREFIX)) {row.names(x.df) <- gsub(SHINYREPS_PREFIX, "", row.names(x.df))}
  x.df$filename <- factor(row.names(x.df))
  
  
  # passing the different factors given in targetsdf to x.df which was created from cutadapt logfile names (if 1 cell per file)
  if(!is.null(colorByFactor) && nrow(x.df) == nrow(targetsdf)) { # if targets object fits in length, add information to x.df
    
    
    targetsdf$samplemod <- gsub(Biobase::lcSuffix(targetsdf$file ), "", targetsdf$file ) # shorten filename suffix
    #if(!is.na(SHINYREPS_PREFIX)) {targetsdf$samplemod  <- gsub(SHINYREPS_PREFIX, "", targetsdf$samplemod)}
    targetsdf$samplemod <- gsub(Biobase::lcPrefix(targetsdf$samplemod ), "", targetsdf$samplemod ) # shorten filename prefix
    
    
    index <- as.numeric(sapply(targetsdf$samplemod, function(x) grep(x, x.df$filename, ignore.case = T))) # grep for sample name in shortened file names
    if(nrow(x.df) != length(index) || any(is.na(index))) {
      stop("\nThere seem to be ambiguous sample names in targets. Can't assign them uniquely to cutadapt logfile names")
    }
    
    targetsdf$filename <- x.df$filename[index]
    x.df <- merge(x.df, targetsdf, by="filename")
    x.df <- x.df[order(rownames(x.df)),, drop=F]
    x.df <- x.df[,!apply(x.df,2, function(x) any(is.na(x))), drop=F] # remove NA columns from unsuccessful matching
    if("sample" %in% colnames(x.df)) {x.df$filename <- x.df$sample} # use sample column as identifier if present
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
    # if colorByFactor == NULL or targets does not fit to number of files
    colorByFactor <- "filename"
  }
  
  # melt data frame for plotting
  x.melt <- reshape2::melt(x.df, measure.vars=c("skipped_reads", "counted_reads"),
                 variable.name="reads")
  #everything which is not a value should be a factor
  
  #now we do a violin plot of the trimmed/too_short/etc. ones and color it
  # according to the different factors given in colorByFactor 
  
  # prepare palette of appropriate length
  colourCount = length(unique(x.melt[,colorByFactor]))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  create.violin <- function(x.melt, color.value){
    ylab <- "% reads"
    p <- ggplot(x.melt, aes_string(x="reads",
                                   y="value",
                                   color=color.value ))+
      geom_quasirandom() +
      scale_color_manual(values=getPalette(colourCount)) + # creates as many colors as needed
      ylab(ylab) +
      xlab("") +
      scale_y_continuous( breaks=seq(0, max(x.melt$value), 10),
                          limits = c(0, max(x.melt$value))) 
    return(p)
  }
  
  # one plot each element of colorByFactor
  violin.list <- lapply(colorByFactor, create.violin, x.melt=x.melt) # "colorByFactor" submitted as color.value
  
  for(i in 1:length(violin.list)){
    plot(violin.list[[i]])
  }
  
  #kable(x.df[,c("total.reads", "trimmed","tooshort")], output=F, format="markdown", align=c("l")) |> kable_styling()
  colnames(x.df)[colnames(x.df)=="skipped_reads"] <- "skipped_reads_perc"
  colnames(x.df)[colnames(x.df)=="counted_reads"] <- "counted_reads_perc"
  DT::datatable(x.df[,c("input_reads_total", "skipped_reads_total","skipped_reads_perc", "counted_reads_total", "counted_reads_perc")])
}



##
## extract the intron/exon and intergenic regions from the qualimap report
##
DEhelper.Qualimap <- function(targetsdf=SHINYREPS_TARGET, ...) {
  
  # logs folder
  if(!file.exists(SHINYREPS_QUALIMAP_LOGS) || length(list.files(SHINYREPS_QUALIMAP_LOGS))==0) {
    return("Qualimap reports not available")
  }  
  
  QC <- SHINYREPS_QUALIMAP_LOGS    
  
  # list files to plot
  plotfiles <- list.files(QC, pattern="rnaseq_qc_results.txt$", recursive=T, full.names=T)
  
  # select subset of plotfiles in case there are too many separate files (select all by default)
  plotfiles <- selectSampleSubset(plotfiles, grepInBasename=T, ...)
  
  names(plotfiles)  <- gsub("_counts_qualimap", "", basename(dirname(plotfiles)))
  names(plotfiles)  <-  gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", names(plotfiles))
  
  plotfiles <- sapply(plotfiles, function(f){
    l <- readLines(f)
    l <- l[grep("exonic|intronic|intergenic", l)] 
    l <- strsplit(gsub(" ", "", l), "=")
    names(l) <- c("exonic","intronic","intergenic")
    l <- sapply(l, function(x){
      as.numeric(gsub("\\%\\)", "", strsplit(x[[2]], "\\(")[[1]][2]))    
    })
  })

  # qualimap outputs an additional category "overlapping_exon", which is part of intron or intergenic,
  # but since it is not specified, which of the two, we will not display this category here.

  # sort according to the groups and samplenames in the targets or alphabetically
  df <- reshape2::melt(plotfiles, value.name="perc")
  colnames(df) <- c("class", "plotfile", "perc")
  
  if(class(targetsdf)=="data.frame" || file.exists(targetsdf)){
    
    # get targets
    if(class(targetsdf)=="data.frame") {
      targets <- targetsdf
    } else {
      targets <- read.delim(targetsdf, comment.char = "#")
    }
    
    targets$file <- gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", targets$file)
    
    # now left_join by file name
    df <- dplyr::left_join(df, targets, by=dplyr::join_by(plotfile==file))
    
    # use sample column as plotfile if no multiple sample names per plotfile present
    # and sample names are shorter than file names (e.g. for ParseBio the sample column 
    # may be longer if it consist of multiple concatenated sample names)
    if("sample" %in% colnames(df) && all(!is.na(df$sample)) && 
       identical(df[!duplicated(df$sample),], df[!duplicated(df$plotfile),]) && 
       max(nchar(df$sample)) <= max(nchar(df$plotfile))) {
      df$plotfile <- df$sample
    }
  }   
  
  # prune prefix from plotfile if present
  if(!is.na(SHINYREPS_PREFIX)) {
    df$plotfile <- gsub(paste0("^",SHINYREPS_PREFIX), "", df$plotfile)
  } 

  # check if targets.txt info is merged properly, factorize and sort.
  df$plotfile <- factor(df$plotfile, levels = rev(sort(unique(df$plotfile))))
  df$class <- factor(df$class, levels=c("intergenic", "intronic", "exonic" )) 
  if("group" %in% colnames(df) && all(!is.na(df$group)) && SHINYREPS_SORT_ALPHA==FALSE) {
    df$group <- factor(df$group, levels = sort(unique(df$group)))
    df <- df[order(df$group, df$plotfile),]
  } else {
    df <- df[order(df$plotfile),]
    }
 
  # plot
  p <- ggplot(df, aes(x=plotfile,y=perc,fill=class)) +
    geom_col(width = 0.8) +
    labs(x = "",
         y = "% of mapped reads") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=1,title="",reverse=TRUE)) +
    scale_fill_brewer(palette = "Dark2") +
    coord_flip()
  return(p)
}




## 
## DEhelper.subchunkify: small function stolen from here 
##http://michaeljw.com/blog/post/subchunkify/
## to dynamically create chunks and adjust their size accordingly.
##
DEhelper.subchunkify <- function(g, fig_height=7, fig_width=5) {
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
## DEhelper.Trackhub: display the UCSC trackhub URL
##
DEhelper.Trackhub <- function() {
    
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


