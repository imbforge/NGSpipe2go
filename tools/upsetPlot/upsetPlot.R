#############################################################################
##
## What: upsetPlot.R
## Who: Frank Rühle
## When: 2021-07-08
##
## prepare combination matrix (time-consuming) and UpSet Plot for peak data
##
## Args: 
## -----
## peakData=                       # path to the xls file result from MACS2
############################################################################
library("ComplexHeatmap")	
library("Cairo")
library("grid")
library("RColorBrewer")


##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  
  if(!is.null(default)) default else do.call(convert,list(NA))
}
run_custom_code <- function(x) {
  eval(parse(text=x))
}


args             <- commandArgs(T)
PEAKDATA         <- parseArgs(args, "peakData=") # subdir already in path
FTARGETS         <- parseArgs(args, "targets=","targets.txt")
OUT              <- parseArgs(args, "out=", "upsetPlot") # subdir already in path
MODE             <- parseArgs(args, "mode=", "intersect") 
PEAKOVERLAPMODE  <- parseArgs(args, "peakOverlapMode=", "peaknumber") 
SETSIZE          <- parseArgs(args, "setsize=", 25, "as.numeric") 
ADDBARANNOTATION <- parseArgs(args, "addBarAnnotation=", FALSE, "as.logical")  

runstr <- paste0("Call with: Rscript upsetPlot.R [peakData=",PEAKDATA,"] [targets=",FTARGETS,"] [out=",OUT,"] [mode=",MODE,"] [peakOverlapMode=",PEAKOVERLAPMODE,"] [setsize=",SETSIZE,"] [addBarAnnotation=",ADDBARANNOTATION,"]")
cat(runstr)


#' Define a group color palette.
#' @param num - number of groups
#' @return a color vector
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

#' Define a replicate color palette.
#' @param num - number of replicates
#' @return a color vector
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
    if(!file.exists(file.path(PEAKDATA,peaksSubdir))) {
      return("MACS2 results not available")
    }
    if(file.exists(file.path(PEAKDATA,peaksSubdir, paste0(targets$IPname, ".vs.", targets$INPUTname,"_macs2_blacklist_filtered_peaks.xls")))[1]) {
      comparisons <- paste0(targets$IPname, ".vs.", targets$INPUTname, "_macs2_blacklist_filtered_peaks.xls")  
    } else {
      comparisons <- paste0(targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls")
    }
    exist <- sapply(paste0(file.path(PEAKDATA,peaksSubdir), "/", comparisons), file.exists)
    targets <- targets[exist, ]
    
    columnNames2replace <- c(seqnames="chr", abs_summit="summit", pileup="tags", X.log10.pvalue.="-log10 pval", fold_enrichment="fold enrichment", X.log10.qvalue.="-log10 FDR")
    
    # and return the tables
    if(file.exists(paste0(file.path(PEAKDATA,peaksSubdir), "/", targets$IPname, ".vs.", targets$INPUTname,"_macs2_blacklist_filtered_peaks.xls"))[1]) {
      peaks <- lapply(paste0(file.path(PEAKDATA,peaksSubdir), "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_blacklist_filtered_peaks.xls"), function(x) {
        if(as=="data.frame") {
          x <- tryCatch(read.delim(x, comment.char="#", stringsAsFactors=F), error=function(e) as.data.frame(matrix(ncol=10)))
          colnames(x) <- dplyr::recode(colnames(x), !!!columnNames2replace)
          x[order(x$chr, x$start, x$end),  !colnames(x) %in% c("-log10 pval", "name")]
        } else {
          x <- tryCatch(ChIPseeker::readPeakFile(x))
        }
      })
    } else {
      peaks <- lapply(paste0(file.path(PEAKDATA,peaksSubdir), "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls"), function(x) {
        if(as=="data.frame") {
          x <- tryCatch(read.delim(x, comment.char="#", stringsAsFactors=F), error=function(e) as.data.frame(matrix(ncol=10)))
          #colnames(x) <- plyr::revalue(colnames(x), replace=columnNames2replace , warn_missing = F)
          colnames(x) <- dplyr::recode(colnames(x), !!!columnNames2replace)
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



#'
#' ChIPhelper.UpsetPlot: shows a upset plot for the peaks called in one branch of the pipeline
#'
#' @param subdir character with sub-directory containing the peak files
#' @param Mode character defining the mode for forming the combination set (one of "distinct", "intersect", "union")
#' @param peakOverlapMode select the value function to calculate size of combination sets ("peaknumber" for number of peaks and/or "bp" for basepairs)
#' @param setsize numeric, maximal number of sets shown
#' @param targetsdf data.frame with targets data. If not NULL, combination sets are highlighted by exclusive sample groups.
#' @param addBarAnnotation logical, whether the intersection sizes are printed on top pf the column annotation 

ChIPhelper.UpSetPlot <- function(subdir="", Mode = "intersect", peakOverlapMode=c("peaknumber"), setsize=25, targetsdf=targets, addBarAnnotation=F){
  
  #create granges from the peaks
  peak.ranges <- ChIPhelper.init("readPeaks", subdir, peaks_as="GRanges")
  if(length(peak.ranges)>15) {return("UpSet plots are available only for max 15 peaksets")}
  upsetReturn <- list()
  upsetReturn[["peak.ranges"]] <- peak.ranges
  
  if("peaknumber" %in% peakOverlapMode) {
    cat(paste0("#### Overlap of peaks per peak number"), fill=T)
    cat("\n", fill=T)
    #this upset plot is based on the number of peaks which are overlapping (value_fun is length)
    upset_matrix <- make_comb_mat(peak.ranges, mode = Mode, value_fun = length)
    #we subset the matrix to only display the top sets
    upset_matrix <- upset_matrix[order(comb_size(upset_matrix), decreasing = T)[1:setsize]]
    upsetReturn[["matrix_peaknumber"]] <- upset_matrix

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
    
    # plot upsetplot_peaknumber
    png(file.path(OUT, paste0("upsetplot_", Mode, "_peaknumber.png")), width = 200, height = 100, units = "mm", res=300)
    
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
    dev.off()

    cat("\n", fill=T)
  }
  
  if("bp" %in% peakOverlapMode) {
    cat("\n", fill=T)
    cat(paste0("#### Overlap of peaks based on bp"), fill=T)
    cat("\n", fill=T)
    cat("\n", fill=T)
    cat("\n", fill=T)
    #this upset plot is based on the number of bp which are overlapping 
    upset_matrix <- make_comb_mat(peak.ranges, mode = Mode)
    #we subset the matrix to only display the top sets
    upset_matrix <- upset_matrix[order(comb_size(upset_matrix), decreasing = T)[1:setsize]]
    upsetReturn[["matrix_bp"]] <- upset_matrix
    
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
    
    # plot upsetplot_bp
    png(file.path(OUT, paste0("upsetplot_", Mode, "_bp.png")), width = 200, height = 100, units = "mm", res=300)
    
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
    dev.off()
    cat("\n", fill=T)
    cat("\n", fill=T)
  }
  
  return(upsetReturn)
}

# load targets file
if(file.exists(FTARGETS)){
  targets <- read.delim(FTARGETS, stringsAsFactors = F)
} else {targets <- NULL}

# run the helper function
upsetPlotList <- ChIPhelper.UpSetPlot(Mode=MODE, peakOverlapMode=PEAKOVERLAPMODE, setsize=SETSIZE, targetsdf=targets, addBarAnnotation=ADDBARANNOTATION)


# save the session Information
writeLines(capture.output(sessionInfo()), paste(OUT, "/upsetPlot_session_info.txt", sep=""))
  
save(upsetPlotList, file=paste0(OUT,"/upsetPlot.RData"))










