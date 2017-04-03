#! /usr/bin/Rscript --vanilla --no-save --no-restore


#### libraries ####
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
theme_set(theme_bw(16))
###############


args <- commandArgs(trailingOnly=T)
data.folder   <- args[1]
output.folder <- args[2]
projectprefix <- args[3]

logs.cutadapt <- list.files(path=data.folder,pattern='*.cutadapt.log$',full.names=TRUE)
logs.quality  <- list.files(path=data.folder,pattern='*.fastq_quality_filter.log$',full.names=TRUE) 
logs.dedup    <- list.files(path=data.folder,pattern='*.dedup_stats.txt$',full.names=TRUE)


## get cutadapt trimming data, read in.
res <- list()
for (i in 1:length(logs.cutadapt)){
   log <- logs.cutadapt[i]
   processed.reads <- system(
      paste("grep \"Total reads processed\"", log, "| awk '{print $4}'"), intern=TRUE
      )
   processed.reads <- gsub(",", "", processed.reads)
   trimmed.reads <- system(
      paste("grep \"Reads with adapters\"", log, "| awk '{print $4}'"), intern=TRUE
      )
   trimmed.reads <- gsub(",", "", trimmed.reads)
   tooshort.reads <- system(
      paste("grep \"Reads that were too short\"", log, "| awk '{print $6}'"), intern=TRUE
      )
   tooshort.reads <- gsub(",", "", tooshort.reads)
   ## Trimmed also includes the reads, which are too short and, thus, not kept. correct that.
   trimmed.reads <- as.numeric(trimmed.reads) - as.numeric(tooshort.reads)

res[[log]] <- as.numeric(c(processed.reads,trimmed.reads))

}

stats.cutadapt <- data.frame(do.call(rbind, res))
stats.cutadapt <- data.frame(
   rownames(stats.cutadapt),
   stats.cutadapt
   )

colnames(stats.cutadapt) <- c("sample", "total", "properlyTrimmed")
#reduce size of file names #
stats.cutadapt$sample <- gsub(".cutadapt.log", "", stats.cutadapt$sample)
stats.cutadapt$sample <- gsub(paste(data.folder,"/",sep=""),"",stats.cutadapt$sample)
stats.cutadapt$sample <- gsub(paste0("^",projectprefix),"",stats.cutadapt$sample)


## get quality trimming data, read in.
res <- list()
for (i in 1:length(logs.quality)){
   log <- logs.quality[i]
   kept.reads <- system(
      paste("grep \"Output\"", log, "| awk '{print $2}'"), intern=TRUE
      )
res[[log]] <- as.numeric(c(kept.reads))

}

stats.quality <- data.frame(do.call(rbind, res))
stats.quality <- data.frame(
   rownames(stats.quality),
   stats.quality
   )

colnames(stats.quality) <- c("sample", "keptAfterQFilter")
#reduce size of file names #
stats.quality$sample <- gsub(".fastq_quality_filter.log|.cutadapt", "", stats.quality$sample)
stats.quality$sample <- gsub(paste(data.folder,"/",sep=""),"",stats.quality$sample)
stats.quality$sample <- gsub(paste0("^",projectprefix),"",stats.quality$sample)

## get dedup data, read in.
res <- list()
for (i in 1:length(logs.dedup)){
   log <- logs.dedup[i]
   kept.reads <- system(
      paste("grep \"unique\"", log, "| awk '{print $1}'"), intern=TRUE
      )
res[[log]] <- as.numeric(c(kept.reads))

}

stats.dedup <- data.frame(do.call(rbind, res))
stats.dedup <- data.frame(
   rownames(stats.dedup),
   stats.dedup
   )

colnames(stats.dedup) <- c("sample", "keptAfterDedup")
#reduce size of file names #
stats.dedup$sample <- gsub(".dedup_stats.txt|.cutadapt|.highQ", "", stats.dedup$sample)
stats.dedup$sample <- gsub(paste(data.folder,"/",sep=""),"",stats.dedup$sample)
stats.dedup$sample <- gsub(paste0("^",projectprefix),"",stats.dedup$sample)

## merge all stats
stats.all <- merge(merge(stats.cutadapt,stats.quality,by="sample"),stats.dedup,by="sample")


## ggplot needs data in a specific layout
MeltedReadCount <- melt(stats.all, id=c('sample'))

names(MeltedReadCount) <- c('Name', 'Category', 'Counts')

calculatePercent <- function(df){
   df$Percentage <- df$Counts / df$Counts[1] * 100
   df$Percentage[1] <- NA
   return(df)
}

stattable <- ddply(MeltedReadCount, .(Name), calculatePercent)
stattable$Category <- factor(stattable$Category,
   levels=c("total", "properlyTrimmed", "keptAfterQFilter", "keptAfterDedup"))
stattable$String <- ifelse(is.na(stattable$Percentage ), "",
   paste(round(stattable$Percentage),"%", sep=""))


##########
## Plot ##

plotter <- ggplot(stattable, aes(x=Name, y=Counts,fill=Category)) +
   geom_bar(stat="identity", position="dodge", width=0.7) +
   ggtitle("All trimming statistics") +
   theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=rel(0.9)), plot.title=element_text(vjust=1, size= rel(0.95)),  axis.title.x=element_text(size=rel(0.9)), axis.title.y=element_text(size=rel(0.9)), axis.text.y=element_text(size=rel(0.9))) +
   xlab("Library ID") +
   ylab("Total Number of Reads") +
   geom_text(size=3, aes(label=String), position=position_dodge(width=0.8), vjust=-0.5) +
   scale_y_continuous(labels=comma) +
   scale_fill_manual(values=c("#984EA3","#4DAF4A","#377EB8","#E41A1C"))

ggsave(paste(output.folder,"/allTrimmingStats.pdf",sep=""), plotter, height=7,width=12)
ggsave(paste(output.folder,"/allTrimmingStats.png",sep=""), plotter, height=7,width=12)


