#! /usr/bin/Rscript --vanilla --no-save --no-restore


#### libraries ####
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(forcats)
theme_set(theme_bw(16))
###############

args <- commandArgs(trailingOnly=T)
data.folder   <- args[1]
output.folder <- args[2]
projectprefix <- args[3]

logs <- list.files(path=data.folder,pattern='*.dedup_stats.txt$',full.names=TRUE)

res <- list()
for (i in 1:length(logs)){
   log <- logs[i]
   processed.reads <- system(
      paste("grep \"highQ\"", log, "| awk '{print $1}'"), intern=TRUE
      )
   kept.reads <- system(
      paste("grep \"unique\"", log, "| awk '{print $1}'"), intern=TRUE
      )
res[[log]] <- as.numeric(c(processed.reads,kept.reads))

}

stats <- data.frame(do.call(rbind, res))
stats <- data.frame(
   rownames(stats),
   stats
   )

colnames(stats) <- c("Sample", "Total", "Kept")
#reduce size of file names #
stats$Sample <- gsub(".dedup_stats.txt|.cutadapt|.highQ", "", stats$Sample)
stats$Sample <- gsub(paste(data.folder,"/",sep=""),"",stats$Sample)
stats$Sample <- gsub(paste0("^",projectprefix),"",stats$Sample)


## df for stacked barplot
stats.stacked <- data.frame(Sample=stats$Sample,Kept=stats$Kept,Removed=stats$Total-stats$Kept)

## ggplot needs data in a specific layout
MeltedReadCount.stacked <- melt(stats.stacked, id=c('Sample'))

names(MeltedReadCount.stacked) <- c('Name', 'Category', 'Counts')


##########
## Plot ##

plotter.stacked <- ggplot(MeltedReadCount.stacked, aes(x=Name, y=Counts,fill=fct_rev(Category))) +
   geom_bar(stat="identity", position="stack", width=0.7) +
   ggtitle("Deduplication Statistics") +
   theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=rel(0.9)), plot.title=element_text(vjust=1, size= rel(0.95)),  axis.title.x=element_text(size=rel(0.9)), axis.title.y=element_text(size=rel(0.9)), axis.text.y=element_text(size=rel(0.9))) +
   labs(x="Library ID", y="Total Number of Reads", fill="Category") +
   scale_y_continuous(labels=comma) +
   scale_fill_manual(values=c("#377EB8","#E41A1C"))

ggsave(paste(output.folder,"/dedupReads.pdf",sep=""), plotter.stacked, height=7,width=10)
ggsave(paste(output.folder,"/dedupReads.png",sep=""), plotter.stacked, height=7,width=10)

