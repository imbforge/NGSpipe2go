#! /usr/bin/Rscript --vanilla --no-save --no-restore


#### libraries ####
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
theme_set(theme_bw(16))
###############

## create plot directory:
dir.create(file.path("figure"), showWarnings = FALSE)


bams <- list.files(pattern='*.bam$')

res <- list()
for (i in 1:length(bams)){
   bam <- bams[i]
   total <- system(
      paste("samtools idxstats", bam, "| awk '{s+=$3+$4} END {print s}'"), intern=TRUE
      )
   mapped <- system(
      paste("samtools idxstats", bam, "| awk '{s+=$3} END {print s}'"), intern=TRUE
      )
   unmmapped <- system(
      paste("samtools view -c -f 4", bam), intern=TRUE
      )
   unique <- system(
      paste("samtools view -c -q 255", bam), intern=TRUE
      )
res[[bam]] <- as.numeric(c(total, mapped, unmmapped, unique))

}

stats <- data.frame(do.call(rbind, res))
stats <- data.frame(
   rownames(stats),
   stats
   )

colnames(stats) <- c("Sample", "Total", "Mapped", "Unmmapped", "Unique")
#reduce size of file names #
stats$Sample <- gsub(".bam|.so|.rg", "", stats$Sample)


## ggplot needs data in a specific layout
MeltedReadCount <- melt(stats[,c(1:3,5)], id=c('Sample'))

names(MeltedReadCount) <- c('Name', 'Reads', 'Counts')

calculatePercentMapped <- function(df){
   df$Percentage <- df$Counts / df$Counts[1] * 100
   df$Percentage[1] <- NA
   return(df)
}

stattable <- ddply(MeltedReadCount, .(Name), calculatePercentMapped)
stattable$Reads <- factor(stattable$Reads,
   levels=c("Total", "Mapped", "Unique"))
stattable$String <- ifelse(is.na(stattable$Percentage ), "",
   paste(round(stattable$Percentage),"%", sep=""))


##########
## Plot ##

outputname <- 'figure/totalReads.pdf'
plotter <- ggplot(stattable, aes(x=Name, y=Counts,fill=Reads)) +
   geom_bar(stat="identity", position="dodge", width=0.7) +
   ggtitle("Mapping Statistics") +
   theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=rel(0.9)), plot.title=element_text(vjust=1, size= rel(0.95)),  axis.title.x=element_text(size=rel(0.9)), axis.title.y=element_text(size=rel(0.9)), axis.text.y=element_text(size=rel(0.9))) +
   xlab("Library ID") +
   ylab("Total Number of Reads") +
   geom_text(size=3, aes(label=String), position=position_dodge(width=0.8), vjust=-0.5) +
   scale_y_continuous(labels=comma) +
   scale_fill_brewer(palette="Set1")

ggsave(outputname, plotter, height=7,width=10)
outputname <- 'figure/totalReads.png'
ggsave(outputname, plotter, height=7,width=10)
