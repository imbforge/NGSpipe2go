#!/usr/bin/env Rscript

##############################################################################
### This script will plot the number of nucleotides around a piRNA read
###############################################################################


libs=c(
   "ggplot2",
   "RColorBrewer",
   "scales"
   )

lapply(libs, require, character.only=T)

# new_theme <- theme_minimal(16) + theme(
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    axis.ticks = element_blank()
#   )
# theme_set(new_theme)

## create plot directory:
dir.create(file.path("figure"), showWarnings = FALSE)

# if (file.exists("pp_freq.txt")){
#    pp_file <- "pp_freq.txt"
# } else {
#    if (file.exists("heterotypical.pp_freq.txt")){
#       pp_file <- "heterotypical.pp_freq.txt"
#    } else {
#       stop("No file with the summary of ping-pong pairs was found")
#    }
# }

testEmptyFile <- function(file_path){
   info = file.info(file_path)
   size = info$size
   if (size == 0){
      TRUE
   } else {
      FALSE
   }
}

plotPP <- function(pp_file){
# '''
# Plots the distribution of 5 prime to 5 prime overlaps in piRNAs.
# Input: a file with the distribution of reads. Col1 is the overlap distance, col2 is the frequency of pairs.
# Ouput: a ggplot
# '''
   pp <- read.table(
   pp_file,
   sep="\t",
   header = FALSE
   )

   colnames(pp) <- c("position", "frequency")
   ## use only pairs with fewer than 30 nts overlap
   pp <- subset(pp, position <= 30)

   ## calculate z-score as in Wasik et al 2015 (Genes Dev)
   pos10 <- pp$frequency[10]
   not10 <- pp$frequency[-10]
   m <- mean(not10)
   std <- sd(not10)
   zscore <- (pos10-m)/std

   zscore_lab <- paste(
      "Z-score =",
      round(zscore, 2)
      )


   ## plot title
   tt <- basename(getwd())
   # tt <- gsub('.pp_freq.txt', '', pp_file)
   pp_plot <- ggplot(pp, aes(x=position, y=frequency)) +
      geom_bar(stat="identity") +
      xlab("\n5'-5' distance") +
      ylab("Read pairs\n") +
      scale_y_continuous(labels = comma, expand = c(0,0)) +
      ggtitle(tt) +
      annotate("text",  x=Inf, y = Inf, label = zscore_lab, vjust=1, hjust=1, size=4.5)+
      theme_bw(16) +
     theme(axis.line = element_line(colour = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank())

   ggsave(paste("figure/", tt, ".ppPlot.pdf", sep=""), pp_plot)
   ggsave(paste("figure/", tt, ".ppPlot.png", sep=""), pp_plot)
}

pp_files <- list.files(pattern='pp_freq.txt')
for (pp_file in pp_files){
   if (testEmptyFile(pp_file)){
      print(paste(pp_file, 'is empty. Skipping'))
   } else {
      print(paste('Plotting', pp_file))
      plotPP(pp_file)
   }
}

