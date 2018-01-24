library('ggplot2')
library("scales")
theme_set(theme_bw(16))
theme_white <- function() {
     theme_update(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      )
 }
theme_white()
# make a palette the default: http://stackoverflow.com/a/16437625/1274242
ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_brewer(palette="Set1")

dir.create(file.path("figure"), showWarnings = FALSE)



stats <- read.delim('dedup.stats.txt', sep=" ", header=FALSE)
stats <- subset(stats, !grepl('\\*', stats$V2))
stats$Condition <- factor(gsub('(\\w+)\\..+', '\\1', stats$V2))
stats$Reads <- factor(ifelse(grepl('unique', stats$V2), 'Unique', 'Original'), levels = c("Original", "Unique"))
colnames(stats)[1:2] <- c('Nreads', 'File')

stats_wide <- tidyr::spread(stats[,c(1,3:4)], Reads, Nreads)
stats_wide$Percentage <- paste(
   round(stats_wide$Unique / stats_wide$Original * 100,1),
   '%',
   sep=''
)

idx <- match( stats$Condition, stats_wide$Condition )
stats$Percentage <- stats_wide$Percentage[ idx ]
stats$Percentage <- ifelse(stats$Reads == 'Original', '', stats$Percentage)

stats_p <- ggplot(stats, aes(y=Nreads, x=Condition, fill=Reads)) +
   geom_bar(stat='identity', position=position_dodge()) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   ggtitle('Putative PCR duplicates') +
   ylab('Number of reads') +
   scale_y_continuous(labels = comma) +
   geom_text(aes(label = Percentage, vjust = -.5), size=4)

ggsave('figure/PCRDuplicates.pdf', stats_p, width=10, height=7)
ggsave('figure/PCRDuplicates.png', stats_p, width=10, height=7)
