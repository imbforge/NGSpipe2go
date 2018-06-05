# ==========================================================================
# Libraries
# ==========================================================================
library('data.table')
library('ggplot2')
library('dplyr')
library('scales')
library('RColorBrewer')
library('argparser')
library("makeitprettier")

# ==========================================================================
# Arguments
# ==========================================================================
# Create a parser
p <- arg_parser("Plots QC results of PCR duplicate removal. Log files should be in the format ")

# Add command line arguments
p <- add_argument(p,
   "--inputDir",
   help="Input directory",
   type="character"
   )

p <- add_argument(p,
   "--outputDir",
   help="Output directory",
   type="character"
   )

p <- add_argument(p,
   "--projectPrefix",
   help="Name to be used in plots/files",
   type="character",
   default = ""
   )


# Parse the command line arguments
argv <- parse_args(p)
input_dir <- argv$inputDir
output_dir <- argv$outputDir
prefix <- argv$projectPrefix

files <- list.files(
  input_dir,
  '.dedup_stats.log',
  recursive=TRUE,
  full.names=TRUE
  )
out_fig_path <- paste0(output_dir, "/figure")
dir.create(file.path(out_fig_path), showWarnings = FALSE)

stats <- lapply(files, fread)
stats <- rbindlist(stats, use.names=TRUE, fill=FALSE, idcol=NULL)

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
   ggtitle('Putative PCR duplicates') +
   ylab('Number of reads') +
   scale_y_continuous(labels = comma) +
   geom_text(aes(label = Percentage, vjust = -.5), size=4) +
   scale_fill_prettier() +
   coord_flip() +
   theme_poster()

ggsave('figure/PCRDuplicates_test.pdf', stats_p, width=10, height=7)
ggsave('figure/PCRDuplicates_test.png', stats_p, width=10, height=7)
