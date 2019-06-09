#! /usr/bin/Rscript --vanilla --no-save --no-restore

# ==========================================================================
# Arguments
# ==========================================================================
library('argparser')

# Create a parser
p <- arg_parser("Plots the proportion small RNA classes per sample.")

# Add command line arguments
p <- add_argument(p,
   "--bams",
   help="Files to count alignments",
   type="character",
   nargs = Inf
   )

p <- add_argument(
   p,
   "--libsizes",
   help = "File with library sizes, number of reads to normalize the read counts. First column should be the file name and the second the number of reads.",
   type = "character",
   default = NA
)

p <- add_argument(
   p,
   "--title",
   help="Plot title.",
   type="character",
   default = ""
   )

p <- add_argument(
   p,
   "--prefix",
   help="Prefix to add to every file and plot saved.",
   type="character",
   default = "small_RNA_classes"
   )

p <- add_argument(
   p,
   "--out",
   help = "Where should the files be saved. Directory is created if it doesn't exist.",
   type = "character",
   default = "./"
)

# Parse the command line arguments
argv <- parse_args(p)
bams <- argv$bams
fac_path <- argv$libsizes
title <- argv$title
prefix <- argv$prefix
out <- argv$out

prefix <- gsub("\\s+", "_", prefix)
prefix <- gsub(
  "[[:punct:]]",
  "_",
  prefix,
  )


# ==========================================================================
# settings
# ==========================================================================
#### libraries ####
library("data.table")
library("ggplot2")
library("makeitprettier")
library('lemon')
dir.create(
   paste0(out, "/", "figure"),
   showWarnings = FALSE
)

# ==========================================================================
# Counts
# ==========================================================================

res <- list()
for (i in 1:length(bams)){
   bam <- bams[i]
   print(bam)
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

colnames(stats) <- c("Bam", "Total", "Mapped", "Unmmapped", "Unique")
stats <- stats[, c(1:2)]
#reduce size of file names #
stats$Bam <- gsub(".unique|.bam|.so|.rg", "", basename(as.character(stats$Bam)))
stats$Sample <- gsub("\\..*", "", basename(stats$Bam))
stats$Class <- gsub("\\w+\\.(\\w+)", "\\1", stats$Bam)

message("Reading normalization factors")
fac <- fread(fac_path)
setnames(fac, c("Sample", "n_reads"))
fac[, Sample:=gsub(".All|.bam", "", basename(Sample))]
setDT(stats)
all_stats <- merge(stats, fac, by.x="Sample", by.y="Sample") 
all_stats[,Bam:=NULL,]
all_stats[,Perc:=round(Total / n_reads * 100, 1)]
all_stats[,RPM:=Total / n_reads * 10^6]
setnames(
  all_stats,
  c("Total", "n_reads"),
  c("counts", "libSize"),
  )


# ==========================================================================
# bar plot
# ==========================================================================
n_classes <- length(unique(all_stats$Class))

p1 <- ggplot(all_stats, aes(x = forcats::fct_rev(Sample), y = RPM, fill = Class)) +
   geom_bar(stat = "identity") +
   lemon::facet_rep_wrap(
      Class ~ .,
      scales = "free_x",
      ncol = 1,
      repeat.tick.labels = 'bottom'
 ) +
   coord_flip() +
 labs(
   x = "", 
   y = "RPM"
 ) +
   expand_limits(y = 0) +
   scale_y_continuous(labels = comma) +
   scale_fill_prettier() +
   theme_redl(base_size = 18) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(
   paste0(out, "/figure/", prefix, ".smallRNAClassesNormalized"),
   p1,
   width = 8,
   height = n_classes * 3.5 + 1
   )



p2 <- ggplot(all_stats, aes(x = forcats::fct_rev(Sample), y = counts, fill = Class)) +
   geom_bar(stat = "identity") +
   lemon::facet_rep_wrap(
      Class ~ .,
      scales = "free_x",
      ncol = 1,
      repeat.tick.labels = 'bottom'
 ) +
   coord_flip() +
 labs(
   x = "", 
   y = "Total small RNAs in library"
 ) +
   expand_limits(y = 0) +
   scale_y_continuous(labels = comma) +
   scale_fill_prettier() +
   theme_redl(base_size = 18) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(
   paste0(out, "/figure/", prefix, ".smallRNAClassesRaw"),
   p2,
   width = 8,
   height = n_classes * 3.5 + 1
   )


# ==========================================================================
# Save table
# ==========================================================================

write.table(
   all_stats,
   file = paste0(out, "/smallRNAClassesCount.txt"),
   quote = FALSE,
   sep = "\t", 
   row.names = FALSE
)