#! /usr/bin/Rscript --vanilla --no-save --no-restore

# ==========================================================================
# Arguments
# ==========================================================================
library('argparser')

# Create a parser
p <- arg_parser("Plots the proportion small RNA classes per sample.")

# Add command line arguments
p <- add_argument(p,
   "--inputs",
   help="Path to files with summary counts",
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
   default = "all_reads"
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
inputs <- argv$inputs
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
library("lemon")
dir.create(
   paste0(out, "/", "figure"),
   showWarnings = FALSE
)

# ==========================================================================
# Counts
# ==========================================================================
# inputs <- list.files(pattern='*.nuc_bias.txt')
# cov_files <- cov_files[1:10]

read_and_summarize <- function(input_file){
	print(input_file)
	cov_tmp <- fread(input_file)
	setnames(cov_tmp, c('Length', 'Nucleotide', 'Count'))
	cov_tmp[, Sample := gsub('.nuc_bias.txt', '', input_file), ]
	return(cov_tmp)
}

covs <- lapply(inputs, read_and_summarize)
counts <- rbindlist(covs, fill=TRUE)
counts[, Nucleotide := gsub("T", "U", Nucleotide)]
counts[, Length := as.integer(Length)]

fac <- fread(fac_path)
setnames(fac, c("Sample", "n_reads"))
fac[, Sample:=gsub(".All|.bam", "", basename(Sample))]

counts <- merge(counts, fac)
counts[, RPM:=Count/n_reads * 10^6, ]

n_samples <- length(unique(counts$Sample))
n_pos <- length(unique(counts$Length))

ggplot(counts, aes(x=Length, y=Count, fill=Nucleotide)) +
   geom_bar(stat='identity') +
   scale_fill_prettier() +
   scale_y_continuous(labels=comma) +
     lemon::facet_rep_wrap(
      Sample ~ .,
      scales = "free_x",
      ncol = 1,
      repeat.tick.labels = 'bottom'
 ) +
   expand_limits(y = 0) +
   theme_redl(base_size = 18) +
   labs(
      x="Read length",
      y="Number of reads"
   )
save_plot(
   "figure/nucleotide_bias_read_length",
   save_data = TRUE,
   width = n_pos * 0.75 + 1,
   height = n_samples * 3
)

ggplot(counts, aes(x=Length, y=RPM, fill=Nucleotide)) +
   geom_bar(stat='identity') +
   scale_fill_prettier() +
   scale_y_continuous(labels=comma) +
     lemon::facet_rep_wrap(
      Sample ~ .,
      scales = "free_x",
      ncol = 1,
      repeat.tick.labels = 'bottom'
 ) +
   expand_limits(y = 0) +
   theme_redl(base_size = 18) +
   labs(
      x="Read length",
      y="RPM"
   )
save_plot(
   "figure/nucleotide_bias_read_length.normalized",
   save_data = TRUE,
   width = n_pos * 0.75 + 1,
   height = n_samples * 3
)

#
# highlight 26G
# --------------------------------------------------------------------------
small_counts <- counts[ Length >= 25 & Length <= 27 ]
n_samples <- length(unique(small_counts$Sample))
n_pos <- length(unique(small_counts$Length))

ggplot(small_counts, aes(x=Length, y=Count, fill=Nucleotide)) +
   geom_bar(stat='identity') +
   scale_fill_prettier() +
   scale_y_continuous(labels=comma) +
     lemon::facet_rep_wrap(
      Sample ~ .,
      scales = "free_x",
      ncol = 1,
      repeat.tick.labels = 'bottom'
 ) +
   expand_limits(y = 0) +
   theme_redl(base_size = 18) +
   labs(
      x="Read length",
      y="Number of reads"
   )
save_plot(
   "figure/nucleotide_bias_read_length.25_27nt",
   width = n_pos * 1.5 + 1,
   height = n_samples * 3
)

ggplot(small_counts, aes(x=Length, y=RPM, fill=Nucleotide)) +
   geom_bar(stat='identity') +
   scale_fill_prettier() +
   scale_y_continuous(labels=comma) +
     lemon::facet_rep_wrap(
      Sample ~ .,
      scales = "free_x",
      ncol = 1,
      repeat.tick.labels = 'bottom'
 ) +
   expand_limits(y = 0) +
   theme_redl(base_size = 18) +
   labs(
      x="Read length",
      y="RPM"
   )
save_plot(
   "figure/nucleotide_bias_read_length.25_27nt.normalized",
   width = n_pos * 1.5 + 1,
   height = n_samples * 3
)
