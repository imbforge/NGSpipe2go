# ==========================================================================
# Arguments
# ==========================================================================
library("argparser")
# Create a parser
p <- arg_parser("Plots summary of piRNAs - length and other biases.")

# Add command line arguments
p <- add_argument(
  p,
   "--inputs",
   help = "Path input counts",
   type = "character",
   nargs = Inf
)

p <- add_argument(
   p,
   "--out",
   help = "Where should the files be saved. Directory is created if it doesn't exist.",
   type = "character",
   default = "./"
)

p <- add_argument(
   p,
   "--libsizes",
   help = "File with library sizes, number of reads to normalize the read counts. First column should be the file name and the second the number of reads.",
   type = "character",
   default = NA
)

# Parse the command line arguments
argv <- parse_args(p)
inputs <- argv$inputs
libsizes <- argv$libsizes
out <- argv$out


# ==========================================================================
# Libraries
# ==========================================================================
library("data.table")
library("magrittr")
library("makeitprettier")
library("ggplot2")
library("scales")
library("RColorBrewer")
library("lemon")

out_fig_path <- paste0(out, "/figure")
dir.create(file.path(out_fig_path), showWarnings = FALSE)

data_file <- paste0(out, "/", "piRNA_quantification.RData")

# inputs <- list.files(pattern=".counts")[1:4]
## check if file is empty
info <- file.info(inputs)
notempty <- rownames(info[!(info$size == 0), ])

counts_l <- list()
for (counts_file in notempty){
   exp <- gsub(".counts", "", basename(counts_file))
   print(exp)
   # counts_tmp$exp <- exp
   counts_l[[exp]] <- fread(counts_file, drop = c(1, 4:6)) %>%
      setnames(c("start_read", "end_read", "chr_re", "start_re", "end_re", "repFamily", "repName", "strand_re", "repClass")) %>%
      .[, readlength := end_read - start_read] %>%
      .[, c("end_read","start_read") := NULL] %>%
      .[, exp := exp] %>%
      splitstackshape:::cSplit("exp", ".", drop = FALSE, direction = "wide") %>%
      setnames(c("exp_1", "exp_2"), c("Sample", "Mapping"))
}

counts <- data.table::rbindlist(counts_l, fill = TRUE)
counts[, repClass:=gsub("\\?", "", repClass)]
# set factor levels
counts[, Mapping:=relevel(Mapping, "sense")]

## group te elements, DNA and RNA
retro_te_classes <- c("LTR", "SINE", "LINE")
counts[, Feature:=ifelse(repClass %in% retro_te_classes, "RNA_repeats", repClass),]
counts[, Feature:=ifelse(Feature %in% c("DNA", "DNA?"), "DNA_repeats", Feature),]
counts[, Feature:=ifelse(grepl("protein_coding", Feature), "protein_coding", Feature),]
counts[, Feature:=ifelse(grepl("nonsense_mediated_decay", Feature), "protein_coding", Feature),]
counts[, Feature:=ifelse(grepl("processed_transcript", Feature), "protein_coding", Feature),]


len <- counts[, list(.N), by = list(Sample, readlength, Mapping)][order(Sample, readlength, Mapping)]
n_samples <- length(unique(len$Sample))

len_ip_p <- ggplot(len, aes(x = factor(readlength), y = N, color = Mapping)) +
   geom_line(aes(group = Mapping)) + geom_point() +
   scale_y_continuous(labels = comma) +
   facet_rep_wrap(Sample ~ ., ncol = 1, repeat.tick.labels = TRUE) +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   xlab("read length") +
   ylab("Number of reads") +
   theme_redl(base_size = 18)

save_plot(
   paste0(out_fig_path, "/AllLibLengthDistributionBySample"),
   len_ip_p,
   width = 10,
   height= 5 * n_samples + 0.5
)

classes <- c("DNA_repeats", "RNA_repeats")
re_only <- counts[Feature %in% classes]

re_only_length <- re_only[, list(.N), by = list(Sample, Mapping, readlength)]

length_ip_class_p <- ggplot(re_only_length, aes(x = factor(readlength), y = N, color = Mapping)) +
   geom_line(aes(group = Mapping)) + geom_point() +
   scale_y_continuous(labels = comma) +
   facet_rep_wrap(Sample ~ ., ncol = 1, repeat.tick.labels = TRUE) +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   xlab("read length") +
   ylab("Number of reads") +
   labs(
         caption = paste(classes, collapse = ", ")
   ) +
   theme_redl(base_size = 18)


save_plot(
   paste0(out_fig_path, "/AllLibLengthDistributionBySample.TE_only"),
   len_ip_p,
   width = 10,
   height= 5 * n_samples + 0.5
)

# ==========================================================================
#  rpm normalization
# ==========================================================================
fac <- fread(libsizes) %>%
  setnames(c("Exp", "libSize")) %>%
  .[, Exp := gsub("\\..+", "", basename(Exp)), ] %>%
  .[]

len <- merge(len, fac, by.x="Sample", by.y="Exp")
len[, RPM:=N/(libSize) * 10^6, ]

len_ip_p <- ggplot(len, aes(x = factor(readlength), y = RPM, color = Mapping)) +
   geom_line(aes(group = Mapping)) + geom_point() +
   scale_y_continuous(labels = comma) +
   facet_rep_wrap(Sample ~ ., ncol = 1, repeat.tick.labels = TRUE) +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   xlab("read length") +
   ylab("Number of reads") +
   theme_redl(base_size = 18)


save_plot(
   paste0(out_fig_path, "/AllLibLengthDistributionBySample.RPM"),
   len_ip_p,
   width = 10,
   height= 5 * n_samples + 0.5
)


#
# REpeat elements only  - RNA and DNA
# --------------------------------------------------------------------------
re_only_length <- re_only[, list(.N), by = list(Sample, Mapping, readlength)]
re_only_length <- merge(re_only_length, fac, by.x="Sample", by.y="Exp")
re_only_length[, RPM:=N/(libSize) * 10^6, ]

length_ip_class_p <- ggplot(re_only_length, aes(x = factor(readlength), y = N, color = Mapping)) +
   geom_line(aes(group = Mapping)) + geom_point() +
   scale_y_continuous(labels = comma) +
   facet_rep_wrap(Sample ~ ., ncol = 1, repeat.tick.labels = TRUE) +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   xlab("read length") +
   ylab("Number of reads") +
   labs(
      caption = paste(classes, collapse = ", ")
   ) +
   theme_redl(base_size = 18)

save_plot(
   paste0(out_fig_path, "/AllLibLengthDistributionBySample.RPM.TE_only"),
   len_ip_p,
   width = 10,
   height= 5 * n_samples + 0.5
)


# ==========================================================================
# Save object
# ==========================================================================
save(counts, file = data_file)
