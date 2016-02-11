
library("data.table")
library("dplyr")
library("ggplot2")
library("scales")
library("RColorBrewer")

theme_set(theme_bw(16))
theme_white <- function() {
     theme_update(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      )
 }
theme_white()
# make a palette the default: http://stackoverflow.com/a/16437625/1274242
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Set1")
dir.create(file.path("figure"), showWarnings = FALSE)


files <- list.files(pattern='.counts')
## check if file is empty
info <- file.info(files)
notempty <- rownames(info[!(info$size == 0), ])

counts_l <- list()
for (counts_file in notempty){
   exp <- gsub(".counts", "", counts_file)
   exp <- gsub('so.rg.', '', exp)
   exp <- gsub('kh-', 'Tdrkh-', exp)
   split_exp <- unlist(strsplit(exp, '\\.'))
   exp <- paste(split_exp[1], split_exp[length(split_exp)], sep='.')
   print(exp)
   counts_tmp <- fread(counts_file) %>%
   setnames(c('chr_read', 'start_read', 'end_read', 'name_read', 'score_read', 'strand_read', 'chr_re', 'start_re', 'end_re', 'repFamily', 'repName', 'strand_re', 'repClass'))
   counts_tmp$exp <- exp
   counts_l[[exp]] <- counts_tmp %>%
      mutate(readlength = end_read - start_read) %>%
      filter(readlength < 40) %>%
      splitstackshape:::cSplit('exp', '.', drop=FALSE, direction='wide') %>%
      setnames(c('exp_1', 'exp_2'), c('Sample', 'Mapping')) %>%
      splitstackshape:::cSplit('Sample', '-', drop=FALSE, direction='wide') %>%
      setnames(c('Sample_1', 'Sample_2', 'Sample_3', 'Sample_4'), c('Genotype1', 'Genotype2', 'Tissue', 'IP')) %>%
      mutate(Genotype=paste(Genotype1, Genotype2, sep='-'))
}

counts <- data.table::rbindlist(counts_l, fill=TRUE)
counts[, repClass:=gsub('\\?', '', repClass)]
# set factor levels
counts[, Mapping:=relevel(Mapping, "sense")]

## group te elements, DNA and RNA
retro_te_classes <- c('LTR', 'SINE', 'LINE')
counts[, Feature:=ifelse(repClass %in% retro_te_classes, "RNA_repeats", repClass),]
counts[, Feature:=ifelse(Feature %in% c("DNA", "DNA?"), "DNA_repeats", Feature),]
counts[, Feature:=ifelse(grepl("protein_coding", Feature), "protein_coding", Feature),]
counts[, Feature:=ifelse(grepl("nonsense_mediated_decay", Feature), "protein_coding", Feature),]
counts[, Feature:=ifelse(grepl("processed_transcript", Feature), "protein_coding", Feature),]

len <- counts[, list(.N), by=list(Sample, Genotype, IP, Tissue, readlength, Mapping)][order(Sample, Genotype, IP, Tissue, readlength, Mapping)]
n_samples <- length(unique(len$Sample))

# colourCount = length(unique(len$Sample))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))

param = "IP ~ Tissue"
len_ip_p <- ggplot(len, aes(x=factor(readlength), y=N, color=Mapping)) +
   geom_line(aes(group=Mapping)) + geom_point() +
   scale_y_continuous(labels = comma) +
   facet_grid(Sample ~ ., scales='free') +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   theme(strip.text = element_text(size = 14),
      # axis.title.y = element_text(size = 10),
      # axis.title.x = element_text(size = 10),
      strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 8),
      strip.background = element_rect(fill="white")) +
   xlab('read length') +
   ylab('Number of reads')

ggsave('figure/AllLibLengthDistributionBySample.pdf', len_ip_p, width=15, height= 2* n_samples)
ggsave('figure/AllLibLengthDistributionBySample.png', len_ip_p, width=15, height= 2* n_samples)



classes <- c("DNA_repeats", "RNA_repeats", "protein_coding")
re_only <- counts[Feature %in% classes]

re_only_length <- re_only[, list(.N), by=list(Feature, Sample, Mapping, readlength)]

length_ip_class_p <- ggplot(re_only_length, aes(x=factor(readlength), y=N, color=Mapping, fill=Mapping)) +
   geom_line(aes(group=Mapping)) +
   geom_point() +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(Sample ~ Feature, scales="free_y") +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   scale_y_continuous(labels = comma, name='Total read counts', breaks=pretty_breaks(n=2)) +
    theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 8),
      strip.background = element_rect(fill="white"))

ggsave('figure/AllLibLengthDistributionBySampleFeature.pdf',
   length_ip_class_p,
   width=20,
   height=2* n_samples)
ggsave('figure/AllLibLengthDistributionBySampleFeature.png',
   length_ip_class_p,
   width=20,
   height=2* n_samples)


## Mapping to Repeat elements
counts_sum <- counts[, list(.N), by=list(Feature, Sample)][order(Feature, Sample)]

colourCount = length(unique(counts_sum$Feature))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

features_p <- ggplot(counts_sum, aes(x=Sample, y=N, fill=Feature)) +
   geom_bar(stat="identity", position='fill') +
   theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
   # facet_grid(Tissue ~ IP) +
   scale_fill_manual(values = getPalette(colourCount)) +
   scale_y_continuous(labels = percent_format(), name='% of reads mapped to Feature')
ggsave('figure/PercentageOfFeature.pdf', features_p)
ggsave('figure/PercentageOfFeature.png', features_p)


## Mapping to RNA elements
rna <- counts[Feature %in% "RNA_repeats"]

rna_sum <- counts[Feature %in% "RNA_repeats"][, list(.N), by=list(repClass, Sample)][order(repClass, Sample)]

colourCount = length(unique(rna_sum$repClass))
features_p2 <- ggplot(rna_sum, aes(x=Sample, y=N, fill=repClass)) +
   geom_bar(stat="identity", position='fill') +
   theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
   # facet_grid(Tissue ~ IP) +
   scale_fill_manual(values = getPalette(colourCount)) +
   scale_y_continuous(labels = percent_format(), name='% of reads mapped to RNA repeats')

ggsave('figure/PercentageOfFeatureRepeats.pdf', features_p2)
ggsave('figure/PercentageOfFeatureRepeats.png', features_p2)
