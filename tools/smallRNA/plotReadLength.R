library('data.table')
library('ggplot2')
library('dplyr')
library('scales')
library('RColorBrewer')

theme_set(theme_bw(16))
theme_white <- function() {
     theme_update(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      )
 }
theme_white()

dir.create(file.path("figure"), showWarnings = FALSE)

len_files <- list.files(pattern='readlength.txt')
len_l <- list()

for (len_file in len_files){
   exp <- unlist(strsplit(basename(len_file), '\\.'))[1]
   print(exp)
   counts_tmp <- fread(len_file) %>%
      setnames(c('Count', 'Length'))
   counts_tmp$Sample <- exp
   len_l[[exp]] <- counts_tmp
}

len <- data.frame(do.call(rbind, len_l))

colourCount = length(unique(len$Sample))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


p1 <- ggplot(len, aes(x=factor(Length), y=Count, color=Sample, fill=Sample)) +
   geom_line(aes(group=Sample)) + geom_point() +
   scale_y_continuous(labels = comma) +
   ggtitle('Read length for all mapped reads') +
   xlab('Read length (bp)') +
   ylab('Read count')
ggsave('figure/AllReadsLengthDistribution.pdf', p1, width=15, height=7)
ggsave('figure/AllReadsLengthDistribution.png', p1, width=15, height=7)

p2 <- ggplot(len, aes(x=factor(Length), y=Count, color=Sample, fill=Sample)) +
   geom_line(aes(group=Sample)) + geom_point() +
   facet_grid(Sample ~ ., scales="free_y") +
   scale_y_continuous(labels = comma) +
   ggtitle('Read length for all mapped reads') +
   xlab('Read length (bp)') +
   ylab('Read count')
ggsave('figure/AllReadsLengthDistributionSample.pdf', p2, width=15, height=7)
ggsave('figure/AllReadsLengthDistributionSample.png', p2, width=15, height=7)

# ggplot(len, aes(x=factor(Length), y=Count, color=Sample, fill=Sample)) +
#    geom_line(aes(group=Sample)) + geom_point() +
#    scale_y_continuous(labels = comma) +
#    facet_grid(IP ~ ., scales="free_y") +
#    scale_colour_manual(values=getPalette(colourCount)) +
#    ggtitle('Read length for all mapped reads')
# ggsave('figure/AllReadsLengthDistributionByIP.pdf', width=15, height=7)
# ggsave('figure/AllReadsLengthDistributionByIP.png', width=15, height=7)

# ggplot(len, aes(x=factor(Length), y=Count, color=Genotype, fill=Genotype)) +
#    geom_line(aes(group=Sample)) + geom_point() +
#    scale_y_continuous(labels = comma) +
#    facet_grid(IP ~ ., scales="free_y")
# ggsave('figure/AllReadsLengthDistributionByIP2.pdf', width=15, height=7)
# ggsave('figure/AllReadsLengthDistributionByIP2.png', width=15, height=7)


# ggplot(len, aes(x=factor(Length), y=Count, color=IP, fill=IP)) +
#    geom_line(aes(group=Sample)) + geom_point() +
#    scale_y_continuous(labels = comma) +
#    facet_grid(Genotype ~ ., scales="free_y") +
#    scale_colour_manual(values=getPalette(colourCount)) +
#    ggtitle('Read length for all mapped reads')
# ggsave('figure/AllReadsLengthDistributionByGenotype.pdf', width=15, height=7)
# ggsave('figure/AllReadsLengthDistributionByGenotype.png', width=15, height=7)
