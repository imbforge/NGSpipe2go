usage = "
   Summarizes RNA biotypes by library and small RNA classes. Outputs are plots"

library('data.table')
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
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Set1")
dir.create(file.path("figure"), showWarnings = FALSE)


args <- commandArgs(trailingOnly = TRUE)

# tab="featureCounts_summary.txt"
# bio_file="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/wormbaseGeneID2biotype.txt"
tab <- args[1]
bio_file <- args[2]

gene_features <- read.table(tab, sep="\t", header=T)
nc <- ncol(gene_features)
colnames(gene_features)[1] <- "wormbase_gene"
## remove cols with strand, start...
countdata <- gene_features[,c(1,7:nc)]

## clean RNA classes
colnames(countdata) <- gsub(
   ".+\\.(\\w+.2\\d\\w).bam",
   "\\1",
   colnames(countdata))

## cleans path and adds "full" tag
colnames(countdata) <- gsub(
   ".+\\.(\\w+).bam",
   "\\1.full",
   colnames(countdata)
   )


# head(countdata)

## biotype:
bio <- read.table(bio_file,
   sep="\t")

colnames(bio) <- c("wormbase_gene", "gene_biotype")
# unique(bio$gene_biotype)

idx <- match(countdata$wormbase_gene, bio$wormbase_gene )
countdata$gene_biotype <- bio$gene_biotype[ idx ]


## count non-structural reads
# select full libs
full_libs <- countdata[grepl('gene_biotype|.full', colnames(countdata))]
struct_rna <- c("rRNA", "tRNA", "snoRNA", "snRNA")
non_struct <- subset(full_libs, (!gene_biotype %in% struct_rna), select=-(gene_biotype))
fac <- colSums(non_struct)
non_struct_counts <- cbind(read.table(text = names(fac)), fac)
colnames(non_struct_counts) <- c("exp", "non_struct_counts")

non_struct_counts$exp <- gsub(".full", "", non_struct_counts$exp)

## summarize with data.table
countdata_m <- melt.data.table(data.table(countdata), id.vars=c('wormbase_gene', 'gene_biotype'))
setnames(countdata_m, c('variable', 'value'), c('Bam', 'Counts'))
# countdata_m

## create column with exp to facilitate merging
countdata_m[,base_exp:=gsub("(\\w+).\\w+", "\\1", Bam),]
countdata_m[,small_RNA_class:=gsub("\\w+.(\\w+)", "\\1", Bam),]
countdata_m[,library_treatment:=gsub(".+_(\\w+)", "\\1", base_exp),]

rpms <- merge(countdata_m, non_struct_counts, by.x="base_exp", by.y="exp", all.x=TRUE)
rpms[,RPM:=Counts/non_struct_counts*10^6]

biotype_summary <- rpms[,list(Counts=sum(RPM)), by=list(gene_biotype, Bam, base_exp, small_RNA_class, library_treatment)][!(gene_biotype %in% struct_rna)]


p <- ggplot(biotype_summary, aes(x=base_exp, y=Counts, fill=gene_biotype))

p1 <- p +
   geom_bar(stat="identity") +
   scale_y_continuous(labels = comma) +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(small_RNA_class ~ library_treatment, scales="free_y") +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 10),
      strip.background = element_rect(fill="white")) +
   guides(fill = guide_legend(reverse=TRUE)) +
   ylab("RPM (non-structural reads)") +
   xlab("")

p2 <- p +
   geom_bar(stat="identity", position="dodge") +
   scale_y_continuous(labels = comma) +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(small_RNA_class ~ library_treatment, scales="free_y") +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 10),
      strip.background = element_rect(fill="white")) +
   guides(fill = guide_legend(reverse=TRUE)) +
   ylab("RPM (non-structural reads)") +
   xlab("")

p3 <- p +
   geom_bar(stat="identity", position='fill') +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(small_RNA_class ~ library_treatment, scales="free_y") +
   scale_colour_manual(values = rev(brewer.pal(3,"Set1")[1:2])) +
   theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 10),
      strip.background = element_rect(fill="white")) +
   guides(fill = guide_legend(reverse=TRUE)) +
   ylab("Relative abundance  (a.u.)") +
   xlab("")

pdf("figure/GeneBiotypeSummary.pdf")
p1
p2
p3
dev.off()
