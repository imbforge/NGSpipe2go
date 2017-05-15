usage = "
   Summarizes the number of non-structural reads for each the libraries. It uses the output from featureCounts.
   Inputs
   1: counts table from featureCounts, featureCounts_summary.txt
   2: table with gene names and corresponding biotypes"


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

colnames(countdata) <- gsub(
   ".+\\.(\\w+).bam",
   "\\1",
   colnames(countdata)
   )

## select full libs
countdata <- countdata[!grepl('\\w+.2\\d\\w', colnames(countdata))]

# head(countdata)

## biotype:
bio <- read.table(bio_file,
   sep="\t")

colnames(bio) <- c("wormbase_gene", "gene_biotype")
# unique(bio$gene_biotype)

idx <- match( countdata$wormbase_gene, bio$wormbase_gene )
countdata$gene_biotype <- bio$gene_biotype[ idx ]

struct_rna <- c("rRNA", "tRNA", "snoRNA", "snRNA")
non_struct <- subset(countdata, (!gene_biotype %in% struct_rna))[, -c(1,ncol(countdata))]
head(non_struct)

fac <- colSums(non_struct[ ,c(1:ncol(non_struct))])
# fac2 <- colSums(countdata[ ,c(2:ncol(non_struct))])

write.table(t(data.frame(as.list(fac))),
   "normlization_factors.txt",
   quote=FALSE,
   sep="\t",
   row.names=TRUE,
   col.names=FALSE
   )
