library(edgeR)
tab="featureCounts_summary.txt"
bio_file="/home/adomingu/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/wormbaseGeneID2biotype.txt"

gene_features <- read.table(tab,
   sep="\t",
   header=T,
   row.names=1,
   stringsAsFactors=FALSE
)

nc <- ncol(gene_features)
## remove cols with strand, start...
countdata <- gene_features[,c(6:nc)]

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

## full lib
col_keep <- colnames(countdata)[c(grep(".full", colnames(countdata)))]

#exclude mi_pi_r_t_sn_sno_RNAs
## biotype:
bio <- read.table(bio_file,
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)

colnames(bio) <- c("wormbase_gene", "gene_biotype")
# unique(bio$gene_biotype)
exclude_bio <- c("miRNA", "piRNA", "rRNA", "tRNA", "snoRNA", "snRNA")
row_keep <- subset(bio, !(bio$gene_biotype %in% exclude_bio))$wormbase_gene

sub <- subset(countdata, rownames(countdata) %in% row_keep, select=col_keep)
groups <- gsub("(^\\w+)_\\w+_r\\d+_\\w+.\\w+", "\\1", colnames(sub))
groups <- gsub("(xf)\\d+", "\\1", groups)
groups


## create object for comparisions
cds <- DGEList(sub, group = groups)
names(cds)
list (cds$samples)
dim(cds)
# filter out genes with no or very low coverage (keep RPM > 10 in 3+ samples)
keep <- rowSums( cpm( cds ) > 10) >=3
cds <- cds[keep,]
dim( cds )
cds$samples$lib.size <- colSums(cds$counts)
cds$samples

# perform TMM normalisation
# later use these norm. factors for the 22G and 26G species analyses
cds <- calcNormFactors( cds )
cds$samples
# diagnostic plot
plotMDS ( cds, col=as.numeric(cds$samples$group))
plotMDS ( cds, method="bcv", col=as.numeric(cds$samples$group))
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds )


# DE analysis of the full libraries
et1 <- exactTest( cds, pair=c("N2","xf"))
de1 <- decideTestsDGE(et1, adjust.method="BH", p.value=0.01)
summary(de1)
de1tags <- rownames(cds)[as.logical(de1)]
plotSmear(et1, de.tags=de1tags)
DEtable1 <- topTags(et1, n=Inf)

write.table (DEtable1, "xf-vs-N__edgeR_analysis.txt", sep="\t", quote=F, col.names=NA)
