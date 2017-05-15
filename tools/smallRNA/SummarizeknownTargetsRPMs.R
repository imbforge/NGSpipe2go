usage = "
   Summarizes RNA biotypes by library and small RNA classes. Outputs are plots"

library('data.table')
library('ggplot2')
library('RColorBrewer')
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
# target_file <- '~/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans/data/annotations_targets/22G/all_targets.txt'
tab <- args[1]
bio_file <- args[2]
target_file <- args[3]

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
   "\\1.All",
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
full_libs <- countdata[grepl('gene_biotype|.All', colnames(countdata))]
struct_rna <- c("rRNA", "tRNA", "snoRNA", "snRNA")
non_struct <- subset(full_libs, (!gene_biotype %in% struct_rna), select=-(gene_biotype))
fac <- colSums(non_struct)
non_struct_counts <- cbind(read.table(text = names(fac)), fac)
colnames(non_struct_counts) <- c("exp", "non_struct_counts")

non_struct_counts$exp <- gsub(".All", "", non_struct_counts$exp)

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

## get rid of columns with unnecessary information
rpms[,gene_biotype:=NULL]
rpms[,Counts:=NULL]
rpms[,non_struct_counts:=NULL]

## genotype
rpms[,genotype:=gsub("(^\\w+)_\\w+_r\\d+_\\w+.\\w+", "\\1", base_exp),]
rpms[,genotype:=gsub("(xf)\\d+", "\\1", genotype),]
## replicates
rpms[,replicate:=gsub("^\\w+_\\w+_(r\\d+)_\\w+.\\w+", "\\1", base_exp),]
rpms
rpms_c <- dcast(rpms, small_RNA_class + library_treatment + wormbase_gene ~ genotype +replicate, value.var = "RPM")

## ratio
rpms_c[ , avg_N2:=rowMeans(.SD), .SDcols = grep("N2", colnames(rpms_c))]
rpms_c[ , avg_xf:=rowMeans(.SD), .SDcols = grep("xf", colnames(rpms_c))]
rpms_c[ , RPM_Ratio:=avg_xf / (avg_xf + avg_N2),]

# clean
rpms_c[, grep("N2|xf", names(rpms_c)) := NULL]

## Target classes:
targets <- read.delim(target_file, header = TRUE, stringsAsFactors=FALSE)
# test <- dcast(targets, geneid ~ exp)
# head(test)
targets$exp <- (gsub(" ", "_", targets$exp))

library(biomaRt)
wormbase=useMart("ensembl", dataset="celegans_gene_ensembl")

ensembl = useMart("ensembl", dataset = "celegans_gene_ensembl")
genemap <- getBM( attributes = c("wormbase_gene", "wormbase_gene_seq_name"),
                  filters = "wormbase_gene_seq_name",
                  values = targets$geneid,
                  mart = ensembl)

idx <- match( targets$geneid, genemap$wormbase_gene_seq_name )
targets$wormbase_gene <- genemap$wormbase_gene[ idx ]

as.data.frame(table(targets$exp))

## cast into wide, where if a gene is a target it will be marked as "1"
## if not target it will be "0"
# tgt <- dcast(targets, geneid + wormbase_gene ~ exp, value.var="geneid")

## merge
# add geneid to rpms
genemap <- getBM( attributes = c("wormbase_gene", "wormbase_gene_seq_name"),
                  filters = "wormbase_gene",
                  values = rpms_c$wormbase_gene,
                  mart = ensembl)

idx <- match(rpms_c$wormbase_gene, genemap$wormbase_gene)
rpms_c$geneid <- genemap$wormbase_gene_seq_name[ idx ]


rpms_targets <- data.table(merge(targets, rpms_c, by="geneid"))

ggplot(rpms_targets, aes(x=exp, y=RPM_Ratio, fill=exp)) +
   geom_boxplot(notch=TRUE) +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(small_RNA_class ~ library_treatment, scales="free_y") +
   scale_fill_manual(values = rev(brewer.pal(7,"Set1"))) +
   theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 10),
      strip.background = element_rect(fill="white")) +
   geom_hline(yintercept=0.5, linetype="dashed", color="gray") +
   xlab("") +
   ylab("Normalized reads in mutant relative to wild-type\n(xf/(xf+N2))\n")

ggsave("figure/known_targets_rpms.boxplot.pdf")
ggsave("figure/known_targets_rpms.boxplot.png")

ggplot(rpms_targets, aes(x=exp, y=RPM_Ratio, fill=exp)) +
   geom_violin() +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(small_RNA_class ~ library_treatment, scales="free_y") +
   scale_fill_manual(values = rev(brewer.pal(7,"Set1"))) +
   theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 10),
      strip.background = element_rect(fill="white")) +
   geom_hline(yintercept=0.5, linetype="dashed", color="gray") +
   xlab("") +
   ylab("Normalized reads in mutant relative to wild-type\n(xf/(xf+N2))\n")

ggsave("figure/known_targets_rpms.violinplot.pdf")
ggsave("figure/known_targets_rpms.violinplot.png")

ggplot(rpms_targets, aes(x=exp, y=RPM_Ratio, fill=exp)) +
   geom_violin() +
   stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) +
   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
   facet_grid(small_RNA_class ~ library_treatment, scales="free_y") +
   scale_fill_manual(values = rev(brewer.pal(7,"Set1"))) +
   theme(strip.text = element_text(size = 14),
      strip.text.y = element_text(size = 10),
      strip.background = element_rect(fill="white")) +
   geom_hline(yintercept=0.5, linetype="dashed", color="gray") +
   xlab("") +
   ylab("Normalized reads in mutant relative to wild-type\n(xf/(xf+N2))\n")

ggsave("figure/known_targets_rpms.violinplotMedian.pdf")
ggsave("figure/known_targets_rpms.violinplotMedian.png")


# ## Sanity checks
# tgt_check <- rpms_c[ geneid %in% targets[targets$exp=="ergo-1",]$geneid & small_RNA_class == "All"]
# summary(tgt_check$RPM_Ratio)
# tgt_check <- rpms_c[ geneid %in% targets[targets$exp=="alg-3_4",]$geneid & small_RNA_class == "All"]
# summary(tgt_check$RPM_Ratio)
# tgt_check <- rpms_c[ geneid %in% targets[targets$exp=="NRDE-3",]$geneid & small_RNA_class == "26G"]
# summary(tgt_check$RPM_Ratio)
# tgt_check


## Overlaps of hits with targets
library(VennDiagram)
library(GeneOverlap)
targets_exps <- unique(targets$exp)
res_files <- list.files(pattern="DESeq2_results.txt")

## background genes
## genes expressed in the germinal line
## http://www.ncbi.nlm.nih.gov/pubmed/25060624
germ_line_genes <- 10754




## plot pval function
# source: https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.text.html
# grid.newpage()
draw.text <- function(txt, just, i, j) {
  grid.text(txt, x=x[j], y=y[i], just=just)
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")

for (i in seq_along(res_files)){ # loop over DESeq2 results
    res_file <- res_files[i]
    print(res_file)
    de_name <- gsub("_DESeq2_results.txt", "", res_file)
    res <- fread(res_file)


    for (j in seq_along(targets_exps)){ # loop over known target experiments
        targets_exp <- targets_exps[j]
        print(targets_exp)
        targets_tmp <- subset(
            targets, exp == targets_exp & !is.na(wormbase_gene)
            )$wormbase_gene
        up <- res[padj < 0.1 & log2FoldChange > 0 ]$wormbase_gene
        down <- res[padj < 0.1 & log2FoldChange < 0 ]$wormbase_gene


        ## significance of overlaps
        over_up <- newGeneOverlap(
            up,
            targets_tmp,
            genome.size=germ_line_genes
            )
        over_down <- newGeneOverlap(
            down,
            targets_tmp,
            genome.size=germ_line_genes
            )
        pval_up <- round(testGeneOverlap(over_up)@pval, 3)
        pval_down <- round(testGeneOverlap(over_down)@pval, 3)
        txt <- paste(
            "Overlap (P-value)", "\n",
            "down-regulated: ", pval_down, "\n",
            "up-regulated: ", pval_up,
            sep="")

        venn.p <- venn.diagram(list(A=targets_tmp,
                    B=up,
                    C=down),
                    fill = c("#461d7e", "#79910f", "#a61a28"),
                    # cex.prop="lin",
                    # alpha = c(0.3, 0.3, 0.3),
        #           cex = 4,
                    main.cex = 2,
                    cat.cex = 1.5,
                    cat.fontface = 2,
        #           lty = 1,
        #           lwd = 4,
                    filename = NULL,
                    category.names = c(targets_exp, "Up-regulated", "Down-regulated"),
                    cat.pos = 0,
                    main = de_name)

        # dev.off()
        # grid.draw(venn.p) # use this if you just want to see plot
        # draw.text(txt, c("left", "top"), 4, 1)
        # ## save pdf
        out_venn <- paste("figure/", targets_exp, "_", de_name, ".venn.pdf", sep="")
        pdf(out_venn)
        grid.draw(venn.p) # use this if you just want to see plot
        draw.text(txt, c("left", "top"), 4, 1)
        dev.off()
    }
    system("rm VennDiagram*.log")
}

system("pdfjam $(ls figure/*21U.venn.pdf) --nup 4x2 --landscape --outfile figure/21U.merged.venn.pdf")
system("pdfjam $(ls figure/*22G.venn.pdf) --nup 4x2 --landscape --outfile figure/22G.merged.venn.pdf")
system("pdfjam $(ls figure/*26G.venn.pdf) --nup 4x2 --landscape --outfile figure/26G.merged.venn.pdf")
system("pdfjam $(ls figure/*All.venn.pdf) --nup 4x2 --landscape --outfile figure/All.merged.venn.pdf")

system("convert -verbose -density 150 -quality 100 -sharpen 0x1.0 figure/21U.merged.venn.pdf figure/21U.merged.venn.png")
system("convert -verbose -density 150 -quality 100 -sharpen 0x1.0 figure/22G.merged.venn.pdf figure/22G.merged.venn.png")
system("convert -verbose -density 150 -quality 100 -sharpen 0x1.0 figure/26G.merged.venn.pdf figure/26G.merged.venn.png")
system("convert -verbose -density 150 -quality 100 -sharpen 0x1.0 figure/All.merged.venn.pdf figure/All.merged.venn.png")
