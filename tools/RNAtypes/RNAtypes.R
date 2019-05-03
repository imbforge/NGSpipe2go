########################################
##
## RNA type barplots based on classes defined in a GTF file
## --
## Who:  Sergi Sayols
## When: 5-nov-2014
## Changed: 28-oct-2015 (Oliver Drechsel)
## --
## Ideally, gencode gtf files should be used, since they contain the most
## comprehensive rna annotation. Check if available for your organism.
## --
## Input:
##   folder=<input folder>
##   pattern=<analyze only files matching this pattern (valid R regular expression)>
##   gtf=<gtf annotation file>
##   out=<output file name>
##   pre=<prefix to be removed from sample name (for plotting)>
##   suf=<suffix to be removed from sample name (for plotting)>
##   paired=<paired end experiment?>
##   stranded=<stranded experiment?>
##   multimap=<how to deal with multimappers (loci and feature)>
##   ftype=<feature type to be counted on the GTF file>
##   ftypecol=<column containing the biotype information in the GTF file>
##   cores=<number of cores to use>
## --
## Run it from bash: double escape special characters
##   $ Rscript RNAtypes.R folder=./test pattern="\\\\.sam$|\\\\.bam$" gtf=./test/test.gtf.gz pre="Sample_imb_richly_2014_06_\\\\d+_" suf="\\\\.bam"
## Run it from bsub: double-double escape special characters
##   $ echo "Rscript RNAtypes.R folder=./test pattern=\\"\\\\\\\\.sam\\$'|'\\\\\\\\.bam\\$\\" gtf=./test/test.gtf.gz pre=\\"Sample_imb_richly_2014_06_\\\\\\\\d+_\\" suf=\\"\\\\\\\\.bam\\"" | bsub
##
########################################

##
## Parse input parms
##
parseArgs <- function(args, string, default=NULL, convert="as.character") {
    if(length(i <- grep(string, args, fixed=T)) == 1) 
        return(do.call(convert, list(gsub(string, "", args[i]))))
    
    if(!is.null(default)) default else do.call(convert, list(NA))
}

args     <- commandArgs(T)
FOLDER   <- parseArgs(args, "folder=", "./")       # folder containing the bam files
PATTERN  <- parseArgs(args, "pattern=", "\\.bam$") # files to be analyzed. Default: all bams
GTF      <- parseArgs(args, "gtf=")                # the GTF file
OUT      <- parseArgs(args, "out=", "RNAtypes")    # output filename
PRE      <- parseArgs(args, "pre=", "")            # pattern to remove from the file name
SUF      <- parseArgs(args, "suf=", "\\.bam$")     # pattern to remove from the file name
PAIRED   <- parseArgs(args, "paired=", "no")       # strand specific assay
STRANDED <- parseArgs(args, "stranded=", "no")     # strand specific assay
MMAPPERS <- parseArgs(args, "multimap=", "NONE")   # multimappers (loci and feature): ALL|UNAMBIGUOUS|NONE|RANDOM
FTYPE    <- parseArgs(args, "ftype=", "exon")      # feature type on the GTF file
FTYPECOL <- parseArgs(args, "ftypecol=", "gene_type")  # column containing the biotype
CORES    <- parseArgs(args, "cores=", 1, "as.numeric") # number of cores to use

print(args)
if(length(args) == 0 | args[1] == "-h" | args[1] == "--help")
    stop(paste("Rscript RNAtypes.R [arguments|-h|--help]\n", 
               "  [folder=./]         : input folder\n", 
               "  [pattern=\"\\\\.bam$\"] : analyze only files matching this pattern (valid R regular expression)\n", 
               "  <gtf=gencode.gtf>   : gtf annotation file. Can be a compressed file\n", 
               "  [out=RNAtypes]      : output filename\n", 
               "  [pre=\"\"]          : prefix to be removed from sample name (for plotting)\n", 
               "  [suf=\"\\\\.bam$\"] : suffix to be removed from sample name (for plotting)\n", 
               "  [paired=no]         : no|yes\n", 
               "  [stranded=no]       : no|yes|reverse\n", 
               "  [multimap=NONE]     : multimappers (loci and feature): ALL|UNAMBIGUOUS|NONE|RANDOM\n", 
               "  [ftype=exon]        : feature type on the GTF file to count on (exon|gene|...)\n", 
               "  [ftypecol=gene_type]: column name in the GTF file containing the biotype info\n", 
               "  [cores=1]           : number of cores to use"))
if(!grepl("ALL|UNAMBIGUOUS|NONE|RANDOM", MMAPPERS)) stop("multimap must be ALL|UNAMBIGUOUS|NONE|RANDOM")
if(is.na(CORES))      stop("cores has to be an integer number")
if(is.na(GTF))        stop("gtf argument is mandatory")
if(!file.exists(GTF)) stop(paste("File", GTF, "does NOT exist"))
if(!file.exists(FOLDER)) stop(paste("Dir", FOLDER, "does NOT exist"))
if(is.na(PAIRED)   | !(grepl("no|yes", PAIRED))) stop("paired has to be no|yes")
if(is.na(STRANDED) | !(grepl("no|yes|reverse", STRANDED))) stop("stranded has to be no|yes|reverse")
files <- list.files(path=FOLDER, pattern=PATTERN)
files <- files[grep(PATTERN, files)]
samples <- gsub(PRE, "", gsub(SUF, "", files))
if(length(files) == 0) stop(paste("Dir", FOLDER, "does not contain files matching the pattern. Nothing to do"))
if(!dir.exists(dirname(OUT))) dir.create(dirname(OUT), showWarnings=FALSE, recursive=TRUE)

##
## load libraries
##
library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(reshape)
library(ggplot2)
library(parallel)

##
## Read in the GTF file
##
gencode <- import.gff(GTF, format="gtf", feature.type=FTYPE)
if(!any(grepl(FTYPECOL, colnames(elementMetadata(gencode))))) stop("GTF file doesn't contain info about the biotype")

# extract unique gene biotypes
gene.biotypes <- elementMetadata(gencode)[, c("gene_id", FTYPECOL)]
gene.biotypes <- unique(gene.biotypes, MARGIN=1)

# flatten GTF file and annotate witht he biotypes
gencode.flat  <- unlist(reduce(split(gencode, elementMetadata(gencode)$gene_id)))
elementMetadata(gencode.flat)$gene_biotype <- gene.biotypes[match(names(gencode.flat), gene.biotypes$gene_id), FTYPECOL]

##
## Read in the bam file and intersect with the annotation
##
d <- mclapply(1:length(files), function(i) {

    # read bam file
    f <- files[i]
    aln <- switch(PAIRED, 
                  no =readGAlignments    (paste0(FOLDER, "/", f), param=ScanBamParam(tag="NH")), 
                  yes=readGAlignmentPairs(paste0(FOLDER, "/", f), param=ScanBamParam(tag="NH")))

    # reverse strand if required
    if(STRANDED == "reverse" & PAIRED == "no" ) aln <- invertStrand(aln)
    if(STRANDED == "reverse" & PAIRED == "yes") strandMode(aln) <- switch(STRANDED, no=0, yes=1, reverse=2)

    # match reads with gencode.flat features
    reads <- switch(MMAPPERS, 
        ALL  = { 
            mm <- findOverlaps(aln, gencode.flat, ignore.strand=(STRANDED == "no"))
            table(gencode.flat$gene_biotype[as.matrix(mm)[, 2]])
        }, 
        UNAMBIGUOUS = {
            hits <- countOverlaps(aln, gencode.flat, ignore.strand=(STRANDED == "no"))
            mm   <- findOverlaps(aln[hits == 1], gencode.flat, ignore.strand=(STRANDED == "no"))
            table(gencode.flat$gene_biotype[as.matrix(mm)[, 2]])
        }, 
        NONE = {
            aln  <- switch(PAIRED, 
                           no =aln[elementMetadata(aln)$NH == 1, ], 
                           yes=aln[elementMetadata(first(aln))$NH == 1, ])    # both mates (first and last) have the same NH
            hits <- countOverlaps(aln, gencode.flat, ignore.strand=(STRANDED == "no"))
            mm   <- findOverlaps(aln[hits == 1], gencode.flat, ignore.strand=(STRANDED == "no"))
            table(gencode.flat$gene_biotype[as.matrix(mm)[, 2]])
        }, 
        RANDOM = {
            mm <- findOverlaps(as(aln, "GRanges"), gencode.flat, ignore.strand=(STRANDED == "no"), 
                               select="arbitrary", type="within")
            table(gencode.flat$gene_biotype[mm[!is.na(mm)]])
        })

    # normalize RPK of feature length
    feat.width <- sapply(names(reads), function(x) sum(width(gencode.flat[gencode.flat$gene_biotype == x])))
    rpk <- round(reads * 10^3 / feat.width, 2)

    # return
    return(rbind(reads, rpk))

}, mc.cores=CORES)

## merge the tables from all the files
rmerge <- function(x, y) { 
    suppressWarnings(
    z <- merge(t(x[grepl(what, rownames(x)), , drop=F]), 
               t(y[grepl(what, rownames(y)), , drop=F]), by=0, all=T)
    )
    rownames(z) <- z[, 1]
    t(as.matrix(z[, -1]))
}

what    <- "reads"
d.reads <- Reduce(rmerge, d)
rownames(d.reads) <- samples

what  <- "rpk"
d.rpk <- Reduce(rmerge, d)
rownames(d.rpk) <- samples

##    
## barplot
##

# plot total counts
p  <- function(x, main, normalize=T) {
    
    # calculate the % if requested, and melt the df into ggplot format
    if(normalize) x <- t(apply(x, 1, function(x) 100 * x / sum(x, na.rm=T)))
    x <- melt(x)
    x$value <- ifelse(is.na(x$value), 0, x$value)

    ggplot(x, aes(y=value, x=X1, fill=X2)) + 
        geom_bar(stat="identity") + 
        ggtitle(main) +
        theme_bw() + 
        theme(text=element_text(size=14), axis.text.x=element_text(angle=90, vjust=0.5)) + 
        labs(x="samples") + 
        guides(fill=guide_legend(ncol=3, title="feature type", title.position="top", title.hjust=.5))
}

# plots
png(paste0(OUT, ".counts.raw.png"), height=700, width=1400)
write.csv(d.reads, file=paste0(OUT, ".counts.raw.csv"))
print(p(d.reads, "Counts (raw)", normalize=F))
dev.off()

png(paste0(OUT, ".counts.per.png"), height=700, width=1400)
print(p(d.reads, "Counts (percentage)", normalize=T))
dev.off()

png(paste0(OUT, ".counts.rpk.png"), height=700, width=1400)
write.csv(d.rpk, file=paste0(OUT, ".counts.rpk.csv"))
print(p(d.rpk  , "RPK", normalize=F))
dev.off()
