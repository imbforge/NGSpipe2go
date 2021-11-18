
## load libraries
library(GenomicRanges)
library(rtracklayer)
options(stringsAsFactors=FALSE)

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
   
   if(length(i <- grep(string,args,fixed=T)) == 1)
      return(do.call(convert,list(gsub(string,"",args[i]))))
   
   if(!is.null(default)) default else do.call(convert,list(NA))
}

## read input arguments
args <- commandArgs(T)
gene.model <- parseArgs(args,"gtf=","")       # gtf gene model
count.file <- parseArgs(args,"input=","")  
output.dir <- parseArgs(args,"outdir=","")

## check input arguments
runstr <- "Rscript extract_miRNA_snoRNA2.R [gtf=] [input=]"
if(!file.exists(gene.model)) stop(paste("GTF File:", gene.model, " does NOT exist. Run with: \n", runstr))
if(!file.exists(count.file)) stop(paste("input count file:", count.file, " does NOT exist. Run with: \n", runstr))
if(!dir.exists(output.dir))  stop(paste("output directory:", output.dir, " does NOT exist. Run with: \n", runstr))

## output file names
mirna.output.file  <- paste0(output.dir,"/",gsub(".readcounts.tsv",".miRNA.readcounts.tsv",basename(count.file)))

## read counts
counts.df <- read.delim(count.file,header=FALSE,sep="\t")
names(counts.df) <- c("gene_id","count")

## read in annotation
gtf <- import.gff(gene.model, format="gtf", feature.type="gene")

## get different name information
gtf.names.df <- data.frame(gene_id=gtf$gene_id,gene_type=gtf$gene_type)

## merge counts.df and gtf.names.df
counts.names.df <- merge(counts.df,gtf.names.df,by="gene_id",sort=F)

## write to file
write.table(counts.names.df[counts.names.df$gene_type=="miRNA",c("gene_id","count")],
            file=mirna.output.file,
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)




