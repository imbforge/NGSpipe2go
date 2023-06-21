#####################################
##
## What: sc_readAggrData.R
## Who : Frank Rühle
## When: 24.04.2023
##
## Script to import 10X multiome assay data.
##
## Args:
## -----
## targets=targets.txt      # file describing the targets 
##
##
######################################

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args               <- commandArgs(T)
projectdir         <- parseArgs(args,"project=") 
ftargets           <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
gene.model         <- parseArgs(args,"gtf=","")       # gtf gene model
resultsdir         <- parseArgs(args,"res=")   
out                <- parseArgs(args,"outdir=") # output folder
org                <- parseArgs(args,"org=")   # organism name needed for cell cycle and GO enrichment 
db                 <- parseArgs(args,"db=")     
MTgenes            <- parseArgs(args,"mtgenes=",NA) 
selgenes           <- parseArgs(args,"selgenes=","") 
cellranger_aggr_id <- parseArgs(args, "cellranger_aggr_id=") 
run_demux          <- parseArgs(args, "run_demux=") 
demux_out          <- parseArgs(args, "demux_out=") 
demuxcluster_out   <- parseArgs(args, "demuxCluster_out=") 
colorbyfactor      <- parseArgs(args, "colorByFactor=", "") 

if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(gene.model)) stop(paste("File",gene.model,"does NOT exist. Run with:\n",runstr))

runstr <- "Rscript sc_readAggrData.R [targets=targets.txt] [gtf=]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_readData/renv.lock"))
print(.libPaths())

library(tidyverse)
library(AnnotationDbi)
library(Biobase)
library(data.table)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(reshape2)
library(scater)
library(scran)
library(scuttle)
library(Seurat)
library(Signac)
library(uwot)

# set options
options(stringsAsFactors=FALSE)

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("ftargets:", ftargets))
print(paste("gene.model:", gene.model))
print(paste("org:", org))
print(paste("db:", db))
print(paste("MTgenes:", MTgenes))
print(paste("selgenes:", selgenes))
print(paste("cellranger_aggr_id:", cellranger_aggr_id))
print(paste("run_demux:", run_demux))
print(paste("demux_out:", demux_out))
print(paste("demuxcluster_out:", demuxcluster_out))
print(paste("colorbyfactor:", colorbyfactor))


# load list of mitochondrial genes (if not given, we will use all genes starting with "MT-")
mito.genes <- if(file.exists(file.path(projectdir,MTgenes))) {
  read.delim(file.path(projectdir,MTgenes))[, 1]
} else {NA}

# load relevant BSgenome package (needed by Signac for motif analysis)
switch(db,
       hg38={ failed_BSgenome <- library("BSgenome.Hsapiens.UCSC.hg38") },
       mm10={ failed_BSgenome <- library("BSgenome.Mmusculus.UCSC.mm10") },
       stop(c("Don't find genome:", db))   
)


aggrcsv <- read.delim(file.path(resultsdir, "aggr.csv"), sep=",") # aggr.csv contains the order of GEM wells
# load targets.txt (data per fastq file). Columns required: "sample", "file", "group", "replicate"
targets_pools <- read.delim(ftargets, sep="\t", comment.char = "#")
# create targets on cell level
cellranger_dir <- file.path(resultsdir, cellranger_aggr_id, "outs")
matrix_dir <- file.path(cellranger_dir, "filtered_peak_bc_matrix")                                    # 10X ATAC
if (!dir.exists(matrix_dir)) {matrix_dir <- file.path(cellranger_dir, "filtered_feature_bc_matrix")}  # 10X Multiome
targets = read.delim(file.path(matrix_dir, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
names(targets) = "cell_id"

if(!is.na(run_demux)) { # sample multiplexing applied?
  
  if(run_demux=="demux_HTO") { # add HTO information if cell hashing applied
    demuxfiles <- lapply(1:nrow(aggrcsv), function(x) {
      read.delim(file.path(demux_out, gsub("_S._L..._R._001", "", aggrcsv$library_id)[x], "Seurat", "demux.txt"), sep="\t")
    })
    demuxfiles <- dplyr::bind_rows(demuxfiles, .id = "GEMwell")
    demuxfiles$cell_id <- paste(demuxfiles$cell_id, demuxfiles$GEMwell, sep="-")
    targets <- merge(targets, demuxfiles[, c("cell_id", "GEMwell", "HTO_maxID", "HTO_classification.global", 
                                             "nCount_HTO", "nFeature_HTO", "HTO_secondID", "HTO_margin")], by="cell_id")
    targets$name_HTO <- gsub("-.*$", "", targets$HTO_maxID)
  }
  
  if(run_demux=="demux_GT") { # add information if sample multiplexing by genetic variance is applied
    
    # load assignSouporcellCluster files to correlate cluster among samples
    clusterassignFilenames <- list.files(demuxcluster_out, pattern = "\\.txt$", full.names=T)
    clusterassign <- lapply(clusterassignFilenames, read.delim, sep="\t", skip=5)
    clustercount  <- as.numeric(gsub("clusters for experiment1 ", "", sapply(clusterassignFilenames, read.delim, skip=2, nrows=1)))
    names(clusterassign) <- names(clustercount) <- gsub("\\.txt$", "", basename(clusterassignFilenames))
    for (i in names(clusterassign)) {
      colnames(clusterassign[[i]])[1:2] <- unlist(strsplit(i, split="_vs_"))
      clusterassign[[i]] <- clusterassign[[i]][1:clustercount[i],1:2]
    }
    
    # load demux_GT files in order defined in aggrcsv
    aggrcsvPrunedSampleId <- gsub("_S._L..._R._001", "", aggrcsv$library_id)
    demuxfiles <- lapply(1:nrow(aggrcsv), function(x) { # files loaded in order defined in aggrcsv
      read.delim(file.path(demux_out, aggrcsvPrunedSampleId[x], "clusters.tsv"), sep="\t")
    })
    
    # print table with status column later
    overviewStatus <- sapply(demuxfiles, function(x) table(x$status))
    colnames(overviewStatus) <- aggrcsvPrunedSampleId
    overviewStatus <- as.data.frame(dplyr::bind_rows(as.data.frame(overviewStatus), .id = "file")) ########### CHECK WHAT's WRONG
    
    # continue with singlets only
    demuxfiles <- lapply(demuxfiles, function(x) {
      x <- x[x$status=="singlet",]
      x$subsample <- x$assignment # column used for re-assignment of souporcell clusters
      return(x)})
    
    # Replace cluster assignments for all samples in demuxfiles except for the first one which is used as reference
    for (j in 2:length(aggrcsvPrunedSampleId)) {
      clusterassignCurrent <- names(clusterassign)[grepl(aggrcsvPrunedSampleId[1], names(clusterassign)) & grepl(aggrcsvPrunedSampleId[j], names(clusterassign))]
      newValues <- clusterassign[[clusterassignCurrent]][,aggrcsvPrunedSampleId[1]]
      oldValues <- clusterassign[[clusterassignCurrent]][,aggrcsvPrunedSampleId[j]]
      demuxfiles[[j]]$subsample <- plyr::mapvalues(demuxfiles[[j]]$assignment, from=oldValues, to=newValues) # new column with adjusted assignment
    }
    
    # bind rows to single dataframe and filter cell ids by merging with targets
    demuxfiles <- dplyr::bind_rows(demuxfiles, .id = "GEMwell")
    demuxfiles$cell_id <- paste(gsub("-.$", "", demuxfiles$barcode), demuxfiles$GEMwell, sep="-")
    targets <- merge(targets, demuxfiles[, c("cell_id", "GEMwell", "subsample", "assignment", "status", "log_prob_singleton", "log_prob_doublet")], by="cell_id")
    
    # If subsample column in targets.txt is not pre-defined by user, it is automatically generated.
    # This is used for merging targets with targets_pools to link the cell barcodes with sample names
    # but note that in that case the assignment of subsamples to the sample ids in targets_pools is arbitrary 
    # (need external genetic information for specific assignment).
    if(is.null(targets_pools$subsample) || all(is.na(targets_pools$subsample))) {
      targets_pools <- targets_pools %>% dplyr::group_by(file)  %>%
        dplyr::mutate(subsample = 0:(dplyr::n_distinct(sample)-1)) %>% # enumerate samples per file (start with 0)
        dplyr::ungroup() %>%
        as.data.frame()
    }
  }
}

# merge targets with targets_pools and sort 
# ensure that GEM wells merged in aggr module are assigned the sample names from targets.txt in correct order
# the -digit suffix of the cell barcode reflects the order that the GEM wells were provided in the aggr.csv
targets$file <- gsub("_S._L..._R._001", "", aggrcsv[as.numeric(gsub("^.*-", "", targets$cell_id)), "library_id"])
targets_pools$file <- gsub("(_S.)*(_L...)*(_R.)*(_001)*", "", targets_pools$file) # in case file name was partly pruned in targets_pools
targets <- merge(targets, targets_pools) # if "name_HTO" column is present in targets, it is used for merging as well (unmapped cells discarded)
rownames(targets) <- targets$cell_id
StartingColumnNames <- c("cell_id", "sample", "file", "group", "replicate")
targets <- targets[order(targets$sample, targets$cell_id), c(StartingColumnNames, colnames(targets)[!colnames(targets) %in% StartingColumnNames])]

# define group.vars
demux_columns <- c("assignment", "status", "log_prob_singleton", "log_prob_doublet",
                   "GEMwell","HTO_maxID","nCount_HTO","HTO_secondID", "HTO_margin", "file_HTO","seq_HTO")
group.vars <- colnames(targets)[!colnames(targets) %in% c("cell_id", "sample", "file", demux_columns)] 
targets[,group.vars] <- lapply(targets[,group.vars], factor)
group.vars <- group.vars[sapply(group.vars, function(x) {length(unique(targets[,x]))>1})] # remove group vars with single value
# in case of hto multiplexing targets_pools is collapsed to (not de-multiplexed) fastq files for the fastq-file level QC plots (targets4plots).
# But 'group' and 'replicate' columns are not kept if they show multiple values per file (e.g. multiplexed samples belong to different groups). 
# This is only relevant for the gene body coverage plots per group or per replicate, which would be misleading in that case.
collapsedColumns <- sapply(c("sample", "group", "replicate"), function(x) {
  all(sapply(unique(targets_pools$file), function(y) {length(unique(targets_pools[targets_pools$file==y, x]))==1}))
})
targets4plots <- targets_pools[!duplicated(targets_pools$file), c("file", names(collapsedColumns)[collapsedColumns])] 
if(!"sample" %in% colnames(targets4plots)) {targets4plots$sample <- targets4plots$file} # "sample" column is needed. Same as "file" if ambiguous.
colorByFactor <- "sample" # default for pipeline overview plots on file level
colorByFactor2 <- if(!is.na(colorbyfactor)) {colorbyfactor} else {"sample"} # default for downstream plots on cell level (up to 2 categories)


#### ######################
# load gene annotation provided in essential.vars.groovy
gtf <- import.gff(gene.model, format="gtf") # , feature.type="exon"

# Signac requires 'gene_biotype' field in annotation GTF, 10X (GENCODE) uses 'gene_type'
if ("gene_type" %in% colnames(mcols(gtf))) {
  gtf$gene_biotype <- gtf$gene_type
}

# Other changes necessary for Signac to properly parse the 10X GTF
genome(gtf) <- db
seqlevelsStyle(gtf) <- "UCSC"
gtf <- keepStandardChromosomes(gtf, pruning.mode = "coarse")

names(gtf) <- gtf$gene_name



#############################
# load the RNA and ATAC data
counts <- Read10X_h5(file.path(cellranger_dir, "filtered_feature_bc_matrix.h5"))
fragpath <- file.path(cellranger_dir, "atac_fragments.tsv.gz")

# library(EnsDb.Hsapiens.v86)
# anno <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"

# filter out genes not annotated in gtf
genes.use <- rownames(counts$`Gene Expression`) %in% names(gtf)
counts$`Gene Expression` <- counts$`Gene Expression`[genes.use,]

# create a Seurat object containing the RNA adata
sobj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  #  annotation = gtf, # why no effect? Can't I load different anno for ATAC peaks?
  assay = "RNA"
)


# create ATAC assay and add it to the object
# blacklist_hg38_unified comes with Signac package
# we'll only use peaks in standard chromosomes and not overlapping with blacklisted regions
grange.counts <- StringToGRanges(rownames(counts$Peaks), sep = c(":", "-"))
grange.counts <- subsetByOverlaps(grange.counts, ranges = blacklist_hg38_unified, invert = TRUE) 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
counts$Peaks <- counts$Peaks[as.vector(grange.use), ]

sobj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks, # creates matrix with read counts not fragments
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = gtf # 
)


if(!all(rownames(targets) %in% colnames(sobj))) {stop("There are cells in your targets list which are not in the Seurat object!")}
#table(colnames(sobj) %in% rownames(targets))

sobj <- sobj[,colnames(sobj) %in% rownames(targets)]
targets <- targets[match(colnames(sobj), rownames(targets)),]

#sobj <- AddMetaData(object = sobj, metadata = targets$sample, col.name = 'sample')
sobj <- AddMetaData(
  object = sobj,
  metadata = targets
)


# Add fragment counts (e.g. for counting reads in peaks)
total_fragments <- CountFragments(fragpath)
rownames(total_fragments) <- total_fragments$CB
sobj$fragments <- total_fragments[colnames(sobj), "frequency_count"]





#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/sc_readAggrData_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(gtf, file=paste0(resultsdir,"/gtf.RData"))
save(targets, targets4plots, colorByFactor, colorByFactor2, fragpath, cellranger_dir, 
     file=paste0(out,"/sc_readAggrData.RData"))

