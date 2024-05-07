#####################################
##
## What: sc_readIndivSamplesAndMerge.R
## Who : Sivarajan Karunanithi
## When: 07.11.2023
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
cellranger_output  <- parseArgs(args,"cellranger_output=","") 
run_demux          <- parseArgs(args, "run_demux=") 
demux_out          <- parseArgs(args, "demux_out=") 
demuxcluster_out   <- parseArgs(args, "demuxCluster_out=") 
colorbyfactor      <- parseArgs(args, "colorByFactor=", "") 

if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(gene.model)) stop(paste("File",gene.model,"does NOT exist. Run with:\n",runstr))

runstr <- "Rscript sc_readIndivSamplesAndMerge.R [targets=targets.txt] [gtf=]"

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
print(paste("cellranger_output:", cellranger_output))
print(paste("run_demux:", run_demux))
print(paste("demux_out:", demux_out))
print(paste("demuxcluster_out:", demuxcluster_out))
print(paste("colorbyfactor:", colorbyfactor))



# load relevant BSgenome package (needed by Signac for motif analysis)
switch(db,
       hg38={ failed_BSgenome <- library("BSgenome.Hsapiens.UCSC.hg38") },
       mm10={ failed_BSgenome <- library("BSgenome.Mmusculus.UCSC.mm10") },
       stop(c("Don't find genome:", db))   
)

# load targets.txt (data per fastq file). Columns required: "sample", "file", "group", "replicate"
targets_pools <- read.delim(ftargets, sep="\t", comment.char = "#")
targets_pools$file <- gsub("(_S.)*(_L...)*(_R.)*(_001)*", "", targets_pools$file) # in case file name was partly pruned in targets_pools

# this will determine the order of samples loaded, and will affect the cell barcode extension (-1, -2, ...)
filenames = unique(sort(targets_pools$file))

if(!is.na(run_demux)) { # sample multiplexing applied?

  # demux_HTO is not tested
  if(run_demux=="demux_HTO") { # add HTO information if cell hashing applied
    demuxfiles <- lapply(1:length(filenames), function(x) {
      read.delim(file.path(demux_out, filenames[x], "Seurat", "demux.txt"), sep="\t")
    })
    demuxfiles <- dplyr::bind_rows(demuxfiles, .id = "GEMwell")
    demuxfiles$cell_id <- paste(demuxfiles$cell_id, demuxfiles$GEMwell, sep="-")
    targets <- demuxfiles[, c("cell_id", "GEMwell", "HTO_maxID", "HTO_classification.global",
                                             "nCount_HTO", "nFeature_HTO", "HTO_secondID", "HTO_margin")]
    targets$name_HTO <- gsub("-.*$", "", targets$HTO_maxID)
  }

  if(run_demux=="demux_GT") {

    # load assignSouporcellCluster files to correlate cluster among samples
    clusterassignFilenames <- list.files(demuxcluster_out, pattern = "\\.txt$", full.names=T)
    clusterassign <- lapply(clusterassignFilenames, read.delim, sep="\t", skip=5)
    clustercount  <- as.numeric(gsub("clusters for experiment1 ", "", sapply(clusterassignFilenames, read.delim, skip=2, nrows=1)))
    names(clusterassign) <- names(clustercount) <- gsub("\\.txt$", "", basename(clusterassignFilenames))
    for (i in names(clusterassign)) {
      colnames(clusterassign[[i]])[1:2] <- unlist(strsplit(i, split="_vs_"))
      clusterassign[[i]] <- clusterassign[[i]][1:clustercount[i],1:2]
    }

    # load demux_GT files in order defined in filenames
    demuxfiles <- lapply(1:length(filenames), function(x) {
      read.delim(file.path(demux_out, filenames[x], "clusters.tsv"), sep="\t")
    })

    # continue with singlets only
    demuxfiles <- lapply(demuxfiles, function(x) {
      x <- x[x$status=="singlet",]
      x$subsample <- x$assignment # column used for re-assignment of souporcell clusters
      return(x)})

    # Replace cluster assignments for all samples in demuxfiles
    # except for the first one which is used as reference
    for (j in 2:length(filenames)) {
      clusterassignCurrent <- names(clusterassign)[grepl(filenames[1], names(clusterassign)) &
                                                     grepl(filenames[j], names(clusterassign))]
      newValues <- clusterassign[[clusterassignCurrent]][,filenames[1]]
      oldValues <- clusterassign[[clusterassignCurrent]][,filenames[j]]
      demuxfiles[[j]]$subsample <- plyr::mapvalues(demuxfiles[[j]]$assignment,
                                                   from=oldValues, to=newValues) # new column with adjusted assignment
    }

    # bind rows to single dataframe
    demuxfiles <- dplyr::bind_rows(demuxfiles, .id = "GEMwell")

    # Need to modify the cellid with a suffix representing which sample the barcode comes from
    # with out that, there will be duplicate barcodes as same cell barcodes can come from different sample.
    # As we need unique barcodes to use it as rownames for a metadata (targets, we call it)
    # which will be added as metadata of the seurat object later.
    demuxfiles$cell_id <- paste(gsub("-.$", "", demuxfiles$barcode), demuxfiles$GEMwell, sep="-")

    # for consistency reasons (with the loading of cellrange aggr results- sc_readAggrData.R)
    # a subset of columns from the demuxfiles will be called targets, which will be later merged
    # based on the cellid with the barcodes used in the analysis.
    # We keep the default barcodes (ending with -1) as well, for merging with cellranger count output
    targets <- demuxfiles[, c("cell_id", "barcode", "GEMwell", "subsample", "assignment",
                              "status", "log_prob_singleton", "log_prob_doublet")]


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
# As the GEMwells are created based on the alphabetical order of filenames
# I use that to populate the file name.
targets$file = filenames[as.numeric(targets$GEMwell)]
targets <- merge(targets, targets_pools) # if "name_HTO" column is present in targets, it is used for merging as well (unmapped cells discarded)
rownames(targets) <- targets$cell_id
StartingColumnNames <- c("cell_id", "sample", "file", "group", "replicate")
targets <- targets[order(targets$sample, targets$cell_id),
                   c(StartingColumnNames, colnames(targets)[!colnames(targets) %in% StartingColumnNames])]

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

###### Reading RNA/ATAC data #######

# Reading again the ftargets as we need the full folder names from the targets file to get the filtered matrices
# targets_filenames <- read.delim(ftargets, sep="\t", comment.char = "#")
# filenames = unique(sort(targets_filenames$file))

# We need the full folder names from the targets file to get the filtered matrices
count_outputfolders=list.dirs(path = cellranger_output, full.names = FALSE, recursive = FALSE)

allsamples_preprocess = lapply( c(1:length(filenames)), function(x){

  # In case the cellranger folder names are different (like including the _S1_L001...) 
  # from the file column of targets
  h5fileloc = file.path(cellranger_output, 
                        count_outputfolders[grepl(filenames[x], count_outputfolders)], 
                        "outs", "filtered_feature_bc_matrix.h5"
                        )
  # Read the counts - contains both peaks and RNA counts
  counts <- Read10X_h5(h5fileloc)

  # Read the counts - contains both peaks and RNA counts
  # counts <- Read10X_h5(file.path(cellranger_output, filenames[x] , "outs", "filtered_feature_bc_matrix.h5"))

  # filter out genes not annotated in gtf
  genes.use <- rownames(counts$`Gene Expression`) %in% names(gtf)
  counts$`Gene Expression` <- counts$`Gene Expression`[genes.use,]

  # create a Seurat object containing the RNA data
  sobj_rna <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )

  # We will merge (reduce) the peaks to create a combined peaks for all the samples later
  peaks_info = StringToGRanges(rownames(counts$Peaks), sep = c(":", "-"))

  list("sobj_rna" = sobj_rna, "peaks" = peaks_info)
})

gc()

# create a merged list of peaks from all samples
reduced_peaks = reduce(unlist(GRangesList(lapply(allsamples_preprocess, "[[", "peaks"))))

# Filter out bad peaks based on length
peakwidths <- width(reduced_peaks)
reduced_peaks <- reduced_peaks[peakwidths  < 10000 & peakwidths > 20]
# blacklist_hg38_unified comes with Signac package
# we'll only use peaks in standard chromosomes and not overlapping with blacklisted regions
reduced_peaks <- subsetByOverlaps(reduced_peaks, ranges = blacklist_hg38_unified, invert = TRUE)
reduced_peaks.use <- seqnames(reduced_peaks) %in% standardChromosomes(reduced_peaks)

# Create the RNA and ATAC merged seurat objects
sample_sobjs = lapply( c(1:length(filenames)), function(x){
  # Read and prepare the fragments for quanitifying the merged peaks
  # In case the cellranger folder names are different (like including the _S1_L001...) 
  # from the file column of targets
  fragpath = file.path(cellranger_output, 
                        count_outputfolders[grepl(filenames[x], count_outputfolders)], 
                        "outs", "atac_fragments.tsv.gz"
  )
  # fragpath <- file.path(cellranger_output, filenames[x] , "outs", "atac_fragments.tsv.gz")
  frags_info <- CreateFragmentObject(path = fragpath, cells = colnames(allsamples_preprocess[[x]]$sobj_rna))
  
  # measure the signal using the merged peak list
  peakcounts <- FeatureMatrix(fragments = frags_info,
                              features = reduced_peaks,
                              cells = colnames(allsamples_preprocess[[x]]$sobj_rna),
                              sep = c("-","-")
  )
  # create an atac seurat object
  sobj_atac <- CreateSeuratObject(CreateChromatinAssay(
    peakcounts[as.vector(reduced_peaks.use), ],
    fragments = frags_info,
    annotation = gtf
  ),
  assay = "ATAC"
  )
  
  # Let us merge the RNA and ATAC assays
  sobj = allsamples_preprocess[[x]]$sobj_rna
  sobj[["ATAC"]] = sobj_atac[["ATAC"]]
  
  # For some reason, if I first add this metadata to ATAC only object and then copy the ATAC to the RNA+ATAC CreateSeuratObject
  # the "fragments" column doesn't get copied.
  # We need to quantify the total fragments for each sample to use later
  total_fragments_info <- CountFragments(fragments=fragpath)
  rownames(total_fragments_info) = total_fragments_info$CB
  # Add additional metadata - fragment frequency (total number of fragments sequenced for the cell)
  sobj = AddMetaData(object=sobj, metadata=total_fragments_info[colnames(sobj), "frequency_count"], col.name='fragments')
  
  # will append the cellid (-1,-2,..) in the same order as processed in the demux section
  sobj = RenameCells(object=sobj, new.names=paste(gsub("-.$", "", colnames(sobj)), x, sep="-"))
}
)

gc()

# Merge all datasets
sobj <- merge(x = sample_sobjs[[1]],
              y = sample_sobjs[2:length(sample_sobjs)])


if(!all(rownames(targets) %in% colnames(sobj))) {stop("There are cells in your targets list which are not in the Seurat object!")}

sobj <- sobj[,colnames(sobj) %in% rownames(targets)]
targets <- targets[match(colnames(sobj), rownames(targets)),]

#sobj <- AddMetaData(object = sobj, metadata = targets$sample, col.name = 'sample')
sobj <- AddMetaData(
  object = sobj,
  metadata = targets
)

# Calculate the percentage of mitochondrial counts (either load MT gene list or match feature pattern (e.g. "MT-"))
if(!is.na(MTgenes)) {
  if(file.exists(file.path(projectdir,MTgenes))) {
    mito.genes <- read.delim(file.path(projectdir,MTgenes))[, 1]
    # store mitochondrial percentage in object meta data
    sobj <- PercentageFeatureSet(sobj, assay = "RNA", features=mito.genes, col.name = "percent.mt")
  } else {
    # store mitochondrial percentage in object meta data
    sobj <- PercentageFeatureSet(sobj, assay = "RNA", pattern = paste0("^", MTgenes), col.name = "percent.mt")
  }
}


#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/sc_readIndivSamplesAndMerge_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(gtf, file=paste0(resultsdir,"/gtf.RData"))
save(targets, targets4plots, colorByFactor, colorByFactor2,
     file=paste0(out,"/sc_readIndivSamplesAndMerge.RData"))

