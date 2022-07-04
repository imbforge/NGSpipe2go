#####################################
##
## What: demux_hto.R
## Who : Frank Rühle
## When: 27-04-2022
##
## Script to perform sample demultiplexing for 10X single cell experiments involving cell hashing 
## with hashtag oligos (HTO). This scripts uses the output generated from CITE-seq-Count and
## outputs a table with respective cellbarcode - hashtag combinations.
##
######################################

options(stringsAsFactors=FALSE)
library(ggplot2)
library(Seurat)


##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
excludeFailedHTO <- parseArgs(args,"excludeFailedHTO=","")
out              <- parseArgs(args,"out=","Seurat") # output filename
hto_matrix_dir   <- parseArgs(args,"hto_matrix_dir=","umi_count") # output from Cite-Seq-Count
rna_matrix_dir   <- parseArgs(args,"rna_matrix_dir=","") # output from Cellranger

runstr <- "Rscript demux_hto.R [excludeFailedHTO=] [out=Seurat] [hto_matrix_dir=umi_count] [rna_matrix_dir=]"
if(!dir.exists(hto_matrix_dir))   stop(paste("Directory",hto_matrix_dir,"does NOT exist. Run with:\n",runstr))
if(!dir.exists(rna_matrix_dir))   stop(paste("Directory",rna_matrix_dir,"does NOT exist. Run with:\n",runstr))


if (!dir.exists(out)) {dir.create(out)}

## read in data with Seurat

cat("\nLoad in hashtag data\n")
hto <- Seurat::Read10X(hto_matrix_dir, gene.column=1, strip.suffix = T)
cat("Counts per Hashtag (if a hashtag did not work at all, it should be removed from downstream analysis):\n")
print(apply(hto, 1, sum))
if(excludeFailedHTO != "") { # remove failed hashtags if given
  cat(paste("Remove hashtags:", paste(excludeFailedHTO, collapse=", "), "\n"))
  hto <- hto[!rownames(hto) %in% excludeFailedHTO,]
}
# Confirm that the HTO have the correct names
cat(paste("hashtag names used in downstream analysis:", paste(rownames(hto), collapse=", "), "\n"))


cat("\nLoad in rna data\n")
rna <- Seurat::Read10X(rna_matrix_dir, gene.column=1, strip.suffix = T) # gene.column=1 ensemblid, 2 symbol
# remove features with too few counts
threshold4features <- 1
htsum <- apply(rna, 1, sum)
htsum0 <- htsum[htsum<=threshold4features]
cat(paste("Remove", sum(htsum0), "features with <=", threshold4features, "read counts.\n"))
rna <- rna[!(rownames(rna) %in% names(htsum0)),]
# remove cell barcodes with too few counts
threshold4bc <- 1
bcsum <- apply(rna, 2, sum)
bcsum0 <- bcsum[bcsum<=threshold4bc]
cat(paste("Remove", sum(bcsum0), "cell barcodes with <=", threshold4bc, "read counts.\n"))
rna <- rna[, !(colnames(rna) %in% names(bcsum0))]


# Select cell barcodes detected by both RNA and HTO in the datasets 
joint.bcs <- intersect(colnames(rna), colnames(hto))
# Subset RNA and HTO counts by joint cell barcodes
rna <- rna[, joint.bcs]
hto <- as.matrix(hto[, joint.bcs])
cat(paste("Keep", length(joint.bcs), "cell barcodes detected by both RNA and HTO in the datasets\n\n"))


# Setup Seurat object
sobj <- CreateSeuratObject(counts = rna,
                           project = "SeuratProject", # Project name for the ‘Seurat’ object
                           assay = "RNA", # Name of the initial assay
                           names.field = 1, # For the initial identity class for each cell, choose this field from the cell's name. E.g. If your cells are named as BARCODE_CLUSTER_CELLTYPE in the input matrix, set 'names.field' to 3 to set the initial identities to CELLTYPE
                           names.delim = "_", # For the initial identity class for each cell, choose this delimiter from the cell's column name.
                           min.cells = 1, # Include features detected in at least this many cells.
                           min.features = 1, # Include cells where at least this many features are detected
                           row.names = NULL # When ‘counts’ is a ‘data.frame’ or ‘data.frame’-derived object: an optional vector of feature names to be used
)

# Normalize RNA data with log normalization
sobj <- NormalizeData(sobj)
# Find and scale variable features
sobj <- FindVariableFeatures(sobj, selection.method = "mean.var.plot")
sobj <- ScaleData(sobj, features = VariableFeatures(sobj))


## Adding HTO data as an independent assay
#You can read more about working with multi-modal data here: https://satijalab.org/seurat/articles/multimodal_vignette.html
# Add HTO data as a new assay independent from RNA
sobj[["HTO"]] <- CreateAssayObject(counts = hto)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
sobj <- NormalizeData(sobj, assay = "HTO", normalization.method = "CLR")


## Demultiplex cells based on HTO enrichment
# Here we use the Seurat function HTODemux() to assign single cells back to their sample origins.
# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
sobj <- HTODemux(sobj, assay = "HTO", positive.quantile = 0.99)

## Visualize demultiplexing results
# Output from running HTODemux() is saved in the object metadata. We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.
print("Global classification results:")
print(table(sobj$HTO_classification.global))

# Visualize enrichment for selected HTOs with ridge plots
# Group cells based on the max HTO signal
Idents(sobj) <- "HTO_maxID"
plot1 <- RidgePlot(sobj, assay = "HTO", features = rownames(sobj[["HTO"]]), ncol = 2) 
ggsave(filename=file.path(out, "RidgePlot_HTOs.pdf"), plot=plot1, device = "pdf", dpi=300,
       units="cm", width=2*15, height = 15*ceiling(length(rownames(sobj[["HTO"]]))/2))

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
plot2 <- FeatureScatter(sobj, feature1 = rownames(sobj[["HTO"]])[1], feature2 = rownames(sobj[["HTO"]])[2])
ggsave(filename=file.path(out, "scatter_HTOs.pdf"), plot=plot2, device = "pdf", dpi=300,
       units="cm", width=15, height = 15)

# Compare number of UMIs for singlets, doublets and negative cells
Idents(sobj) <- "HTO_classification.global"
plot3 <- VlnPlot(sobj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
ggsave(filename=file.path(out, "violinPlot_doublets.pdf"), plot=plot3, device = "pdf", dpi=300,
       units="cm", width=15, height = 15)


# Generate a two dimensional tSNE embedding for HTOs. Here we are grouping cells by singlets and doublets for simplicity.
# First, we will remove negative cells from the object
sobj.subset <- subset(sobj, idents = "Negative", invert = TRUE)
# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = sobj.subset, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix
sobj.subset <- RunTSNE(sobj.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
plot4 <- DimPlot(sobj.subset)
ggsave(filename=file.path(out, "tSNE_doublets.pdf"), plot=plot4, device = "pdf", dpi=300,
       units="cm", width=15, height = 15)


# You can also visualize the more detailed classification result by running Idents(object) <-
# 'HTO_classification' beofre plotting. Here, you can see that each of the small clouds on the
# tSNE plot corresponds to one of the 28 possible doublet combinations.

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.
# To increase the efficiency of plotting, you can subsample cells using the ncells argument (default 5000). 
# This avoids having to draw exceptionally large heatmaps.
plot5 <- HTOHeatmap(sobj, assay = "HTO", ncells = 5000)
ggsave(filename=file.path(out, "heatmap_HTOs.pdf"), plot=plot5, device = "pdf", dpi=300,
       units="cm", width=15, height = 15)


# Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.
# Extract the singlets
sobj.singlet <- subset(sobj, idents = "Singlet")
# Select the top 1000 most variable features
sobj.singlet <- FindVariableFeatures(sobj.singlet, selection.method = "mean.var.plot")
# Scaling RNA data, we only scale the variable features here for efficiency
sobj.singlet <- ScaleData(sobj.singlet, features = VariableFeatures(sobj.singlet))
# Run PCA
sobj.singlet <- RunPCA(sobj.singlet, features = VariableFeatures(sobj.singlet))
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
sobj.singlet <- FindNeighbors(sobj.singlet, reduction = "pca", dims = 1:10)
sobj.singlet <- FindClusters(sobj.singlet, resolution = 0.6, verbose = FALSE)
sobj.singlet <- RunTSNE(sobj.singlet, reduction = "pca", dims = 1:10)
# Projecting singlet identities on TSNE visualization
plot6 <- DimPlot(sobj.singlet, group.by = "HTO_classification")
ggsave(filename=file.path(out, "tSNE_singlets.pdf"), plot=plot6, device = "pdf", dpi=300,
       units="cm", width=15, height = 15)

meta <- sobj.singlet@meta.data
write.table(data.frame(cell_id=rownames(meta), meta), file=file.path(out, "demux.txt"), sep="\t", quote = F, row.names = F)

