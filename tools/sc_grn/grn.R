#####################################
##
## What: grn.R
## Who : Frank Rühle
## When: 06.06.2023
##
## Script to identify gene regulatory networks with Pando.
##
## Args:
## -----
## projectdir      # project directory
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
run_custom_code <- function(x) {
  eval(parse(text=x))
}

args <- commandArgs(T)
projectdir    <- parseArgs(args,"project=") 
resultsdir    <- parseArgs(args,"res=")   
out           <- parseArgs(args,"outdir=") # output folder
peak_assay    <- parseArgs(args,"peak_assay=")     
rna_assay     <- parseArgs(args,"rna_assay=")     
db            <- parseArgs(args,"db=")     
methodModel   <- parseArgs(args,"methodModel=")     
genes2use     <- parseArgs(args,"genes2use=", convert="as.character")
pval_thresh   <- parseArgs(args,"pval_thresh=", convert="as.numeric")
min_genes     <- parseArgs(args,"min_genes=", convert="as.numeric")
features4graph <- parseArgs(args,"features4graph=", convert="as.character")
umap_method   <- parseArgs(args,"umap_method=")     
n_neighbors   <- parseArgs(args,"n_neighbors=", convert="as.numeric")   

runstr <- "Rscript grn.R [projectdir=projectdir]"

if(length(genes2use)==1) {
  if(is.na(genes2use)) {genes2use <- NULL} else {
    if(file.exists(genes2use)) {genes2use <- read.delim(genes2use, header=F)[,1]} else {
      cat("\npath in genes2use not found. genes2use is set to NULL\n")
      genes2use <- NULL
    }
  }
}
if(length(features4graph)==1) {
  if(is.na(features4graph)) {features4graph <- NULL} else {
    if(file.exists(features4graph)) {features4graph <- read.delim(features4graph, header=F)[,1]} else {
      cat("\npath in features4graph not found. features4graph is set to NULL\n")
      features4graph <- NULL
    }
  }
}

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_grn/renv.lock"))
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
library(openxlsx)
library(motifmatchr)
library(JASPAR2020)
library(Pando)
library(grr)

# set options
options(stringsAsFactors=FALSE)
addTaskCallback(function(...) {set.seed(100);TRUE})

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("peak_assay:", peak_assay))
print(paste("rna_assay:", rna_assay))
print(paste("db:", db))
print(paste("methodModel:", methodModel))
print(paste("genes2use:", paste(genes2use, collapse=" ")))
print(paste("pval_thresh:", pval_thresh))
print(paste("min_genes:", min_genes))
print(paste("features4graph:", paste(features4graph, collapse=" ")))
print(paste("umap_method:", umap_method))
print(paste("n_neighbors:", n_neighbors))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- peak_assay


# load relevant BSgenome package (needed by Signac for motif analysis)
switch(db,
       hg38={ failed_BSgenome <- library("BSgenome.Hsapiens.UCSC.hg38") 
       BSgenome <- BSgenome.Hsapiens.UCSC.hg38
       JasparTaxGroup <- "vertebrates"
       },
       mm10={ failed_BSgenome <- library("BSgenome.Mmusculus.UCSC.mm10") 
       BSgenome <- BSgenome.Mmusculus.UCSC.mm10
       JasparTaxGroup <- "vertebrates"
       },
       stop(c("Don't find genome:", db))   
)


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = JasparTaxGroup, all_versions = FALSE)
)


# Select variable features
sobj <- Seurat::FindVariableFeatures(sobj, assay=rna_assay)

# Initiate GRN object and select candidate regions. 
# This will create a RegulatoryNetwork object inside the Seurat object 
# and select candidate regulatory regions. By default, Pando will consider 
# all peaks as putative regulatory regions, but the set of candidate 
# regions can be constrained by providing a GenomicRanges object in the 
# regions argument. Pando ships with a set of conserved regions 
# (phastConsElements20Mammals.UCSC.hg38) as well as predicted regulatory 
# elements from ENCODE (SCREEN.ccRE.UCSC.hg38) for the human genome (hg38).
# data('phastConsElements20Mammals.UCSC.hg38')
# data('SCREEN.ccRE.UCSC.hg38') # e.g. use it in initiate_grn as:
# regions = union(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38)
# When working with genomes other than human you can use the peaks that are correlated 
# with the expression of nearby genes identified by the Seurat function LinkPeaks (see above).
# Using only these genes is quite strict, as only a small fraction of genes will be linked, 
# but the developers found that this does result in quite robust GRNs. Use as:
# regions = StringToGRanges(Links(sobj[["ATAC"]])$peak)
sobj <- initiate_grn(sobj,
                     regions = NULL, # if specified, intersections with peak ranges from peak_assay are used
                     peak_assay = peak_assay,
                     rna_assay = rna_assay,
                     exclude_exons = TRUE) # new object is now called SeuratPlus
# GetGRN(sobj) # get RegulatoryNetwork object 
# regions <- NetworkRegions(sobj) # We can also inspect the candidate regulatory regions that we have selected:
# regions@ranges
# regions@peaks

# Scan candidate regions for TF binding motifs
# This uses motifmatchr to pair up TFs with their putative binding sites. 
# Pando provides a custom motif database (motifs) compiled from JASPAR and CIS-BP, 
# but in principle any PFMatrixList object can be provided here. 
# A data frame with motif-to-TF assignments (motif IDs in 1st column and TF names in 2nd column) 
# can be provided in the motif_tfs argument.
# data(motifs)
# data('motif2tf')
sobj <- find_motifs(
  sobj,
  pfm = pfm,
  genome = BSgenome,
  motif_tfs = NULL
)
# The Regions object has now gotten a new slot containing the Motif object. 
# This object stores a sparse peak x motif matrix with the matches and some other information.
# regions <- NetworkRegions(sobj)
# regions@motifs@data[1:5,1:5]

# Infer gene regulatory network by fitting regression models for the expression of each gene. 
# Here, we first select regions near genes, either by simply considering a distance 
# upstream and/or downstream of the gene (peak_to_gene_method='Signac') or by also 
# considering overlapping regulatory regions as is done by GREAT (peak_to_gene_method='GREAT')
sobj <- infer_grn(sobj,
                  peak_to_gene_method = 'Signac', # One of 'Signac' or 'GREAT'
                  method = methodModel,
                  genes = genes2use, # select a subset of genes that we want to use for GRN inference
                  upstream = 1e+05,
                  downstream = 0,
                  tf_cor = 0.1,
                  peak_cor = 0) 
# GetNetwork(sobj) # access the inferred Network object

# Once the models are fit, model coefficients can be inspected with coef().
# This returns a dataframe coefficients (estimate) and p-values (pval) for each TF/region-target gene pair.
# coef(sobj)

# Module extraction (i.e. the set of genes that are regulated by each transcription factor).
# Based on the model coefficients, we can construct a network between TFs and target genes. 
# This can be further summarized to construct gene and regulatory modules with the set of 
# target genes and regulatory regions for each TF. 
sobj <- find_modules(sobj, p_thresh = pval_thresh,
                     rsq_thresh = 0.1,
                     nvar_thresh = 10,
                     min_genes_per_module = min_genes)

# Print modules
modules <- NetworkModules(sobj) # The meta slot holds a dataframe with module information.
modules_df <- as.data.frame(modules@meta)
write.table(modules_df, file =file.path(out, "GRN_TF_modules.txt"), quote = F, row.names = F,  sep="\t")

cat("Plot goodness-of-fit metrics of the fitted models. This shows us the explained variance of the model, as well as the number of variables in each model. We can observe that with more variables, the models will have a better fit. The dashed lines show the applied thresholds for in the find_modules function.\n")
pgof <- plot_gof(sobj, point_size=3) 
ggsave(plot=pgof, filename=file.path(out, "GRN_goodness_of_fit_metrics.pdf"))
ggsave(plot=pgof, filename=file.path(out, "GRN_goodness_of_fit_metrics.png"))

cat("Plot size of the modules with respect to the number of target genes as well as the number of peaks and regulating TFs per target gene.\n")
pmm <- plot_module_metrics(sobj) 
ggsave(plot=pmm, filename=file.path(out, "GRN_module_metrics.pdf"))
ggsave(plot=pmm, filename=file.path(out, "GRN_module_metrics.png"))

# create the graph to be visualized and optionally a UMAP embedding for the nodes 
sobj <- get_network_graph(sobj, 
                          rna_assay = rna_assay, 
                          graph_name = "module_graph", # we can access the graph with this name later
                          umap_method = umap_method, # method to compute edge weights for UMAP
                          n_neighbors = n_neighbors, # The size of local neighborhood (in terms of number of neighboring sample points) used for manifold
                          # approximation. Larger values result in more global views of the manifold, while smaller values 
                          # result in more local data being preserved. In general values should be in the range 2 to 100 (default 15).
                          features = features4graph) # plot a smaller subgraph of the GRN by selecting a set of genes as features
cat("Visualization of the network graph. The nodes are colored and sized based on their centrality in the graph and the edges are colored by the direction of the regulation (inhibitory: grey; activating: orange)")
pnetg <- plot_network_graph(sobj, 
                            graph = "module_graph",
                            layout = "umap") + ggtitle("Network graph")
ggsave(plot=pnetg, filename=file.path(out, "GRN_network_graph.pdf"))
ggsave(plot=pnetg, filename=file.path(out, "GRN_network_graph.png"))


# # Optional: Sometimes it is useful to plot a subgraph of the GRN centered around one TF (e.g. NEUROD4), to reveal a hierarchy of regulation.
# sobj <- get_tf_network(sobj, tf="NEUROD4", graph="module_graph", keep_all_edges=F)
# plot_tf_network(sobj, tf='NEUROD4', circular=T)
# # further custom plotting of dedicated sub graphs if generated
# library(ggraph)
# subgraph <- NetworkGraph(sobj, graph='module_graph')
#  ggraph(subgraph, layout='drl') +
#     geom_edge_hive(aes(color=estimate), size=2, fill='grey') +
#     geom_node_point(aes(size=centrality), shape=21, fill='grey') +
#     geom_node_text(aes(label=name), size=4, shape=21, repel=T) +
#     theme_void()


#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/grn_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj_grn.RDS"))
save(modules_df, pgof, pmm, pnetg, file=paste0(out,"/grn.RData"))

