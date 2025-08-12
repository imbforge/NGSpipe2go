#####################################
##
## What: norm_bioc.R
## Who : Frank Rühle, Patrick Hüther
## When: 11.06.2025
##
## Normalize gene expression data in singlecellexperiment object and calculate HVGs and reduced dims.
## This module stores the sce object as HDF5 object. From now on, only metadata will be modified.
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## spikein_norm     # use spike-in genes for normalization
## maxgenes2plot    # max number of genes to plot (top expressed genes and HVGs)
## org              # either "human" or "mouse". Organism name needed for cell cycle phase scores. Leave empty to skip.
## explanatory_vars # any other categorical factors you may inspect for potential batch effects
## hvg_prop         # proportion of genes to report as HVGs according to biological component of variance
## block_var        # name of block var to account for uninteresting factor when decomposing variance for HVGs. Leave empty to skip.
## perplexity       # perplexity for tSNE. Can have a large effect on the results. calculateTSNE: Should not be bigger than (nrow(X)-1)/3. Reasonable default value: cellcount/5, capped at a maximum of 50.
## n_neighbors      # number of nearest neighbors for UMAP
## annocat_plot     # category used for plotting
## annocat_plot2    # 2nd category used for plotting
## plot_pointsize   # dot size used in plots
## plot_pointalpha  # dot transparency used in plots
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

as.char.vector <- function(x) { # allows input "c(A, B)" as well as "A, B". Avoids need to escape quotation marks in module header.
  char.vec <- gsub("(^c\\()|(\\)$)", "", x)  
  trimws(strsplit(char.vec, ",")[[1]])
}

args <- commandArgs(T)
resultsdir       <- parseArgs(args,"res=")   
seqtype          <- parseArgs(args,"seqtype=")   
outdir           <- parseArgs(args,"outdir=") # output folder
pipeline_root    <- parseArgs(args,"pipeline_root=") 
spikein_norm     <- parseArgs(args,"spikein_norm=", convert="as.logical", default = FALSE)
maxgenes2plot    <- parseArgs(args,"maxgenes2plot=", convert="as.numeric", default = 15)
org              <- parseArgs(args,"org=")
explanatory_vars <- parseArgs(args,"explanatory_vars=", convert="as.char.vector")
hvg_prop         <- parseArgs(args,"hvg_prop=", convert="as.numeric", default = 0.1)
block_var        <- parseArgs(args,"block_var=")
perplexity       <- parseArgs(args,"perplexity=", convert="run_custom_code") 
n_neighbors      <- parseArgs(args,"n_neighbors=", convert="run_custom_code")
  if(is.na(perplexity)) {perplexity <- NULL}
  if(is.na(n_neighbors)) {n_neighbors <- NULL}
annocat_plot     <- parseArgs(args,"annocat_plot=", default = "group")
annocat_plot2    <- parseArgs(args,"annocat_plot2=", default = "group")
plot_pointsize   <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha  <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  

# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_norm", "bioc_3.16.lock")
print(paste("load renv:", env.path))
renv::use(lockfile=env.path)

library(ggplot2)

# set options
addTaskCallback(function(...) {set.seed(100);TRUE})
options(stringsAsFactors=FALSE)

# check parameter
print(paste("resultsdir:", resultsdir))
print(paste("seqtype:", seqtype))
print(paste("outdir:", outdir))
print(paste("pipeline_root:", pipeline_root))
print(paste("spikein_norm:", spikein_norm))
print(paste("maxgenes2plot:", maxgenes2plot))
print(paste("org:", org))
print(paste("explanatory_vars:", paste0(explanatory_vars, collapse=", ")))
print(paste("hvg_prop:", hvg_prop))
print(paste("block_var:", block_var))
print(paste("perplexity:", paste0(perplexity, collapse=", ")))
print(paste("n_neighbors:", paste0(n_neighbors, collapse=", ")))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- readr::read_rds(file.path(resultsdir, "sce.RDS"))
print(sce)


# Normalization
if(spikein_norm) {
  # use spike-ins for normalization
  print("Summary size factors for ERCC spike-ins (must be >0)")
  sce <- scuttle::computeSpikeFactors(sce, spikes="spikein")
  print(summary(SingleCellExperiment::sizeFactors(sce)))
  } else {
    # deconvolve size factors
    print("deconvolve size factors (must be >0)")
    sce <- scran::computeSumFactors(sce,clusters=scran::quickCluster(sce))
    print(summary(SingleCellExperiment::sizeFactors(sce)))
  }

# calculate per-feature average counts on filtered cells
SummarizedExperiment::rowData(sce)$ave_count <- scuttle::calculateAverage(sce) # size factor-adjusted average count
sce <- scuttle::addPerFeatureQCMetrics(sce) # gives mean counts and percentage of observations above threshold (0).

# add normalized logcounts matrix
sce <- scater::logNormCounts(sce)

# plot size factors
print("plot size factors")
plot_sf <- ggplot(SummarizedExperiment::colData(sce), aes(x=sizeFactor, y=sum/1e3, color=!!sym(annocat_plot))) +
  geom_point(size = plot_pointsize, alpha=plot_pointalpha) + 
  ylab("Library size (thousands)") + 
  scale_color_hue(l=55) +
  guides(color=guide_legend(override.aes = list(size=1, alpha=1)))

ggsave(plot=plot_sf, filename= file.path(outdir, "size_factors.png"), device="png", width=7, height=4, bg = "white")
ggsave(plot=plot_sf, filename= file.path(outdir, "size_factors.pdf"), device="pdf", width=7, height=4)


## plot genes with highest expression (gene filter does not effect order of colData)
# determine ideal number of panel columns in violin plots when using facet_wrap to avoid too many violin plots per row.
print("plot genes with highest expression")
maxViolPerRow = max(nlevels(factor(SummarizedExperiment::colData(sce)[,annocat_plot])), 10) # i.e. max #genes x #groups. 
nGridCol <- max(1,floor(maxViolPerRow/nlevels(factor(SummarizedExperiment::colData(sce)[,annocat_plot])))) # number of cols in violin plot
nGridRow <- ceiling(maxgenes2plot / nGridCol) # number of rows in violin plot

top_expressed_genes <- SummarizedExperiment::rowData(sce) |>
  tibble::as_tibble() |>
  dplyr::arrange(dplyr::desc(ave_count)) |>
  dplyr::slice(1:maxgenes2plot)

top_expressed_data <- SingleCellExperiment::logcounts(sce)[top_expressed_genes$feature_id,] |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
  dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot]) |>
  dplyr::mutate(!!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2]) |>
  tidyr::pivot_longer(cols=!any_of(c(annocat_plot, annocat_plot2))) |>
  dplyr::mutate(symbol=SummarizedExperiment::rowData(sce)$feature_symbol[match(name, rownames(sce))])

top_expressed_plot <- ggplot(top_expressed_data, aes(x=!!dplyr::sym(annocat_plot),y=value)) +
  geom_violin() +
  ggbeeswarm::geom_quasirandom(aes(color=!!dplyr::sym(annocat_plot2)), size = plot_pointsize, alpha=plot_pointalpha) +
  scale_color_hue(l=55) +
  geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA) +
  facet_wrap(as.formula(paste("~", "symbol")), ncol=nGridCol) +
  ylab("expression (logcounts)") +
  xlab(element_blank()) +
  ggtitle(paste0("Top ", maxgenes2plot, " expressed genes")) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), legend.position = "bottom") + 
  guides(color=guide_legend(override.aes = list(size=1, alpha=1)))


ggsave(plot=top_expressed_plot, width=7, height=1+2*nGridRow, # legend + 2 inch per panel row
       filename= file.path(outdir, paste0("top_", maxgenes2plot, "_expressed_genes.png")), device="png", bg = "white")
ggsave(plot=top_expressed_plot, width=7, height=1+2*nGridRow, # legend + 2 inch per panel row
       filename= file.path(outdir, paste0("top_", maxgenes2plot, "_expressed_genes.pdf")), device="pdf")



## Estimate highly variable genes
# Model the variance of the log-expression profiles for each gene, decomposing it into technical and biological components based on a 
# fitted mean-variance trend. Use of block is the recommended approach for accounting for any uninteresting categorical factor of variation.
if(spikein_norm) {
  print(paste("Estimate highly variable genes to be used for reduced dimensions based on spike-ins", if(!is.na(block_var)) {paste("separately within each level of", block_var)}))
  decVar <- scran::modelGeneVarWithSpikes(sce, spikes = "spikein", block = if(!is.na(block_var)) {SummarizedExperiment::colData(sce)[,block_var]} else {NULL})
} else {
  print(paste("Estimate highly variable genes to be used for reduced dimensions", if(!is.na(block_var)) {paste("separately within each level of", block_var)}))
  decVar <- scran::modelGeneVar(sce, assay.type = "logcounts", block = if(!is.na(block_var)) {SummarizedExperiment::colData(sce)[,block_var]} else {NULL})
}
hvg <- scran::getTopHVGs(decVar, var.field = "bio", prop=hvg_prop) # already ordered
SingleCellExperiment::rowSubset(sce, field = "HVGs") <- hvg # add 'HVGs' column (logical) to rowData
print(paste("Determined", length(hvg), "HVGs"))

decVar_ordered <- as.data.frame(decVar) |> # Ordering by most interesting genes for inspection
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(symbol = SummarizedExperiment::rowData(sce)[gene_id, "feature_symbol"]) |>
  dplyr::filter(gene_id %in% hvg) |>
  dplyr::relocate(gene_id, symbol) |>
  dplyr::arrange(dplyr::desc(bio)) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, digits = 4)))

write.table(decVar_ordered, file=file.path(outdir, paste0("highly_variable_genes_bio", hvg_prop, ".txt")), sep="\t", quote = F, row.names = F)


# plot top HVGs (gene filter does not effect order of colData)
# plot genes with highest expression
top_hvg_data <- SingleCellExperiment::logcounts(sce)[hvg[1:maxgenes2plot],] |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
  dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot]) |>
  dplyr::mutate(!!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2]) |>
  tidyr::pivot_longer(cols=!any_of(c(annocat_plot, annocat_plot2))) |>
  dplyr::mutate(symbol=SummarizedExperiment::rowData(sce)$feature_symbol[match(name, rownames(sce))])

top_hvg_plot <- ggplot(top_hvg_data, aes(x=!!dplyr::sym(annocat_plot),y=value)) +
  geom_violin() +
  ggbeeswarm::geom_quasirandom(aes(color=!!dplyr::sym(annocat_plot2)), size = plot_pointsize, alpha=plot_pointalpha) +
  scale_color_hue(l=55) +
  geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA) +
  facet_wrap(as.formula(paste("~", "symbol")), ncol=nGridCol) +
  ylab("expression (logcounts)") +
  xlab(element_blank()) +
  ggtitle(paste0("Top ", maxgenes2plot, " highly variable genes")) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), legend.position = "bottom") + 
  guides(color=guide_legend(override.aes = list(size=1, alpha=1)))


ggsave(plot=top_hvg_plot, width=7, height=1+2*nGridRow, # legend + 2 inch per panel row
       filename= file.path(outdir, paste0("top_", maxgenes2plot, "_highly_variable_genes.png")), device="png", bg = "white")
ggsave(plot=top_hvg_plot, width=7, height=1+2*nGridRow, # legend + 2 inch per panel row
       filename= file.path(outdir, paste0("top_", maxgenes2plot, "_highly_variable_genes.pdf")), device="pdf")


## Check for potential confounder variables

# Determine cell cycle phase
if(!is.null(org) && org %in% c("human", "mouse")) {
  print("Determine Cell cycle phases")
  # load gene.pairs: "mouse_cycle_markers.rds" or "human_cycle_markers.rds"
  gene.pairs <- readRDS(system.file("exdata", paste0(org, "_cycle_markers.rds"), package="scran")) 
  
  # in "human_cycle_markers.rds" ensembl names are given without ".versionnumber" at the end,  
  # therefore, we remove them from the sce rownames if present.
  assignments <- scran::cyclone(sce, pairs=gene.pairs, gene.names=sub("\\..*$","",SummarizedExperiment::rowData(sce)$feature_id ), assay.type="logcounts") 
  write.table(assignments, file=file.path(outdir, "cell_cycle_phase_assignments.txt"), sep="\t", quote=F, row.names=F)
  
  sce$cc_phases <- assignments$phases
  cat(sum(is.na(sce$cc_phases)), "out of", length(sce$cc_phases),
      "cells could not be assigned to a cell cycle phase.", fill=TRUE)
  
  if(any(!is.na(sce$cc_phases))) {
    
    # print cell cycle summary table
    ccptable <- SummarizedExperiment::colData(sce) |>
      tibble::as_tibble() |>
      dplyr::group_by(!!dplyr::sym(annocat_plot)) |>
      dplyr::summarize(G1 = paste0(sum(cc_phases=="G1"), " (", sprintf("%1.1f%%", 100*sum(cc_phases=="G1")/length(cc_phases)), ")"),
                       G2M = paste0(sum(cc_phases=="G2M"), " (", sprintf("%1.1f%%", 100*sum(cc_phases=="G2M")/length(cc_phases)), ")"),
                       S = paste0(sum(cc_phases=="S"), " (", sprintf("%1.1f%%", 100*sum(cc_phases=="S")/length(cc_phases)), ")")
      ) |>
      dplyr::ungroup() |>
      tibble::column_to_rownames(annocat_plot) |>
      t() |>
      as.data.frame()
    
      ccptable$total <- c(paste0(sum(sce$cc_phases=="G1"), " (", sprintf("%1.1f%%", 100*sum(sce$cc_phases=="G1")/ncol(sce)), ")"),
                          paste0(sum(sce$cc_phases=="G2M"), " (", sprintf("%1.1f%%", 100*sum(sce$cc_phases=="G2M")/ncol(sce)), ")"),
                          paste0(sum(sce$cc_phases=="S"), " (", sprintf("%1.1f%%", 100*sum(sce$cc_phases=="S")/ncol(sce)), ")"))
    
    write.table(data.frame(phase=rownames(ccptable), ccptable), file=file.path(outdir, "cell_cycle_phase_summary_table.txt"), sep="\t", quote=F, row.names = F)
  
    # plot cell cycle phases
    ccpscores <- assignments$scores |>
      dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot]) 
      
    ccp_plot <- ggplot(ccpscores, aes(x=G1, y=G2M, color=!!sym(annocat_plot))) +
      geom_point(size = plot_pointsize, alpha=plot_pointalpha) + 
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      xlab("G1 score") + 
      ylab("G2/M score") + 
      scale_color_hue(l=55) +
      guides(color=guide_legend(override.aes = list(size=1, alpha=1))) +
      annotate(geom="segment", x = 0.5, xend = 0.5, y = 0,   yend = 0.5, linetype = "dashed", color = "black") +
      annotate(geom="segment", x = 0,   xend = 0.5, y = 0.5, yend = 0.5, linetype = "dashed", color = "black") +
      annotate(geom="segment", x = 0.5, xend = 1, y = 0.5,   yend = 1,   linetype = "dashed", color = "black") +
      annotate(geom="text", x = 0.75, y = 0.25, label = "G1") +
      annotate(geom="text", x = 0.25, y = 0.75, label = "G2") +
      annotate(geom="text", x = 0.25, y = 0.25, label = "S") + 
      ggtitle("Cell Cycle Phase Scores") 
      
    ggsave(plot=ccp_plot, filename= file.path(outdir, "cell_cycle_phase_score_plot.png"), device="png", width=7, height=5, bg = "white")
    ggsave(plot=ccp_plot, filename= file.path(outdir, "cell_cycle_phase_score_plot.pdf"), device="pdf", width=7, height=5)
   
  }
}

# inspect other potential confounder variables
print("inspect potential confounder variables")
explanatory_vars <- intersect(c(annocat_plot, annocat_plot2, explanatory_vars, "cc_phases"), colnames(SummarizedExperiment::colData(sce)))

explVar <- scater::plotExplanatoryVariables(sce, variables=explanatory_vars, exprs_values = "logcounts") +    
  scale_color_hue(l=55, name = element_blank()) +
  theme(legend.position = "bottom" ) +
  ggtitle("Explanatory variables")
ggsave(plot=explVar, filename= file.path(outdir, paste0("explanatory_variables_density_plot.png")), device="png", width=7, height=5, bg = "white")
ggsave(plot=explVar, filename= file.path(outdir, paste0("explanatory_variables_density_plot.pdf")), device="pdf", width=7, height=5)

# calculate PCA based on hvg for scatter plots (could also be done for spike-ins with altexp = "spikein")
sce <- scater::runPCA(sce, name = "PCA", ncomponents=50, subset_row=hvg, exprs_values="logcounts")

explVarScatter <- lapply(explanatory_vars, function(x) {
  p <- scater::plotPCASCE(sce, ncomponents = 4, point_size=plot_pointsize, point_alpha=plot_pointalpha, colour_by=x) + 
    scale_color_hue(l=55) +
    ggtitle(paste("PCA plots colored by", x)) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha=1)))  
    
  ggsave(plot=p, filename= file.path(outdir, paste0("explanatory_variables_scatter_plot_", x, ".png")), device="png", width=7, height=8, bg = "white")
  ggsave(plot=p, filename= file.path(outdir, paste0("explanatory_variables_scatter_plot_", x, ".pdf")), device="pdf", width=7, height=8)
  return(p)
})


## calculate reduced dimensions TSNE and UMAP (PCA done above)
print("calculate reduced dimensions")
for (p in perplexity) {
  sce <- scater::runTSNE(sce, name = if(length(perplexity)==1) {"TSNE"} else {paste0("TSNE_p", p)}, ncomponents=2, subset_row=hvg,
                         perplexity=p, dimred = NULL, exprs_values="logcounts") # ncomponents for TSNE must be either 1, 2 or 3
  }
for (n in n_neighbors) {
  sce <- scater::runUMAP(sce, name = if(length(n_neighbors)==1) {"UMAP"} else {paste0("UMAP_n", n)}, ncomponents=2, subset_row=hvg,
                         n_neighbors = n, dimred = NULL, exprs_values="logcounts")
}

# reduced dimension plots
for (i in SingleCellExperiment::reducedDimNames(sce)) {
  
  plot_data <- SingleCellExperiment::reducedDim(sce, i) |>
    as.data.frame() |>
    dplyr::select(1:2) |> 
    stats::setNames(c(paste0(gsub("_.*$", "", i), " 1"), paste0(gsub("_.*$", "", i), " 2"))) |>
    dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot],
                  !!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2])
  
  rdplot_base <- ggplot(plot_data, aes(x = !!sym(names(plot_data)[1]), y = !!sym(names(plot_data)[2]))) +
    scale_color_hue(l=55) +
    ggtitle(paste("Plot ", i)) + 
    theme(legend.position = "bottom") + 
    guides(color=guide_legend(override.aes = list(size=1, alpha=1)))
  
  rdplot <- rdplot_base +
    geom_point(aes(color = !!sym(annocat_plot)), size = plot_pointsize, alpha=plot_pointalpha) 
  
  rdplot_anno1 <- rdplot + 
    facet_wrap(as.formula(paste("~", annocat_plot2)))
  
  rdplot_anno2 <- rdplot_base + 
    geom_point(aes(color = !!sym(annocat_plot2)), size = plot_pointsize, alpha=plot_pointalpha) +
    facet_wrap(as.formula(paste("~", annocat_plot)))
  
  ggsave(plot=rdplot, filename= file.path(outdir, paste0("dimred_plot_", i, ".png")), device="png", width=7, height=8, bg = "white")
  ggsave(plot=rdplot, filename= file.path(outdir, paste0("dimred_plot_", i, ".pdf")), device="pdf", width=7, height=8)
  
  plotlayout <- ggplot2::ggplot_build(rdplot_anno1)$layout$layout
  ggsave(plot=rdplot_anno1, filename= file.path(outdir, paste0("dimred_plot_", i, "_split_by_", annocat_plot2, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
  ggsave(plot=rdplot_anno1, filename= file.path(outdir, paste0("dimred_plot_", i, "_split_by_", annocat_plot2, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
  
  plotlayout <- ggplot2::ggplot_build(rdplot_anno2)$layout$layout
  ggsave(plot=rdplot_anno2, filename= file.path(outdir, paste0("dimred_plot_", i, "_split_by_", annocat_plot, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
  ggsave(plot=rdplot_anno2, filename= file.path(outdir, paste0("dimred_plot_", i, "_split_by_", annocat_plot, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))

}


# doublet score plots illustrated by all reduced dimensions and incl plots split by annotation categories
print("create doublet score plots")
purrr::map(SingleCellExperiment::reducedDimNames(sce), function(dimred) {
  
  plot_data <- SingleCellExperiment::reducedDim(sce, dimred) |>
    as.data.frame() |>
    dplyr::select(1:2) |> 
    stats::setNames(c(paste0(gsub("_.*$", "", dimred), " 1"), paste0(gsub("_.*$", "", dimred), " 2"))) |>
    dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot],
                  !!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2],
                  doublet_score = SummarizedExperiment::colData(sce)[,"scDblFinder.score"])
  
  rdplot <- ggplot(plot_data, aes(x = !!dplyr::sym(names(plot_data)[1]), y = !!dplyr::sym(names(plot_data)[2]), color = doublet_score)) +
    geom_point(size = plot_pointsize, alpha=plot_pointalpha) +
    viridis::scale_color_viridis() + 
    ggtitle(paste("scDblFinder doublet score")) + 
    theme(legend.position = "bottom") 
  
  rdplot_anno1 <- rdplot + 
    facet_wrap(as.formula(paste("~", annocat_plot2)))
  
  rdplot_anno2 <- rdplot + 
    facet_wrap(as.formula(paste("~", annocat_plot)))
  
  ggsave(plot=rdplot, filename= file.path(outdir, paste0("doublet_score_by_", dimred, ".png")), device="png", width=7, height=8, bg = "white")
  ggsave(plot=rdplot, filename= file.path(outdir, paste0("doublet_score_by_", dimred, ".pdf")), device="pdf", width=7, height=8)
  
  plotlayout <- ggplot2::ggplot_build(rdplot_anno1)$layout$layout
  ggsave(plot=rdplot_anno1, filename= file.path(outdir, paste0("doublet_score_by_", dimred, "_split_by_", annocat_plot2, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
  ggsave(plot=rdplot_anno1, filename= file.path(outdir, paste0("doublet_score_by_", dimred, "_split_by_", annocat_plot2, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
  
  plotlayout <- ggplot2::ggplot_build(rdplot_anno2)$layout$layout
  ggsave(plot=rdplot_anno2, filename= file.path(outdir, paste0("doublet_score_by_", dimred, "_split_by_", annocat_plot, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
  ggsave(plot=rdplot_anno2, filename= file.path(outdir, paste0("doublet_score_by_", dimred, "_split_by_", annocat_plot, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
})




#############################
# save the sessionInformation and R image
print("store results")
writeLines(capture.output(sessionInfo()),paste0(outdir, "/norm_bioc_session_info.txt"))
readr::write_rds(sce, file = file.path(resultsdir, "sce.RDS"))
HDF5Array::saveHDF5SummarizedExperiment(sce, dir=file.path(resultsdir, "HDF5"), prefix="sce_")
save(hvg, top_hvg_data, file=paste0(outdir,"/norm_bioc.RData"))

