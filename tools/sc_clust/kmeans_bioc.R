#####################################
##
## What: kmeans_bioc.R
## Who : Frank Rühle, Patrick Hüther
## When: 13.08.2025
##
## Apply k-means clustering on singlecellexperiment object
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## data2clust_pre   # data used for clustering. Can be any reduced dimension slot or 'logcounts'.
## kmeans           # number of centers for k-means clustering
## n_neighbors      # number of nearest neighbors to consider during graph construction
## weights          # type of weighting scheme to use for shared neighbors
## algorithm        # clustering algorithm ('walktrap', 'louvain', 'infomap', 'fast_greedy', 'label_prop', 'leading_eigen')
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
data2clust_pre   <- parseArgs(args,"data2clust=", convert="as.char.vector") 
k                <- parseArgs(args,"kmeans=", convert="as.numeric")
annocat_plot     <- parseArgs(args,"annocat_plot=", default = "group")
annocat_plot2    <- parseArgs(args,"annocat_plot2=", default = "sample")
plot_pointsize   <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha  <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  

# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_clust", "bioc_3.16.lock")
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
print(paste("data2clust:", paste0(data2clust_pre, collapse=", ")))
print(paste("kmeans:", k))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- HDF5Array::loadHDF5SummarizedExperiment(dir=file.path(resultsdir, "HDF5"), prefix="sce_")
print(sce)


# set up combinations of clustering parameter (catch e.g. TSNE subtypes as well)
params <- tidyr::expand_grid( 
  data2clust = unique(data2clust_pre, unlist(sapply(data2clust_pre, grep, SingleCellExperiment::reducedDimNames(sce), value=T))), 
  k = k
)


if(any(is.na(params))) { # skip entirely if cluster setting not specified
  cl <- list()
} else {

  cl <- purrr::pmap(params, function(data2clust, k) {

      name_clustering <- paste0("cluster_kmeans_", data2clust, "_k", k)
      
      if(file.exists(file.path(outdir, name_clustering, paste0(name_clustering, ".tsv")))) {
        
        print(paste0(name_clustering, ".tsv already exists. Clustering for this setting is skipped"))
        cluster <- NULL
        
      } else {
        
        print(paste("Processing", name_clustering))
        if (!file.exists(file.path(outdir, name_clustering))) {dir.create(file.path(outdir, name_clustering), recursive=T) }
        
        cluster <- scran::clusterCells(sce[if(data2clust == "logcounts") SummarizedExperiment::rowData(sce)$HVGs else TRUE,],
                                       use.dimred = switch(data2clust != "logcounts",data2clust,NULL),
                                       assay.type = switch(data2clust == "logcounts",data2clust,NULL),
                                       BLUSPARAM=bluster::KmeansParam(centers=k))

        cluster_tsv <- data.frame(cell_id=colnames(sce), cluster=cluster) |>
          purrr::set_names("cell_id", name_clustering) |> 
          readr::write_tsv(file.path(outdir, name_clustering, paste0(name_clustering, ".tsv")))
        
        cells_per_cat <- data.frame(sample=SummarizedExperiment::colData(sce)[,annocat_plot2], cluster=cluster) |>
          table() |> 
          as.data.frame.matrix() |>
          tibble::rownames_to_column(var=annocat_plot2) |>
          readr::write_tsv(file.path(outdir, name_clustering, paste0(name_clustering, "_cellCounts_per_", annocat_plot2, ".txt")))
        
        # create cluster plots illustrated by all reduced dimensions and incl plots split by annotation categories
        purrr::map(SingleCellExperiment::reducedDimNames(sce), function(dimred) {
        
          plot_data <- SingleCellExperiment::reducedDim(sce, dimred) |>
            as.data.frame() |>
            dplyr::select(1:2) |> 
            stats::setNames(c(paste0(gsub("_.*$", "", dimred), " 1"), paste0(gsub("_.*$", "", dimred), " 2"))) |>
            dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot],
                          !!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2],
                          cluster = factor(cluster))
        
          rdplot <- ggplot(plot_data, aes(x = !!dplyr::sym(names(plot_data)[1]), y = !!dplyr::sym(names(plot_data)[2]), color = cluster)) +
            geom_point(size = plot_pointsize, alpha=plot_pointalpha) +
            scale_color_hue(l=55) +
            ggtitle(paste("Cluster kmeans", data2clust, paste0("k", k), "by", dimred)) + 
            theme(legend.position = "bottom") + 
            guides(color=guide_legend(override.aes = list(size=1, alpha=1)))
          
          rdplot_anno1 <- rdplot + 
            facet_wrap(as.formula(paste("~", annocat_plot2)))
          
          rdplot_anno2 <- rdplot + 
            facet_wrap(as.formula(paste("~", annocat_plot)))
          
          ggsave(plot=rdplot, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_by_", dimred, ".png")), device="png", width=7, height=8, bg = "white")
          ggsave(plot=rdplot, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_by_", dimred, ".pdf")), device="pdf", width=7, height=8)
          
          plotlayout <- ggplot2::ggplot_build(rdplot_anno1)$layout$layout
          ggsave(plot=rdplot_anno1, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_by_", dimred, "_split_by_", annocat_plot2, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
          ggsave(plot=rdplot_anno1, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_by_", dimred, "_split_by_", annocat_plot2, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
          
          plotlayout <- ggplot2::ggplot_build(rdplot_anno2)$layout$layout
          ggsave(plot=rdplot_anno2, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_by_", dimred, "_split_by_", annocat_plot, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
          ggsave(plot=rdplot_anno2, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_by_", dimred, "_split_by_", annocat_plot, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
        })
        
        
        # Approximate silhouette width for evaluating cluster separation
        if(length(unique(cluster))>1) {
          print(paste("Approximate silhouette width for", name_clustering))
          sil.approx <- bluster::approxSilhouette(x=if(data2clust != "logcounts") {SingleCellExperiment::reducedDim(sce, data2clust)} else {t(logcounts(sce)[SummarizedExperiment::rowData(sce)$HVGs,])}, 
                                                  clusters=cluster) |>
            as.data.frame() |>
            dplyr::mutate(closest=factor(ifelse(width > 0, cluster, other)))
          
          bplot <- ggplot(sil.approx, aes(x=cluster, y=width)) +
            geom_violin() +
            ggbeeswarm::geom_quasirandom(aes(colour=closest), size = plot_pointsize, alpha=plot_pointalpha) +
            scale_color_hue(l=55) +
            ggtitle(paste("Approx. silhouette width for", name_clustering)) + 
            theme(legend.position = "bottom") + 
            guides(color=guide_legend(override.aes = list(size=1, alpha=1)))
          
          ggsave(plot=bplot, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_approx_silhouette_width.png")), 
                 device="png", width=7, height=5, bg = "white") 
          ggsave(plot=bplot, filename= file.path(outdir, name_clustering, paste0(name_clustering, "_approx_silhouette_width.pdf")), 
                 device="pdf", width=7, height=5)
        } else {print(paste("skip approximate silhouette width for", name_clustering))}
        
        
        # Doublet detection by cluster (findDoubletClusters needs at least 3 clusters to detect doublet clusters)
        if(length(unique(cluster)) >=3) {
          print(paste("Cluster doublet detection for", name_clustering))
          dbl.out <- scDblFinder::findDoubletClusters(sce, clusters=cluster) |>
            as.data.frame() |> 
            tibble::rownames_to_column("cluster") |>
            dplyr::mutate(symbol=SummarizedExperiment::rowData(sce)[best, "feature_symbol"]) |>
            dplyr::mutate(dplyr::across(c(p.value, lib.size1, lib.size2, prop), \(x) signif(x, digits=3))) |>
            dplyr::mutate(problematic=scater::isOutlier(num.de, nmads=3, type="lower", log=TRUE)) |>
            dplyr::relocate(cluster, source1, source2, num.de, median.de, best, symbol) |>
            readr::write_tsv(file.path(outdir, name_clustering, paste0(name_clustering, "_doublet_detection_by_cluster.txt")))
        } else {
          print(paste("No doublet detection possible for", name_clustering, "because there are only", 
                length(unique(cluster)), "cluster (min 3 needed)."))
        }
        
      } 
  return(cluster)
  })
  
}
  

#############################
# save the sessionInformation and R image (no update of sce object necessary)
print("store results")
writeLines(capture.output(sessionInfo()),file.path(outdir, "kmeans_bioc_session_info.txt"))
save(params, cl, file=file.path(outdir,"kmeans_bioc.RData"))
