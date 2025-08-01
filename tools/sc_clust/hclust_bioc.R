#####################################
##
## What: hclust_bioc.R
## Who : Frank Rühle, Patrick Hüther
## When: 18.07.2025
##
## Apply hierarchical clustering on singlecellexperiment object
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## data2clust_pre   # data used for clustering. Can be any reduced dimension slot or 'logcounts'.
## hclust_method    # clustering method ("ward.D", "ward.D2", "single", "complete", average", "mcquitty", "median" or "centroid")
## ds               # deepSplit (range 0 to 4). The higher the value, the more and smaller clusters will be produced.
## minClusterSize   # minimum cluster size
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
hclust_method    <- parseArgs(args,string="hclust_method=", convert="as.char.vector") 
ds               <- parseArgs(args,string="deepSplit=",convert="run_custom_code", default=1)
minClusterSize   <- parseArgs(args,string="minClusterSize=", convert="as.numeric", default = 50)
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
print(paste("hclust_method:", paste0(hclust_method, collapse=", ")))
print(paste("ds:", paste0(ds, collapse=", ")))
print(paste("minClusterSize:", paste0(minClusterSize, collapse=", ")))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- HDF5Array::loadHDF5SummarizedExperiment(dir=file.path(resultsdir, "HDF5"), prefix="sce_")
print(sce)


# set up combinations of clustering parameter (catch e.g. TSNE subtypes as well)
data2clust_ext = unique(data2clust_pre, unlist(sapply(data2clust_pre, grep, SingleCellExperiment::reducedDimNames(sce), value=T))) 
params <- tidyr::expand_grid( 
  ds = ds,
  hclust_method = hclust_method
)


if(any(is.na(params), is.na(data2clust_ext))) { # skip entirely if cluster setting not specified
  cl <- list()
} else {

 cl <- purrr::map(data2clust_ext, function(data2clust) { 
    # create distance matrix separately from other params to avoid unnecessary re-creation for all parameter combinations 
    print(paste("Create dist matrix for", data2clust))
    dst <- if(data2clust != "logcounts") {
      SingleCellExperiment::reducedDim(sce,data2clust) |> dist()
    } else {
      SingleCellExperiment::logcounts(sce) |> t() |> dist()
    } 

    cl <- purrr::pmap(params, function(ds, hclust_method) {
        print(paste0("Process hierachical clustering based on ", data2clust,": ", hclust_method, " ds", ds))
  
  
      name_clustering <- paste0("cluster_hclust_", data2clust, "_", hclust_method, "_ds", ds)
      
      if(file.exists(file.path(outdir, name_clustering, paste0(name_clustering, ".tsv")))) {
        
        print(paste0(name_clustering, ".tsv already exists. Clustering for this setting is skipped"))
        cluster <- NULL
        
      } else {
      
        print(paste("Processing", name_clustering))
        if (!file.exists(file.path(outdir, name_clustering))) {dir.create(file.path(outdir, name_clustering), recursive=T) }
        
        cluster <- dynamicTreeCut::cutreeDynamic(hclust(dst,method=hclust_method),distM=as.matrix(dst),method="hybrid",minClusterSize=minClusterSize,deepSplit=ds,verbose=0)
        
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
            scale_fill_hue(l=55) +
            ggtitle(paste("Cluster hclust", data2clust, hclust_method, paste0("ds", ds), "by", dimred)) + 
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
        
        # Doublet detection by cluster (findDoubletClusters needs at least 3 clusters to detect doublet clusters)
        if(length(unique(cluster)) >=3) {
          print(paste("Cluster doublet detection for", name_clustering))
          dbl.out <- scDblFinder::findDoubletClusters(sce, clusters=cluster) |>
            as.data.frame() |> 
            tibble::rownames_to_column("cluster") |>
            dplyr::mutate(symbol=SummarizedExperiment::rowData(sce)[best, "feature_symbol"]) |>
            dplyr::mutate(dplyr::across(c(p.value, lib.size1, lib.size2, prop), signif, digits=3)) |>
            dplyr::mutate(suspicious=scater::isOutlier(num.de, nmads=3, type="lower", log=TRUE)) |>
            dplyr::relocate(cluster, source1, source2, num.de, median.de, best, symbol) |>
            readr::write_tsv(file.path(outdir, name_clustering, paste0(name_clustering, "_doublet_detection_by_cluster.txt")))
        } else {
          print(paste("No doublet detection possible for", name_clustering, "because there are only", 
                      length(unique(cluster)), "cluster (min 3 needed)."))
        }
        
      }  
    return(cluster)
  }) 
 }) 
}
  

#############################
# save the sessionInformation and R image (no update of sce object necessary)
print("store results")
writeLines(capture.output(sessionInfo()),file.path(outdir, "hclust_bioc_session_info.txt"))
save(params, cl, file=file.path(outdir,"hclust_bioc.RData"))
