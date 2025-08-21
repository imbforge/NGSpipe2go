#####################################
##
## What: scType_bioc.R
## Who : Sivarajan Karunanithi, Frank Rühle
## When: 21.08.2025
##
## Cell type annotation with scType
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## selected_clustering # name(s) of clustering setting to use for downstream analysis
## dbfile           # path to Excel file with marker genes in scType format
## tissueType       # tissue type to use in dbfile
## ctcolumn         # column name in dbfile with cell type names to use for results
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
selected_clustering  <- parseArgs(args,"selected_clustering=", convert="as.char.vector") 
dbfile           <- parseArgs(args,"dbfile=") 
tissueType       <- parseArgs(args,"tissueType=", convert="as.char.vector") 
ctcolumn         <- parseArgs(args,"ctcolumn=") 
annocat_plot     <- parseArgs(args,"annocat_plot=", default = "group")
annocat_plot2    <- parseArgs(args,"annocat_plot2=", default = "sample")
plot_pointsize   <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha  <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  


# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/CTanno", "bioc_3.16.lock")
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
print(paste("selected_clustering:", paste0(selected_clustering, collapse=", ")))
print(paste("dbfile:", dbfile))
print(paste("tissueType:", paste0(tissueType, collapse=", ")))
print(paste("ctcolumn:", ctcolumn))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- HDF5Array::loadHDF5SummarizedExperiment(dir=file.path(resultsdir, "HDF5"), prefix="sce_")
print(sce)

# load cell type annotation and helper function
source(file.path(getwd(), pipeline_root, "tools/CTanno", "scType_helper.R"))
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# set up combinations of clustering parameter (catch e.g. TSNE subtypes as well)
params <- tidyr::expand_grid( 
  selclust = selected_clustering, 
  ttype = tissueType
)


if(any(is.na(params))) { # skip entirely if cluster setting not specified
  ct <- list()
} else {

  ct <- purrr::pmap(params, function(selclust, ttype) {

      name_ct <- paste0("anno_scType_", selclust, "_tissue_", ttype)
      print(paste("Processing", name_ct, "using marker list in", dbfile))
      
      if(file.exists(file.path(outdir, name_ct, paste0(name_ct, ".tsv"))) | !selclust %in% colnames(SummarizedExperiment::colData(sce))) {
        
        print(paste0(name_ct, ".tsv already exists or clustering not in sce object. Cell type annotation for this setting is skipped"))
        ct_coldata <- NULL
        
      } else {
        
        if (!file.exists(file.path(outdir, name_ct))) {dir.create(file.path(outdir, name_ct), recursive=T) }
        
        # prepare gene sets
        gs_list = gene_sets_prepare(path_to_db_file=dbfile, tissue_type=ttype, name_type=ctcolumn)
        
        scRNAseqData = SingleCellExperiment::logcounts(sce) |>
          as.data.frame() |>
          tibble::rownames_to_column(var="feature_id") |>
          dplyr::left_join(SummarizedExperiment::rowData(sce) |> as.data.frame() |> dplyr::select(feature_id, feature_symbol), by="feature_id") |>
          tibble::column_to_rownames(var="feature_symbol")|>
          dplyr::select(!contains("feature_id")) |>
          as.matrix() # sctype_score needs matrix
        
        # get cell-type by cell matrix
        es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = FALSE,
                              gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
        
        # Extract colData and cluster assignments
        clust_assign <- as.vector(SummarizedExperiment::colData(sce)[, selclust])
        names(clust_assign) <- colnames(sce)
        # For each cluster, summarize top scores and combine
        cL_results <- purrr::map_dfr(unique(clust_assign), function(cl) {
          cells_in_cl <- names(clust_assign[clust_assign == cl])
          es_max_cl <- rowSums(es.max[, cells_in_cl])
          es_max_cl_sorted <- sort(es_max_cl, decreasing = TRUE)
          head(tibble::tibble(cluster = cl,
                 type = names(es_max_cl_sorted), 
                 scores = es_max_cl_sorted, 
                 ncells = length(cells_in_cl)),
               min(10, length(es_max_cl_sorted)))
        })

        # store cell type annotation overview
        sctype_scores = cL_results |> 
          dplyr::group_by(cluster) |> 
          dplyr::top_n(n = 1, wt = scores) |>
          dplyr::mutate(type= ifelse(scores < ncells/4, "Unknown", type)) |> # set low-confident (low ScType score) clusters to "unknown"
          dplyr::rename(!!dplyr::sym(selclust) := cluster, !!dplyr::sym(name_ct) := type) |>
          readr::write_tsv(file.path(outdir, name_ct, paste0(name_ct, "_overview.txt")))
        
        # store cell type annotations tsv
        ct_coldata <- SummarizedExperiment::colData(sce)[, c("cell_id", "sample", "group", selclust)] |>
          as.data.frame() |>
          dplyr::left_join(sctype_scores, by=selclust) |>
          readr::write_tsv(file.path(outdir, name_ct, paste0(name_ct, ".tsv")))
        
        # store cell counts per celltype and sample
        cells_per_sample <- as.data.frame.matrix(table(ct_coldata$sample, ct_coldata[,name_ct])) |>
          t() |>
          as.data.frame() |>
          tibble::rownames_to_column(var="celltype") |> 
          readr::write_tsv(file.path(outdir, name_ct, paste0(name_ct, "_cell_counts_per_celltype.txt")))
        

        # create cluster plots illustrated by all reduced dimensions and incl plots split by annotation categories
        purrr::map(SingleCellExperiment::reducedDimNames(sce), function(dimred) {
          
          plot_data <- SingleCellExperiment::reducedDim(sce, dimred) |>
            as.data.frame() |>
            dplyr::select(1:2) |> 
            stats::setNames(c(paste0(gsub("_.*$", "", dimred), " 1"), paste0(gsub("_.*$", "", dimred), " 2"))) |>
            dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot],
                          !!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2],
                          celltype = factor(ct_coldata[,name_ct]))
          
          rdplot <- ggplot(plot_data, aes(x = !!dplyr::sym(names(plot_data)[1]), y = !!dplyr::sym(names(plot_data)[2]), color = celltype)) +
            geom_point(size = plot_pointsize, alpha=plot_pointalpha) +
            scale_color_hue(l=55) +
            ggtitle(paste(gsub("_", " ", name_ct))) + 
            theme(legend.position = "bottom", ) + 
            guides(color=guide_legend(override.aes = list(size=1, alpha=1), title.position = "top"))
          
          rdplot_anno1 <- rdplot + 
            facet_wrap(as.formula(paste("~", annocat_plot2)))
          
          rdplot_anno2 <- rdplot + 
            facet_wrap(as.formula(paste("~", annocat_plot)))
          
          ggsave(plot=rdplot, filename= file.path(outdir, name_ct, paste0(name_ct, "_by_", dimred, ".png")), device="png", width=7, height=8, bg = "white")
          ggsave(plot=rdplot, filename= file.path(outdir, name_ct, paste0(name_ct, "_by_", dimred, ".pdf")), device="pdf", width=7, height=8)
          
          plotlayout <- ggplot2::ggplot_build(rdplot_anno1)$layout$layout
          ggsave(plot=rdplot_anno1, filename= file.path(outdir, name_ct, paste0(name_ct, "_by_", dimred, "_split_by_", annocat_plot2, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
          ggsave(plot=rdplot_anno1, filename= file.path(outdir, name_ct, paste0(name_ct, "_by_", dimred, "_split_by_", annocat_plot2, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
          
          plotlayout <- ggplot2::ggplot_build(rdplot_anno2)$layout$layout
          ggsave(plot=rdplot_anno2, filename= file.path(outdir, name_ct, paste0(name_ct, "_by_", dimred, "_split_by_", annocat_plot, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
          ggsave(plot=rdplot_anno2, filename= file.path(outdir, name_ct, paste0(name_ct, "_by_", dimred, "_split_by_", annocat_plot, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
        })
      } 
  return(ct_coldata)
  })
  
}
  

#############################
# save the sessionInformation and R image (no update of sce object necessary)
print("store results")
writeLines(capture.output(sessionInfo()),file.path(outdir, "scType_bioc_session_info.txt"))
save(args, params, ct, file=file.path(outdir,"scType_bioc.RData"))
