#####################################
##
## What: collectCT_bioc.R
## Who : Frank Rühle
## When: 21.08.2025
##
## Collects all cell type annotation tsv files and merges them to the sce object. 
## It checks if columns have already been added in previous runs.
## Additionally, creating boxplots per clustering and gene expression plots
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## seqtype          # sequencing type     
## resultsdir       # result directory of project
## ctanno_dir       # directory with cell type annotation resullts 
## selected_CTanno  # name(s) of CT annotation setting to use for downstream analysis
## exprPlotDir      # directory for plotting selected genes  
## selected_genes   # gene symbol(s) for expression plots (optional)
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
pipeline_root    <- parseArgs(args,"pipeline_root=") 
seqtype          <- parseArgs(args,"seqtype=")   
resultsdir       <- parseArgs(args,"res=")   
ctanno_dir       <- parseArgs(args,"ctanno_dir=") 
selected_CTanno  <- parseArgs(args,"selected_CTanno=", convert="as.char.vector") 
exprPlotDir      <- parseArgs(args,"exprPlotDir=")
selected_genes   <- parseArgs(args,"selected_genes=", convert="as.char.vector") 
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
print(paste("pipeline_root:", pipeline_root))
print(paste("seqtype:", seqtype))
print(paste("resultsdir:", resultsdir))
print(paste("ctanno_dir:", ctanno_dir))
print(paste("selected_CTanno:", paste0(selected_CTanno, collapse=", ")))
print(paste("exprPlotDir:", exprPlotDir))
print(paste("selected_genes:", paste0(selected_genes, collapse=", ")))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))



# load sce from previous module
sce <- HDF5Array::loadHDF5SummarizedExperiment(dir=file.path(resultsdir, "HDF5"), prefix="sce_")
print(sce)


# collect all cluster .tsv files in clusterdir sub-directories
files_tsv <- list.files(ctanno_dir, pattern = "anno.*\\.tsv$", full.names = TRUE, recursive = TRUE)
files_tsv_subdirs <- files_tsv[dirname(files_tsv) != normalizePath(ctanno_dir)]
print(paste("collect all cluster tsv files:", paste0(basename(files_tsv_subdirs), collapse=", ")))

merged_tsv <- files_tsv_subdirs |>
    purrr::set_names(basename) |>
    purrr::map(\(x) readr::read_tsv(x) |> dplyr::select(cell_id, starts_with("anno"))) |>
    purrr::reduce(dplyr::left_join, by="cell_id") |>
    readr::write_tsv(file.path(ctanno_dir,"summary_CTanno_files.txt"))
  

# add cluster columns to sce object
  merged_tsv_pruned <- merged_tsv |>
    dplyr::select(!any_of(colnames(SummarizedExperiment::colData(sce)) |> setdiff("cell_id"))) # skip annotations already present

  print(paste("add new clusterings to sce object:", paste(colnames(merged_tsv_pruned) |> setdiff("cell_id"), collapse=", ")))
  if(!all(colnames(sce)==merged_tsv_pruned$cell_id)) {stop("Cell IDs from sce object don't match Cell IDs from cell type annotation tsv files!")}
    SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce),dplyr::select(merged_tsv_pruned,-cell_id))

      


##### Expression plots for selected genes

# remove gene ids or symbols not in the data-set
rdat_filt <- SingleCellExperiment::rowData(sce)[,c("feature_id", "feature_symbol")] |>
  as.data.frame() |>
  dplyr::filter(feature_id %in% selected_genes | feature_symbol %in% selected_genes)

if(nrow(rdat_filt)>0) {
    
  if (!file.exists(file.path(exprPlotDir))) {dir.create(file.path(exprPlotDir), recursive=T) }

  list_boxplot <- purrr::map(selected_CTanno, function(selCT) {
    
      # store table with mean logcounts per celltype (all genes)
      meanLogcountsByCluster <- data.frame(celltype=SingleCellExperiment::colData(sce)[,selCT], t(SingleCellExperiment::logcounts(sce))) |>
        dplyr::group_by(celltype) |>
        dplyr::summarize_all(mean, na.rm = T) |>
        dplyr::ungroup() |>
        tibble::column_to_rownames('celltype') |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("gene") |>
        dplyr::mutate(symbol=SingleCellExperiment::rowData(sce)$feature_symbol) |>
        dplyr::relocate(gene, symbol) |>
        readr::write_tsv(file.path(exprPlotDir, paste0(selCT, "_mean_logcounts_all_genes.txt")))
      
      # get logcounts for all boxplots    
      data_boxplot <- t(SingleCellExperiment::logcounts(sce)[rdat_filt$feature_id,]) |>
        as.data.frame() |>
        dplyr::mutate(!!annocat_plot  := factor(SummarizedExperiment::colData(sce)[,annocat_plot]),
                      !!annocat_plot2 := factor(SummarizedExperiment::colData(sce)[,annocat_plot2]),
                      !!selCT      := factor(SummarizedExperiment::colData(sce)[,selCT])) |>
        tidyr::pivot_longer(cols=!any_of(c(selCT, annocat_plot, annocat_plot2)))
      
      # box plots of selected genes per celltype
      print(paste("Create expression boxplots for", selCT))
      
      # determine ideal number of panel columns in violin plots when using facet_wrap to avoid too many violin plots per row.
      maxViolPerRow = max(nlevels(factor(SummarizedExperiment::colData(sce)[,annocat_plot])), 10) # i.e. max #genes x #groups. 
      nGridCol <- max(1,floor(maxViolPerRow/nlevels(factor(SummarizedExperiment::colData(sce)[,annocat_plot])))) # number of cols in violin plot
      nGridRow <- ceiling(nlevels(factor(SummarizedExperiment::colData(sce)[,selCT])) / nGridCol) # number of rows in violin plot
      
      list_bplot <- purrr::pmap(rdat_filt , function(feature_id, feature_symbol) {
  
        data_boxplot_current <- data_boxplot |>
          dplyr::filter(name==feature_id)
        
        bplot <- ggplot(data_boxplot_current, aes(x=!!dplyr::sym(annocat_plot),y=value)) +
          geom_violin() +
          ggbeeswarm::geom_quasirandom(aes(color=!!dplyr::sym(annocat_plot2)), size = plot_pointsize, alpha=plot_pointalpha) +
          scale_color_hue(l=55) +
          geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA) +
          facet_wrap(as.formula(paste("~", selCT)), ncol=nGridCol) + 
          ylab("expression (logcounts)") +
          xlab(element_blank()) +
          ggtitle(paste0(feature_symbol, " expression by\n", selCT)) + 
          theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), legend.position = "bottom") + 
          guides(color=guide_legend(override.aes = list(size=1, alpha=1)))

        ggsave(plot=bplot, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expression_boxplot_by_", selCT, ".png")), 
               device="png", width=7, height=1+2*nGridRow, bg = "white") # legend + 2 inch per panel row
        ggsave(plot=bplot, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expression_boxplot_by_", selCT, ".pdf")), 
               device="pdf", width=7, height=1+2*nGridRow) # legend + 2 inch per panel row
        
        return(bplot)
      })
      return(list_bplot)
  })   
      
  
  # Create expression plots illustrated by all reduced dimensions and split by annotation categories 
  # These plots are independent from cell type annotation
  print("Create expression plots by all reduced dimensions and split by annotation categories")
  list_exprplot <- purrr::pmap(rdat_filt , function(feature_id, feature_symbol) {
    
    rdplots <- purrr::map(SingleCellExperiment::reducedDimNames(sce), function(dimred) {
      
      plot_data <- SingleCellExperiment::reducedDim(sce, dimred) |>
        as.data.frame() |>
        dplyr::select(1:2) |> 
        stats::setNames(c(paste0(gsub("_.*$", "", dimred), " 1"), paste0(gsub("_.*$", "", dimred), " 2"))) |>
        dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot],
                      !!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2],
                      !!feature_id := SingleCellExperiment::logcounts(sce)[feature_id,])
  
      rdplot <- ggplot(plot_data, aes(x = !!dplyr::sym(names(plot_data)[1]), y = !!dplyr::sym(names(plot_data)[2]), color = !!dplyr::sym(feature_id))) +
        geom_point(size = plot_pointsize, alpha=plot_pointalpha) +
        viridis::scale_color_viridis() + 
        ggtitle(paste("Log expression", feature_symbol, "by", dimred)) + 
        theme(legend.position = "bottom")
      
      rdplot_anno1 <- rdplot + 
        facet_wrap(as.formula(paste("~", annocat_plot2)))
      
      rdplot_anno2 <- rdplot + 
        facet_wrap(as.formula(paste("~", annocat_plot)))
      
      ggsave(plot=rdplot, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expressionplot_by_", dimred, ".png")), device="png", width=7, height=8, bg = "white")
      ggsave(plot=rdplot, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expressionplot_by_", dimred, ".pdf")), device="pdf", width=7, height=8)
      
      plotlayout <- ggplot2::ggplot_build(rdplot_anno1)$layout$layout
      ggsave(plot=rdplot_anno1, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expressionplot_by_", dimred, "_split_by_", annocat_plot2, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
      ggsave(plot=rdplot_anno1, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expressionplot_by_", dimred, "_split_by_", annocat_plot2, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
      
      plotlayout <- ggplot2::ggplot_build(rdplot_anno2)$layout$layout
      ggsave(plot=rdplot_anno2, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expressionplot_by_", dimred, "_split_by_", annocat_plot, ".png")), device="png", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)), bg = "white")
      ggsave(plot=rdplot_anno2, filename= file.path(exprPlotDir, paste0(feature_symbol, "_", feature_id, "_expressionplot_by_", dimred, "_split_by_", annocat_plot, ".pdf")), device="pdf", width=2.5*length(unique(plotlayout$COL)), height=2+2.5*length(unique(plotlayout$ROW)))
      return(rdplot)
    })
    return(rdplots)
  })

}
    
    

#############################
# save the sessionInformation and R image
print("store results")
writeLines(capture.output(sessionInfo()),file.path(ctanno_dir, "collectCT_bioc_session_info.txt"))
HDF5Array::quickResaveHDF5SummarizedExperiment(sce)
save(args, merged_tsv, file=file.path(ctanno_dir,"collectCT_bioc.RData"))
