#####################################
##
## What: findmarkers_bioc.R
## Who : Frank Rühle
## When: 14.08.2025
##
## Identify marker genes per cluster and run GO enrichment analysis
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## selected_clustering # data used for clustering. Can be any reduced dimension slot or 'logcounts'.
## rank_effectsize  # effect size to rank marker genes per cluster (e.g. "mean.AUC", "mean.logFC.cohen").
## block_var        # name of block var to account for uninteresting factor. Leave empty to skip.
## maxgenes2plot    # max number of genes to plot (expression plots).
## maxgenes_hm      # max number of f genes to plot (heatmaps).
## org              # either "human" or "mouse". Organism name needed for GO annotation. Leave empty to skip.
## top_genes_for_GO # number of top marker genes to use in GO enrichment per cluster.
## p_threshold_GO   # p-value threshold for over-representation of the GO term in the set.
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
resultsdir          <- parseArgs(args,"res=")   
seqtype             <- parseArgs(args,"seqtype=")   
outdir              <- parseArgs(args,"outdir=") # output folder
pipeline_root       <- parseArgs(args,"pipeline_root=") 
selected_clustering <- parseArgs(args,"selected_clustering=", convert="as.char.vector") 
rank_effectsize     <- parseArgs(args,"rank_effectsize=", convert="as.char.vector") 
block_var           <- parseArgs(args,"block_var=")
maxgenes2plot       <- parseArgs(args,"maxgenes2plot=", convert="as.numeric") 
maxgenes_hm         <- parseArgs(args,"maxgenes_hm=", convert="as.numeric")
org                 <- parseArgs(args,"org=")
top_genes_for_GO    <- parseArgs(args,"top_genes_for_GO=", convert="as.numeric")
p_threshold_GO      <- parseArgs(args,"maxgenes_hm=", convert="as.numeric")
annocat_plot        <- parseArgs(args,"annocat_plot=", default = "group")
annocat_plot2       <- parseArgs(args,"annocat_plot2=", default = "sample")
plot_pointsize      <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha     <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  

# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_findmarkers", "bioc_3.16.lock")
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
print(paste("rank_effectsize:", paste0(rank_effectsize, collapse=", ")))
print(paste("block_var:", block_var))
print(paste("maxgenes2plot:", maxgenes2plot))
print(paste("maxgenes_hm:", maxgenes_hm))
print(paste("org:", org))
print(paste("top_genes_for_GO:", top_genes_for_GO))
print(paste("p_threshold_GO:", p_threshold_GO))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- HDF5Array::loadHDF5SummarizedExperiment(dir=file.path(resultsdir, "HDF5"), prefix="sce_")
print(sce)


# set up combinations of clustering parameter
params <- tidyr::expand_grid( 
  selclust = selected_clustering, 
  es = rank_effectsize
)


if(any(is.na(params))) { # skip entirely if cluster setting not specified
  print("Parameter missing! Marker detection skipped")
  mg <- list()
} else {

  mg <- purrr::pmap(params, function(selclust, es) {

      name_markerlist <- paste0("marker_genes_", selclust, "_ranked_by_", es)
      SingleCellExperiment::colLabels(sce) <- SummarizedExperiment::colData(sce)[,selclust]
      
      if(file.exists(file.path(outdir, selclust, paste0(name_markerlist, "_cluster1.txt"))) | nlevels(factor(SingleCellExperiment::colLabels(sce))) < 2) {
        
        print(paste0(selclust, " already exists or has got just one cluster. Marker gene detection for this setting is skipped"))
        markers_l <- NULL
        
      } else {
        
        print(paste("Processing", name_markerlist))
        if (!dir.exists(file.path(outdir, selclust))) {dir.create(file.path(outdir, selclust), recursive=T) }
        
          markers_l <- scran::scoreMarkers(sce, groups=SingleCellExperiment::colLabels(sce), 
                                           block = if(!is.na(block_var)) {SummarizedExperiment::colData(sce)[,block_var]} else {NULL},
                                           assay.type = "logcounts",
                                           row.data=SummarizedExperiment::rowData(sce)[,"feature_symbol", drop=F]) 
          
          
          purrr::map(names(markers_l), function(c) {
            
            print(paste("Processing cluster", c))
            
            markers_per_cluster <- markers_l[[c]] |> 
              as.data.frame() |>
              tibble::rownames_to_column("feature_id") |>
              dplyr::arrange(dplyr::desc(!!dplyr::sym(es))) |>
              readr::write_tsv(file.path(outdir, selclust, paste0(name_markerlist, "_cluster", c, ".txt")))
            
            print(paste("Expression plots for top marker genes for", selclust, "ranked by", es))
            maxViolPerRow = max(nlevels(factor(SingleCellExperiment::colLabels(sce))), 15) # i.e. max #genes x #cluster 
            nGridCol <- max(1,floor(maxViolPerRow/nlevels(factor(SingleCellExperiment::colLabels(sce))))) # number of cols in violin plot
            nGridRow <- ceiling(maxgenes2plot / nGridCol) # number of rows in violin plot
            
            top_marker_data <- SingleCellExperiment::logcounts(sce)[markers_per_cluster$feature_id[1:maxgenes2plot],] |>
              as.matrix() |>
              t() |>
              as.data.frame() |>
              dplyr::mutate(!!annocat_plot := SummarizedExperiment::colData(sce)[,annocat_plot]) |>
              dplyr::mutate(!!annocat_plot2 := SummarizedExperiment::colData(sce)[,annocat_plot2]) |>
              dplyr::mutate(!!selclust := factor(SummarizedExperiment::colData(sce)[,selclust])) |>
              tidyr::pivot_longer(cols=!any_of(c(annocat_plot, annocat_plot2, selclust))) |>
              dplyr::mutate(symbol=SummarizedExperiment::rowData(sce)$feature_symbol[match(name, rownames(sce))])
            
            top_marker_plot <- ggplot(top_marker_data, aes(x=!!dplyr::sym(selclust),y=value)) +
              geom_violin() +
              ggbeeswarm::geom_quasirandom(aes(color=!!dplyr::sym(annocat_plot2)), size = plot_pointsize, alpha=plot_pointalpha) +
              scale_color_hue(l=55) +
              geom_boxplot(color = "darkgrey", alpha = 0.2, outlier.shape = NA) +
              facet_wrap(as.formula(paste("~", "symbol")), ncol=nGridCol) +
              ylab("expression (logcounts)") +
              xlab(selclust) +
              ggtitle(paste("Top", maxgenes2plot, "marker genes ranked by", es, "for cluster", c)) + 
              theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), legend.position = "bottom") + 
              guides(color=guide_legend(override.aes = list(size=1, alpha=1)))
            
            ggsave(plot=top_marker_plot, width=7, height=1+2*nGridRow, # legend + 2 inch per panel row
                   filename= file.path(outdir, selclust, paste0(name_markerlist, "_expression_plot_cluster", c, ".png")), device="png", bg = "white")
            # ggsave(plot=top_marker_plot, width=7, height=1+2*nGridRow, # legend + 2 inch per panel row
            #        filename= file.path(outdir, selclust, paste0(name_markerlist, "_expression_plot_cluster", c, ".pdf")), device="pdf")
            
            
            print(paste("Expression heatmaps for top marker genes for", selclust, "ranked by", es))
            hm <- scater::plotGroupedHeatmap(sce, features=markers_per_cluster$feature_symbol[1:maxgenes_hm], 
                                             group=selclust,
                                             block=if(!is.na(block_var)) {SummarizedExperiment::colData(sce)[,block_var]} else {NULL},
                                             center=T, swap_rownames="feature_symbol",
                                             main=paste("Top", maxgenes_hm, "marker genes (centered mean expression) for cluster", c, "\n", selclust, "ranked by", es, ""))
              
            ggsave(plot=hm, width=7, height=1.5+0.3*maxgenes_hm, # dendrogram + 0.3 inch per gene
                   filename= file.path(outdir, selclust, paste0(name_markerlist, "_heatmap_cluster", c, ".png")), device="png", bg = "white")
            # ggsave(plot=hm, width=7, height=1.5+0.3*maxgenes_hm, # dendrogram + 0.3 inch per gene
            #        filename= file.path(outdir, selclust, paste0(name_markerlist, "_heatmap_cluster", c, ".pdf")), device="pdf")
            
            
            ## GO enrichment
            if(!is.na(org) && org %in% c("human", "mouse")) {
              print(paste("GO enrichment analysis (BP) with top", top_genes_for_GO, "genes for", name_markerlist)) 
                
              switch(org,
                     human={
                       library(org.Hs.eg.db)
                       orgdb <- org.Hs.eg.db
                       species="Hs"
                     },
                     mouse={
                       library(org.Mm.eg.db)
                       orgdb <- org.Mm.eg.db
                       species="Mm"
                     })
              
              markers_per_cluster$entrez_ids <- AnnotationDbi::mapIds(orgdb, keys=gsub("\\..*$", "", markers_per_cluster$feature_id), 
                                                                      column="ENTREZID", keytype="ENSEMBL")
              markers_top <- markers_per_cluster[1:top_genes_for_GO,]
              print(paste(sum(is.na(markers_top$entrez_ids)), "marker genes skipped because no Entrez IDs found:", paste(markers_top$feature_symbol[is.na(markers_top$entrez_ids)], collapse=", ")))
              
              go_out <- limma::goana(unique(na.omit(markers_top$entrez_ids)), species=species,
                              universe=unique(na.omit(markers_per_cluster$entrez_ids)))
              
              go_out_filt <- go_out |> # Only keeping BP terms that are not overly general and which are significantly enriched.
                dplyr::filter(Ont=="BP" & N<=500 & P.DE < p_threshold_GO) |> # p-value for over-representation of the GO term in the set.
                tibble::rownames_to_column("GOID") |>
                dplyr::arrange(P.DE) |>
                dplyr::mutate(P.DE = signif(P.DE, digits=4)) |>
                readr::write_tsv(file.path(outdir, selclust, paste0("GOenrichment_top", top_genes_for_GO, "_", name_markerlist, "_cluster", c, ".txt")))
              
            }
          return(markers_per_cluster)  
          })
      } 
      return(markers_l)
      })
  }
  
          
#############################
# save the sessionInformation and R image (no update of sce object necessary)
print("store results")
writeLines(capture.output(sessionInfo()),file.path(outdir, "findmarkers_bioc_session_info.txt"))
save(args, params, mg, file=file.path(outdir,"findmarkers_bioc.RData"))
