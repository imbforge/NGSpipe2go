#####################################
##
## What: sc_bioc_qc.R
## Who : Frank Rühle, Patrick Hüther
## When: 04.06.2025
##
## Script to perform quality control metrics of scRNA-Seq data.
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## mito.genes       # list with gene_ids of mitochondrial genes. If empty, standard mito chr names used
## annocat_plot     # category used for plotting
## annocat_plot2    # 2nd category used for plotting
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

args <- commandArgs(T)
resultsdir    <- parseArgs(args,"res=")   
seqtype       <- parseArgs(args,"seqtype=")   
outdir        <- parseArgs(args,"outdir=") # output folder
pipeline_root <- parseArgs(args,"pipeline_root=") 
mito.genes    <- parseArgs(args,"mito_genes=", default = "")
annocat_plot  <- parseArgs(args,"annocat_plot=", default = "group")
annocat_plot2 <- parseArgs(args,"annocat_plot2=", default = "group")
plot_pointsize   <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha  <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  


# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_qc", "bioc_3.16.lock")
print(paste("load renv:", env.path))
renv::use(lockfile=env.path)

library(ggplot2)

# set options
options(stringsAsFactors=FALSE)

# check parameter
print(paste("resultsdir:", resultsdir))
print(paste("seqtype:", seqtype))
print(paste("outdir:", outdir))
print(paste("pipeline_root:", pipeline_root))
print(paste("mito.genes:", mito.genes))
print(paste("annocat_plot:", annocat_plot))
print(paste("annocat_plot2:", annocat_plot2))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- readr::read_rds(file.path(resultsdir, "sce.RDS"))

# get mitochondrial genes
print("get mitochondrial genes")
if(file.exists(file.path(mito.genes))) {
  # use predefined list with mitochondrial genes
  mito.genes <- readr::read_tsv(file=file.path(mito.genes), col_names =F) |> dplyr::pull(X1)
  is.mito  <- row.names(sce) %in% mito.genes 
  print(paste(sum(is.mito), "mitochondrial genes identified by pre-defined mitochondrial gene list."))
} else {
    # use default chr names (mind that some assays have prefixes like 'mm39_' in chr names)
    default_mito_chr_names <- c("chrM", "chrMT", "M", "MT")
    is.mito <- gsub(".*_", "", SummarizedExperiment::rowData(sce)$feature_chrom) %in% default_mito_chr_names 
    print(paste(sum(is.mito), "mitochondrial genes identified on default mitochrondrial chromosome names:", paste(default_mito_chr_names, collapse = ", ")))
}
names(is.mito) <- row.names(sce)
write.table(is.mito, file=file.path(outdir, "is.mito.txt"), sep="\t", quote=F, row.names = T, col.names = F)

# calculate QC metrics
print("calculate QC metrics")
sce <- scuttle::addPerCellQCMetrics(sce,percent_top=2,subsets=list(Mito=is.mito))

# get colData table from sce object to store qc flags
qc.frame <- SummarizedExperiment::colData(sce)

# mark control cells if present
qc.frame <- tidyr::as_tibble(qc.frame) |> dplyr::mutate(controls=if('cells' %in% colnames(qc.frame)) cells!='1c' else F)
write.table(qc.frame, file=file.path(outdir, "qc.frame.txt"), sep="\t", quote=F, row.names = F)



# violinplots
print("create violinplots")
qc.plots.violin <- lapply(c("sum", "detected", "subsets_Mito_percent"), function(to.plot){ 

  p <- ggplot(qc.frame, aes(!!sym(annocat_plot2),!!sym(to.plot)))+ 
    geom_violin() +
    ggbeeswarm::geom_quasirandom(aes(color = !!sym(annocat_plot)), size = plot_pointsize, alpha=plot_pointalpha) +
    scale_color_hue(l=55) +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
    guides(color=guide_legend(override.aes = list(size=1, alpha=1))) +
    labs(title= if(to.plot=="sum") {"Library Size"} else {if(to.plot=="detected") {"Detected Genes"} else {"Reads Mapping to Mitochondrial Genes in Percent"}})
  
  ggsave(plot=p, width=7, height=5, filename= file.path(outdir, paste0("qc_violin_plot_", to.plot, ".png")), device="png", bg = "white")
  ggsave(plot=p, width=7, height=5, filename= file.path(outdir, paste0("qc_violin_plot_", to.plot, ".pdf")), device="pdf")
  return(p)
})



## Top 2% biggest libraries (based on mapped reads on features)
print("Top 2% biggest libraries")
highest.lib.size <- qc.frame |>
  dplyr::filter(sum > quantile(sum, 0.98)) |>
  dplyr::select(cell_id, sample, group, sum, detected, subsets_Mito_percent) |>
  dplyr::arrange(desc(sum)) |>
  dplyr::mutate(dplyr::across(c(sum, subsets_Mito_percent), round, 2))

write.table(highest.lib.size, file=file.path(outdir, "highest.lib.size.txt"), sep="\t", quote=F)


## Plot count distribution per plate position (for Smart-Seq)
if(seqtype %in% c("SmartSeq")) {
  print("Plot count distribution per plate position")
  
  sce$plate_position <- paste0(sce$row, sce$col) # column "plate_position" needed for plotPlatePosition
  
  for (size in c("sum", "detected", "subsets_Mito_percent")) {
    cat(paste("\nPlotting", size, "as spot size\n\n"))
    plates <- list()
    for (p in as.character(sort(unique(sce$plate)))) {
      plates[[p]] <- scater::plotPlatePosition(sce[, sce$plate==p], colour_by="cells",size_by=size, shape_by=annocat_plot,
                                               by_exprs_values = "counts", theme_size = 10,  point_alpha = 0.6, point_size = 3, add_legend = F) + 
        ggtitle(label=paste0("Plate ", p, ": ", size, " (range: ", paste(signif(range(SummarizedExperiment::colData(sce[, sce$plate==p])[,size]),2), collapse=" - "), ")")) +
        theme(plot.title = element_text(color="black", size=8))
    } # plotPlatePosition uses by_exprs_values = "logcounts" by default. But if not available, uses "counts" instead
  
    p <- gridExtra::grid.arrange(grobs=plates, layout_matrix=matrix(c(1:ceiling(length(plates))), ncol=2, byrow=TRUE))
    ggsave(plot=p, width=7, height=4*ceiling(lengh(plates[[p]])/2), filename= file.path(outdir, paste0("plate_distribution_plot_by_", size, ".png")), device="png", bg = "white")
    ggsave(plot=p, width=7, height=4*ceiling(lengh(plates[[p]])/2), filename= file.path(outdir, paste0("plate_distribution_plot_by_", size, ".pdf")), device="pdf")
  }
}


## Correlation plots for different features
print("create scatter plots")

  scatterplot.lib.size <- ggplot(qc.frame, aes(x=detected, y=sum)) +
    geom_point(aes(color=!!sym(annocat_plot)), size = plot_pointsize, alpha=plot_pointalpha) +
    guides(color=guide_legend(override.aes = list(size=1, alpha=1))) +
    scale_color_hue(l=55) +
    ylab("Library size") +
    xlab("Number of expressed genes")
  
  scatterplot.mito.perc <- scatterplot.lib.size + aes(y=subsets_Mito_percent) + ylab("% mitochondrial reads")
  
  ggsave(plot=scatterplot.lib.size, width=7, height=3, filename= file.path(outdir, paste0("scatterplot_lib_size.png")), device="png", bg = "white")
  ggsave(plot=scatterplot.lib.size, width=7, height=3, filename= file.path(outdir, paste0("scatterplot_lib_size.pdf")), device="pdf")
  ggsave(plot=scatterplot.mito.perc, width=7, height=3, filename= file.path(outdir, paste0("scatterplot_mito_perc.png")), device="png", bg = "white")
  ggsave(plot=scatterplot.mito.perc, width=7, height=3, filename= file.path(outdir, paste0("scatterplot_mito_perc.pdf")), device="pdf")
  


#############################
# save the sessionInformation and R image
print("store data")
writeLines(capture.output(sessionInfo()),paste0(outdir, "/sc_bioc_qc_session_info.txt"))
readr::write_rds(sce, file = file.path(resultsdir, "sce.RDS"))
save(qc.frame, highest.lib.size, file=paste0(outdir,"/sc_bioc_qc.RData"))

