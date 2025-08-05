#####################################
##
## What: filtCells_bioc.R
## Who : Frank Rühle, Patrick Hüther
## When: 04.06.2025
##
## Script to apply quality control thresholds to scRNA-Seq assay data.
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## samples2exclude            # string with sample names if entire samples shall be excluded from analysis before QC filtering.         
## type_of_threshold          # "absolute" or "relative". For "relative" NMADs are calculated per category.    
## threshold_total_counts_min # min read counts per cell (for absolute threshold only).
## threshold_total_counts_max # max read counts per cell (for absolute threshold only).
## threshold_total_detected   # min genes detected (for absolute threshold only).
## threshold_pct_counts_Mt    # max percentage of mitochondrial gene counts (for absolute threshold only).
## threshold_pct_counts_spikein # max percentage of spike-in gene counts (for absolute threshold only).
## NMADS                      # median absolute deviations (MAD) to define outliers (for relative threshold only).     
## category_NMADS             # grouping of cells used to calculate MAD. If empty, no grouping applied (for relative threshold only).        
## threshold_low_abundance    # threshold cell portion to exclude low abundance genes. E.g. 0.01 means filter out genes with no expression in 99% of cells.
## annocat_plot               # category used in plots and tables    
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
resultsdir                   <- parseArgs(args,"res=")   
outdir                       <- parseArgs(args,"outdir=") # output folder
pipeline_root                <- parseArgs(args,"pipeline_root=") 
samples2exclude              <- parseArgs(args,string="samples2exclude=", convert="run_custom_code", default=NULL)
type_of_threshold            <- parseArgs(args,string="type_of_threshold=")
threshold_total_counts_min   <- parseArgs(args,string="threshold_total_counts_min=",convert="as.numeric")
threshold_total_counts_max   <- parseArgs(args,string="threshold_total_counts_max=",convert="as.numeric")
threshold_total_detected     <- parseArgs(args,string="threshold_total_detected=",convert="as.numeric")
threshold_pct_counts_Mt      <- parseArgs(args,string="threshold_pct_counts_Mt=",convert="as.numeric")
threshold_pct_counts_spikein <- parseArgs(args,string="threshold_pct_counts_spikein=",convert="as.numeric", default=100)
NMADS                        <- parseArgs(args,string="NMADS=",convert="as.numeric")
category_NMADS               <- parseArgs(args,string="category=", default=NULL)
threshold_low_abundance      <- parseArgs(args,string="threshold_low_abundance=",convert="as.numeric")
annocat_plot                 <- parseArgs(args,"annocat_plot=", default="sample")
plot_pointsize               <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha              <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  


# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_qc", "bioc_3.16.lock")
print(paste("load renv:", env.path))
renv::use(lockfile=env.path)

library(ggplot2)

# set options
options(stringsAsFactors=FALSE)

# check parameter
print(paste("resultsdir:", resultsdir))
print(paste("outdir:", outdir))
print(paste("pipeline_root:", pipeline_root))
print(paste("samples2exclude:", paste0(samples2exclude, collapse=", ")))
print(paste("type_of_threshold:", type_of_threshold))
print(paste("threshold_total_counts_min:", threshold_total_counts_min))
print(paste("threshold_total_counts_max:", threshold_total_counts_max))
print(paste("threshold_total_detected:", threshold_total_detected))
print(paste("threshold_pct_counts_Mt:", threshold_pct_counts_Mt))
print(paste("threshold_pct_counts_spikein:", threshold_pct_counts_spikein))
print(paste("threshold_low_abundance:", threshold_low_abundance))
print(paste("NMADS:", NMADS))
print(paste("category_NMADS:", category_NMADS))
print(paste("annocat_plot:", annocat_plot))
print(paste("plot_pointsize:", plot_pointsize))
print(paste("plot_pointalpha:", plot_pointalpha))


# load sce from previous module
sce <- readr::read_rds(file.path(resultsdir, "sce_raw.RDS"))


# exclude entire samples if requested
if(!is.null(samples2exclude) && !is.na(samples2exclude)) {
cat(paste("Sample(s)", paste(samples2exclude, collapse=", "), "with", sum(SummarizedExperiment::colData(sce)$sample %in% samples2exclude), "cells removed from dataset.\n"))
write.table(data.frame(samples_excluded=samples2exclude,
                       cells=sapply(samples2exclude, function(x) {sum(SummarizedExperiment::colData(sce)$sample==x)})), 
            file= file.path(outdir, "samples_excluded.txt"), sep="\t", quote=F, row.names = F)              
sce <- sce[, !(sce$sample %in% samples2exclude)]
}


# get colData table from sce object to store qc flags
# mark SMART-Seq control wells if present (plate wells with 10 cells ('10c') instead of 1 cell ('1c') in column 'cells'))
qc.drop <- SummarizedExperiment::colData(sce) |> 
  tidyr::as_tibble() |>
  dplyr::mutate(control_wells=if('cells' %in% colnames(SummarizedExperiment::colData(sce))) cells!='1c' else F) 


## apply filtering thresholds
if(type_of_threshold=="absolute") {
  print("apply absolute thresholds")
  
  # overview thresholds:
  qc_thresholds <- data.frame(criterion=c(
    paste("type of threshold:", type_of_threshold), 
    paste("total counts >", threshold_total_counts_min), 
    paste("total counts <", threshold_total_counts_max), 
    paste("detected genes >", threshold_total_detected), 
    paste("MT count percent <", threshold_pct_counts_Mt), 
    if('subsets_spikein_percent' %in% colnames(qc.drop)) {paste("spikein count percent <", threshold_pct_counts_spikein)}  
  )) |>
    readr::write_tsv(file.path(outdir, "qc_thresholds.txt"))
  

  qc.drop <- qc.drop |>
    dplyr::mutate(libsize=!dplyr::between(sum,threshold_total_counts_min,threshold_total_counts_max)) |>
    dplyr::mutate(features=(detected < threshold_total_detected)) |>
    dplyr::mutate(mito=(subsets_Mito_percent > threshold_pct_counts_Mt)) |>
    dplyr::mutate(spikein=if('subsets_spikein_percent' %in% colnames(qc.drop)) subsets_spikein_percent > threshold_pct_counts_spikein else F) |>
    dplyr::mutate(pass = rowSums(dplyr::across(c(libsize, features, mito, spikein, control_wells))) == 0)

} else if(type_of_threshold=="relative") {
  print("apply relative thresholds")
  # ignore all cells which hardly have any counts to avoid bias for relative thresholds
  min_readcount <- 100
  count.drop <- (qc.drop$sum < min_readcount)
  batch_isOutlier <- if(!is.null(category_NMADS) && !is.na(category_NMADS)) {factor(dplyr::pull(qc.drop,category_NMADS))} else {NULL}
  
  # overview thresholds:
  qc_thresholds <- data.frame(criterion=c(
    paste("type of threshold:", type_of_threshold), 
    paste(sum(count.drop), "cells skipped from MAD calc with counts <", min_readcount), 
    paste("number of median absolute deviations (MAD) to define outliers per sample:", NMADS), 
    paste("grouping of cells used to calculate MAD:", if(is.null(category_NMADS) || is.na(category_NMADS)) "none" else category_NMADS) 
  )) |>
    readr::write_tsv(file.path(outdir, "qc_thresholds.txt"))
  
  qc.drop <- qc.drop |>
    dplyr::mutate(libsize=scater::isOutlier(sum,nmads=NMADS,type="lower",log=TRUE,subset=!count.drop,batch=batch_isOutlier)) |>
    dplyr::mutate(features=scater::isOutlier(detected,nmads=NMADS,type="lower",log=TRUE,subset=!count.drop,batch=batch_isOutlier)) |>
    dplyr::mutate(mito=scater::isOutlier(subsets_Mito_percent,nmads=NMADS,type="higher",subset=!count.drop,batch=batch_isOutlier)) |>
    dplyr::mutate(spikein=if('subsets_spikein_percent' %in% colnames(qc.drop)) scater::isOutlier(subsets_spikein_percent,nmads=NMADS,type="higher",subset=!count.drop,batch=batch_isOutlier) else F) |>
    dplyr::mutate(pass = rowSums(dplyr::across(c(libsize, features, mito, spikein, control_wells))) == 0)
  
} else {
  print("Filtering is skipped because no filtering specified!")
  qc_thresholds <- "none"
  qc.drop <- qc.drop |>
    dplyr::mutate(libsize=F,
                  features=F,
                  mito=F,
                  spikein=F,
                  pass=T) 
  }

write.table(qc.drop, file= file.path(outdir, "qc.drop.txt"), sep="\t", quote=F, row.names = F)              


## prepare overview table of filtered cells
print("prepare overview table of filtered cells")
format_pct <- function(n_pass, n_total) { # define formatting function
  paste0(n_pass, " (", scales::percent(n_pass / n_total, accuracy = 0.1), ")")
}

qcfailed <- qc.drop |>
  dplyr::group_by(!!dplyr::sym(annocat_plot)) |>
  dplyr::summarize(
    "counts unfilt." = format_pct(dplyr::n(), ncol(sce)),
    libsize          = format_pct(sum(libsize), dplyr::n()),
    "detected genes" = format_pct(sum(features), dplyr::n()),
    "MT count percent"  = format_pct(sum(mito), dplyr::n()),
    "spikein count percent" = format_pct(sum(spikein), dplyr::n()),
    "control wells"  = format_pct(sum(control_wells), dplyr::n()),
    remaining        = format_pct(sum(pass), dplyr::n()),
    .groups = "drop"
  ) |>
  dplyr::add_row(
    !!dplyr::sym(annocat_plot) := "cell count total",
    "counts unfilt." = format_pct(nrow(qc.drop), ncol(sce)),
    libsize          = format_pct(sum(qc.drop$libsize), nrow(qc.drop)),
    "detected genes" = format_pct(sum(qc.drop$features), nrow(qc.drop)),
    "MT count percent"  = format_pct(sum(qc.drop$mito), nrow(qc.drop)),
    "spikein count percent" = format_pct(sum(qc.drop$spikein), nrow(qc.drop)),
    "control wells"  = format_pct(sum(qc.drop$control_wells), nrow(qc.drop)),
    remaining        = format_pct(sum(qc.drop$pass), nrow(qc.drop))
  ) |>
  tibble::column_to_rownames(annocat_plot) |>
  t() |> 
  as.data.frame() |>
  tibble::rownames_to_column(var = "criterion")

if(!"subsets_spikein_percent" %in% colnames(SummarizedExperiment::colData(sce))) { # exclude control wells if not applicable
  qcfailed <- qcfailed |> dplyr::filter(criterion != "spikein count percent")
  qc.drop$spikein <- NULL
}
if(!"cells" %in% colnames(SummarizedExperiment::colData(sce))) { # exclude control wells if not applicable
  qcfailed <- qcfailed |> dplyr::filter(criterion != "control wells")
  qc.drop$control_wells <- NULL
}
write.table(qcfailed, file= file.path(outdir, "qcfailed_overview.txt"), sep="\t", quote=F, row.names = F)              



# violin plots with indicated thresholds
print("create violin plots with indicated thresholds")
qc_metrics <- c("sum", "detected", "subsets_Mito_percent", "subsets_spikein_percent")
qc_metrics <- qc_metrics[qc_metrics %in% colnames(qc.drop)]

qc.plots.violin <- lapply(qc_metrics, function(to.plot){ 
  column_applied_threshold <- switch(to.plot, "sum" = "libsize", "detected" = "features", "subsets_Mito_percent" = "mito", "subsets_spikein_percent" = "spikein")
  
  p <- ggplot(qc.drop, aes(!!sym(annocat_plot),!!sym(to.plot)))+ 
    geom_violin() +
    ggbeeswarm::geom_quasirandom(aes(color = !!sym(column_applied_threshold)), size = plot_pointsize, alpha=plot_pointalpha) +
    scale_color_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "orange3"),
      labels = c("FALSE" = "kept", "TRUE" = "removed"),
      name = element_blank()
    ) +
    xlab(element_blank()) +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
    guides(color=guide_legend(override.aes = list(size=1, alpha=1))) +
    ggtitle(dplyr::case_match(to.plot, "sum" ~ "Library Size",
                              "detected" ~ "Detected Genes",
                              "subsets_Mito_percent" ~ "Reads Mapping to Mitochondrial Genes in Percent",
                              "subsets_spikein_percent" ~ "Reads Mapping to spikein Genes in Percent"))
  
  ggsave(plot=p, width=7, height=5, filename= file.path(outdir, paste0("qc_violin_plot_", to.plot, "_filtered.png")), device="png", bg = "white")
  ggsave(plot=p, width=7, height=5, filename= file.path(outdir, paste0("qc_violin_plot_", to.plot, "_filtered.pdf")), device="pdf")
  return(p)
})


# filter sce object for cells failing QC
print("apply filtering to sce object")
if(length(qc.drop$pass) == ncol(sce)) { # if chunk is executed multiple times
  sce <- sce[, qc.drop$pass]
} else {stop("length(qc.drop$pass) != ncol(sce). Check QC filtering!")}



## exclude low abundance genes
print("exclude low abundance genes")
SummarizedExperiment::rowData(sce)$expressed_cells <- scater::nexprs(sce, byrow=TRUE) # Counting the number of non-zero counts in each row (per feature)
# get genes expressed in min threshold_low_abundance of cells
genes2keep <- SummarizedExperiment::rowData(sce)$expressed_cells > ceiling(threshold_low_abundance * ncol(sce)) 
# The size factor-adjusted average count is defined by dividing each count by the size factor and taking the average across cells.
# If no size factors given, they are calculated or extracted from sce object.
ave_count_before_low_abundance_filt <- scater::calculateAverage(sce) 

png(file=file.path(outdir, "low_abundance.png"), width=7, height=4, units="in", res=300) 
  smoothScatter(log10(ave_count_before_low_abundance_filt), SummarizedExperiment::rowData(sce)$expressed_cells,
                xlab=expression("Log10 average count"), ylab= "Number of expressing cells")
  # is.ercc <- isSpike(sce, type="ERCC")
  # points(log10(ave_count_before_low_abundance_filt[is.ercc]), numcells[is.ercc], col="red", pch=16, cex=0.5)
dev.off()

# store excluded genes and logical vector genes2keep
write.table(SummarizedExperiment::rowData(sce)[!genes2keep,], file =file.path(outdir, "low_abundance_excluded_genes.txt"), sep="\t", quote = F, row.names = F)
write.table(data.frame(c("total genes", "genes kept", "threshold"), 
                       c(round(length(genes2keep),0), round(sum(genes2keep),0), threshold_low_abundance)), 
            file =file.path(outdir, "low_abundance_count.txt"), sep="\t", quote = F, row.names = F, col.names = F)

# remove low abundance genes from sce object
sce <- sce[genes2keep, ]


#############################
# save the sessionInformation and R image
print("store data")
writeLines(capture.output(sessionInfo()),paste0(outdir, "/filtCells_bioc_session_info.txt"))
readr::write_rds(sce, file = file.path(resultsdir, "sce.RDS"))
save(qc.drop, qcfailed, samples2exclude, qc_thresholds, genes2keep, threshold_low_abundance, file=paste0(outdir,"/filtCells_bioc.RData"))

