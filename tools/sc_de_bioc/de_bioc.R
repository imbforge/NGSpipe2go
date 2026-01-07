#####################################
##
## What: de_bioc.R
## Who : Frank Rühle
## When: 19.08.2025
##
## Differential gene expression analysis per cluster and run GO enrichment analysis
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## contrasts        # path to contrasts.txt     
## selected_clustering # names of clustering setting to use.
## selected_CTanno  # names of CT annotation setting to use.
## minClusterSize   # minimum cells per pseudo-bulk sample in DE analysis. Leave empty to skip.
## lfc_threshold_DE # log-fold change threshold for hypothesis testing (not a post-hoc filter) (default: "0").
## FDR_threshold_DE # FDR threshold for filtering DE gene result table. Leave empty to skip.
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
outdir              <- parseArgs(args,"outdir=") 
pipeline_root       <- parseArgs(args,"pipeline_root=") 
contrasts           <- parseArgs(args,"contrasts=") 
selected_clustering <- parseArgs(args,"selected_clustering=", convert="as.char.vector") 
selected_CTanno     <- parseArgs(args,"selected_CTanno=", convert="as.char.vector") 
minClusterSize      <- parseArgs(args,"minClusterSize=", convert="as.numeric") 
lfc_threshold_DE    <- parseArgs(args,"lfc_threshold_DE=", convert="as.numeric", default=0) 
FDR_threshold_DE    <- parseArgs(args,"FDR_threshold_DE=", convert="as.numeric", default=1)
org                 <- parseArgs(args,"org=")
top_genes_for_GO    <- parseArgs(args,"top_genes_for_GO=", convert="as.numeric")
p_threshold_GO      <- parseArgs(args,"p_threshold_GO=", convert="as.numeric")
annocat_plot        <- parseArgs(args,"annocat_plot=", default = "group")
annocat_plot2       <- parseArgs(args,"annocat_plot2=", default = "sample")
plot_pointsize      <- parseArgs(args,"plot_pointsize=", convert="as.numeric", default = 0.6) 
plot_pointalpha     <- parseArgs(args,"plot_pointalpha=", convert="as.numeric", default = 0.6)  


# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_de_bioc", "bioc_3.16.lock")
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
print(paste("contrasts:", contrasts))
print(paste("selected_clustering:", paste0(selected_clustering, collapse=", ")))
print(paste("selected_CTanno:", paste0(selected_CTanno, collapse=", ")))
print(paste("minClusterSize:", minClusterSize))
print(paste("lfc_threshold_DE:", lfc_threshold_DE))
print(paste("FDR_threshold_DE:", FDR_threshold_DE))
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

# load contrasts for pair-wise group comparison
conts <- readr::read_tsv(contrasts, comment="#") |>
  tibble::column_to_rownames("contrast.name")
print("contrasts:")
print(conts)

# set up combinations of clustering parameter
params <- tidyr::expand_grid( 
  selclust = na.omit(c(selected_clustering, selected_CTanno)), 
  cont_name = rownames(conts)
)
print("params:")
print(params)


if(any(is.na(params))) { # skip entirely if cluster setting not specified
  print("Parameter missing! Marker detection skipped")
  de <- list()
} else {
  
  de <- purrr::pmap(params, function(selclust, cont_name) {
    
    name_de_results <- paste0("de_", selclust, "_contrast_", cont_name)
    SingleCellExperiment::colLabels(sce) <- SummarizedExperiment::colData(sce)[,selclust]

    print(paste("Processing", name_de_results))
    
    if(length(list.files(file.path(outdir, selclust), pattern=name_de_results)) > 0 | !selclust %in% colnames(SummarizedExperiment::colData(sce)) | nlevels(factor(SingleCellExperiment::colLabels(sce))) < 2) {
      
      print(paste0(name_de_results, " output already exists or is missing in sce object or has got just one level. Differential gene expression analysis for this setting is skipped"))
      de.results <- NULL
      
    } else {
      
    if (!dir.exists(file.path(outdir, selclust))) {dir.create(file.path(outdir, selclust), recursive=T) }
      
    cont <- conts[cont_name, "contrast"]
    mmatrix  <- conts[cont_name, "mmatrix"]
    groupvar <- tail(all.vars(as.formula(mmatrix)), 1)
    factors   <- gsub("(^\\s+|\\s+$)", "", unlist(strsplit(cont,"\\W")))
    factors   <- factors[factors != ""]
    if(length(factors) != 2) {
      warning(paste(cont,"cannot deal with designs other than pairwise comparisons!"))
      return(NA)
    }

    # creating pseudo-bulk samples from clusters
    summed <- scuttle::aggregateAcrossCells(sce, 
                                            id=SummarizedExperiment::colData(sce)[,c(selclust, "sample")],
                                            use.assay.type = "counts") 
    
    summed.filt <- summed[,SummarizedExperiment::colData(summed)[,groupvar] %in% factors] # use only samples of current contrast
    if(!is.na(minClusterSize)) {
      summed.filt <- summed.filt[,summed.filt$ncells >= minClusterSize]
      if((dim(summed.filt)[2]==0)) {stop("no samples left after filtering for minClusterSize")}
    }
    # store number of cells per sample-cluster combination
    write.table(SummarizedExperiment::colData(summed.filt)[,c("sample", groupvar, selclust, "ncells")], 
                     file.path(outdir, selclust, paste0(name_de_results, "_overview_cells.txt")), sep="\t", row.names = F, quote = F)
    
    # tidy the variables used in the formula (make vector and drop unused levels):
    vars <- all.vars(as.formula(mmatrix))
    for (v in vars) {
      if (is.character(SummarizedExperiment::colData(summed.filt)[[v]])) SummarizedExperiment::colData(summed.filt)[[v]] <- factor(SummarizedExperiment::colData(summed.filt)[[v]])  # characters -> factors
      if (is.factor(SummarizedExperiment::colData(summed.filt)[[v]]))    SummarizedExperiment::colData(summed.filt)[[v]] <- droplevels(SummarizedExperiment::colData(summed.filt)[[v]]) # drop unused levels
    }
    SummarizedExperiment::colData(summed.filt)[,groupvar] <- relevel(SummarizedExperiment::colData(summed.filt)[,groupvar], factors[2]) # relevel to make factors[2] reference level in contrast

    # For each label, pseudoBulkDGE applies abundance filtering using filterByExpr prior to further analysis.
    # Genes that are filtered out will still show up in the DataFrame for that label, but with all statistics set to NA.
    de.results <- scran::pseudoBulkDGE(summed.filt, 
                                       label=SummarizedExperiment::colData(summed.filt)[,selclust],
                                       design=as.formula(mmatrix),
                                       coef=paste0(groupvar,factors[1]), # as factors[2] is reference level
                                       condition=SummarizedExperiment::colData(summed.filt)[,groupvar],
                                       lfc = lfc_threshold_DE, # log-fold change threshold to use
                                       assay.type = "counts"
    )

    # information for log file
    print(paste(selclust, "cluster succeeded", paste0(names(de.results), collapse=", ")))
    failed_cluster <- S4Vectors::metadata(de.results)$failed
    if(length(failed_cluster)>0) {
      print(paste(selclust, "cluster failed:", paste0(failed_cluster, collapse=", ")))
    }
    if(length(de.results)==0) {return(NULL)}
    
    # create summary table for contrast
    FDR_df <- sapply(de.results, function(res) res$FDR) |> # columns = cluster
      as.data.frame() |>
      tibble::rownames_to_column("Gene") |>
      tidyr::pivot_longer(-Gene, names_to = "Cluster", values_to = "FDR")
    
    logFC_df <- sapply(de.results, function(res) res$logFC) |> # columns = cluster
      as.data.frame() |>
      tibble::rownames_to_column("Gene") |>
      tidyr::pivot_longer(-Gene, names_to = "Cluster", values_to = "logFC")
    
    combined_df <- dplyr::left_join(FDR_df, logFC_df, by = c("Gene", "Cluster"))
    
    summary_df <- combined_df |>
      dplyr::group_by(Cluster) |>
      dplyr::summarize(
        Tested = sum(!is.na(FDR)),
        NotSig = sum(!is.na(FDR) & FDR >= FDR_threshold_DE),
        Up     = sum(!is.na(FDR) & FDR < FDR_threshold_DE & logFC > 0),
        Down   = sum(!is.na(FDR) & FDR < FDR_threshold_DE & logFC < 0)
      ) |> 
      readr::write_tsv(file.path(outdir, selclust, paste0(name_de_results, "_overview_sign_genes.txt")))
    

    # process results per cluster
    purrr::map(names(de.results), function(c) {
        
        # create biological coefficient of variation plot per cluster
        png(file=file.path(outdir, selclust, paste0(name_de_results, "_plotBCV_cluster", c, ".png")), width=7, height=4, units="in", res=300) 
          edgeR::plotBCV(S4Vectors::metadata(de.results[[c]])$y)
        dev.off()
        
        # write DE table per cluster
        cur.results <- de.results[[c]] |>
          as.data.frame() |>
          tibble::rownames_to_column("feature_id") |>
          dplyr::filter(!is.na(PValue)) |>
          dplyr::arrange(PValue) |>
          dplyr::mutate(symbol=SummarizedExperiment::rowData(sce)$feature_symbol[match(feature_id, rownames(sce))]) |>
          dplyr::relocate(feature_id, symbol) |>
          readr::write_tsv(file.path(outdir, selclust, paste0(name_de_results, "_cluster", c, ".txt")))
        
        ## GO enrichment
        if(!is.na(org) && org %in% c("human", "mouse")) {
          print(paste("GO enrichment analysis (BP) with top", top_genes_for_GO, "genes for", name_de_results, "cluster", c)) 
          
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
          
          cur.results$entrez_ids <- AnnotationDbi::mapIds(orgdb, keys=gsub("\\..*$", "", cur.results$feature_id), 
                                                                  column="ENTREZID", keytype="ENSEMBL")
          de_top <- cur.results |> 
            dplyr::filter(FDR < FDR_threshold_DE) |>
            dplyr::slice(1:min(top_genes_for_GO, nrow(cur.results)))

          if (nrow(de_top) > 0) {
            print(paste(sum(is.na(de_top$entrez_ids)), "genes skipped because no Entrez IDs found:", paste(de_top$symbol[is.na(de_top$entrez_ids)], collapse=", ")))
            
            go_out <- limma::goana(unique(na.omit(de_top$entrez_ids)), species=species,
                                   universe=unique(na.omit(cur.results$entrez_ids)))
            
            go_out_filt <- go_out |> # Only keeping BP terms that are not overly general and which are significantly enriched.
              dplyr::filter(Ont=="BP" & N<=500 & P.DE < p_threshold_GO) |> # p-value for over-representation of the GO term in the set.
              tibble::rownames_to_column("GOID") |>
              dplyr::arrange(P.DE) |>
              dplyr::mutate(P.DE = signif(P.DE, digits=4)) |>
              readr::write_tsv(file.path(outdir, selclust, paste0("GOenrichment_top", top_genes_for_GO, "_", name_de_results, "_cluster", c, ".txt")))
          } else {print(paste("GO enrichment skipped for", name_de_results, "cluster", c, "because no significantly diff expressed genes"))}
        }
      })
    } 
    return(de.results)
  })
}


#############################
# save the sessionInformation and R image (no update of sce object necessary)
print("store results")
writeLines(capture.output(sessionInfo()),file.path(outdir, "de_bioc_session_info.txt"))
save(args, params, de, file=file.path(outdir,"de_bioc.RData"))
