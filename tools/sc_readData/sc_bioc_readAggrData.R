#####################################
##
## What: sc_bioc_readAggrData.R
## Who : Frank Rühle, Patrick Hüther
## When: 17.04.2025
##
## Script to import aggregated single cell count data to bioconductor singlecellexpression object for various assays.
##
## Args:
## -----
## pipeline_root    # pipeline directory   
## outdir           # output directory of this module   
## resultsdir       # result directory of project
## seqtype          # sequencing type     
## targets          # user defined sample information in targets.txt 
## gtf.file         # gtf file path
## aggr_data_dir    # path to aggregated count data
## run_demux        # if sample demultiplexing was applied
## demux_out        # output dir of demux module
## demuxCluster_out # output dir of assignSouporcellCluster module            
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
pipeline_root      <- parseArgs(args,"pipeline_root=") 
outdir             <- parseArgs(args,"outdir=") 
resultsdir         <- parseArgs(args,"res=")   
seqtype            <- parseArgs(args,"seqtype=")   
gtf.file           <- parseArgs(args,"gtf=","")       
ftargets           <- parseArgs(args,"targets=","targets.txt")     
aggr_data_dir      <- parseArgs(args,"aggr_data_dir=", "aggr") 
run_demux          <- parseArgs(args,"run_demux=") 
demux_out          <- parseArgs(args,"demux_out=") 
demuxcluster_out   <- parseArgs(args,"demuxCluster_out=") 

# load R environment
env.path <- file.path(getwd(), pipeline_root, "tools/sc_readData", "bioc_3.16.lock")
print(paste("load renv:", env.path))
renv::use(lockfile=env.path)

# check parameter
print(paste("pipeline_root:", pipeline_root))
print(paste("outdir:", outdir))
print(paste("resultsdir:", resultsdir))
print(paste("seqtype:", seqtype))
print(paste("gtf.file:", gtf.file))
print(paste("ftargets:", ftargets))
print(paste("aggr_data_dir:", aggr_data_dir))
print(paste("run_demux:", run_demux))
print(paste("demux_out:", demux_out))
print(paste("demuxcluster_out:", demuxcluster_out))


# check if the tenX setting is used to analyse the RNA-Seq part of a multiome assay and adjust aggr_data_dir accordingly.
if (!dir.exists(aggr_data_dir)) {
  if(seqtype=="tenX" & dir.exists(gsub("count/", "", aggr_data_dir))) {
    aggr_data_dir <- gsub("count/", "", aggr_data_dir)
    warning(paste("The specified aggr_data_dir wasn't found. Assuming you have RNA-Seq data of a 10X Multiome assay, the aggr_data_dir is changed to:", aggr_data_dir))
  } else {
    stop(paste("The specified aggr_data_dir wasn't found:", aggr_data_dir))
  }
}

# load gtf and targets.txt
gtf <- rtracklayer::import.gff(gtf.file, format="gtf", feature.type="exon") |> tidyr::as_tibble()
targets_pools <- read.delim(ftargets, sep="\t", comment.char = "#") # load targets.txt (data per fastq file).
targets_pools$file <- gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", targets_pools$file) # remove file suffix required for 10X cellranger


# load barcodes, features and count data
print(paste("Load barcodes, features and count data for sequencing type:", seqtype))
switch(seqtype,
  tenX = {
    barcodes <- readr::read_tsv(file.path(aggr_data_dir, "barcodes.tsv.gz"),col_names=c("cell_id")) |> tidyr::separate_wider_delim(cell_id,names=c("barcode","GEMwell"),delim="-",cols_remove=F)
    features <- readr::read_tsv(file.path(aggr_data_dir, "features.tsv.gz"),col_names=c("gene_id","gene_name","feature_type"))
    counts   <- Matrix::readMM(file=file.path(aggr_data_dir, "matrix.mtx.gz"))
    aggrcsv <- read.delim(file.path(resultsdir, "aggr.csv"), sep=",") # aggr.csv contains the order of GEM wells
    # Remark: Change in aggr CSV column headers: For Cell Ranger v6.0 (and later) and Loupe Browser v5.1.0 (and later), 
    # the header should be sample_id,molecule_h5. For prior software versions, it should be library_id,molecule_h5.
  },
  tenXmultiome = {
    barcodes <- readr::read_tsv(file.path(aggr_data_dir, "barcodes.tsv.gz"),col_names=c("cell_id")) |> tidyr::separate_wider_delim(cell_id,names=c("barcode","GEMwell"),delim="-",cols_remove=F)
    features <- readr::read_tsv(file.path(aggr_data_dir, "features.tsv.gz"),col_names=c("gene_id","gene_name","feature_type"))
    counts   <- Matrix::readMM(file=file.path(aggr_data_dir, "matrix.mtx.gz"))
    aggrcsv <- read.delim(file.path(resultsdir, "aggr.csv"), sep=",") # aggr.csv contains the order of GEM wells
  },
  ParseBio = {
    barcodes <- readr::read_csv(file.path(aggr_data_dir, "cell_metadata.csv")) |> dplyr::rename(cell_id = 1) 
    features <- readr::read_csv(file.path(aggr_data_dir, "all_genes.csv"))
    counts <- Matrix::readMM(file=file.path(aggr_data_dir, "count_matrix.mtx")) |> t() # for WT assays, count_matrix.mtx contains cells in rows, genes in columns!
  },
  ScaleBio = { # one folder per sample
    barcodes <- purrr::map_dfr(fs::dir_ls(aggr_data_dir, type = "directory"), ~{
      readr::read_csv(file.path(.x, "barcodes.tsv"), col_names = FALSE) |>
      dplyr::rename(cell_id = X1) |>
      dplyr::mutate(file = stringr::str_remove(basename(.x), "\\.filtered\\.matrix"))
    })
    feature_list <- purrr::map(fs::dir_ls(aggr_data_dir, type = "directory"), ~{
      readr::read_tsv(file.path(.x, "features.tsv"), col_names = c("gene_id","gene_name","feature_type"))
      })
      # Check if all elements in the list are equal to the first. Remark: identical() would also compare object pointer attributes.
    if (purrr::every(feature_list, ~ isTRUE(all.equal(.x, feature_list[[1]])))) {
      features <- feature_list[[1]] # If identical, return the first as a dataframe
    } else {
      stop("Not all features.tsv files are identical.")
    }
    counts <- fs::dir_ls(aggr_data_dir, type = "directory") |> 
      purrr::map(~ Matrix::readMM(file.path(.x, "matrix.mtx"))) |>
      purrr::reduce(cbind)
    },
  SmartSeq = {
    counts_list <- fs::dir_ls(aggr_data_dir, type = "file", regexp = "\\.readcounts\\.tsv$") |>
      purrr::map(function(path) {
        df <- readr::read_tsv(path, col_names = c("gene_id", "count"))
        colname <- basename(path) |> stringr::str_remove("\\.readcounts\\.tsv$")
        df |> dplyr::rename(!!colname := count)
      })
    if (!purrr::every(counts_list, ~ isTRUE(all.equal(.x$gene_id, counts_list[[1]]$gene_id)))) {
      stop("Gene ids are not identical for all countfiles.")
    }
    # Reduce by joining on 'gene_id'
    counts <- purrr::reduce(counts_list, dplyr::left_join, by = "gene_id") |>
      tibble::column_to_rownames("gene_id")
    barcodes <- data.frame(file=colnames(counts)) |> 
      dplyr::mutate(cell_id=file) # naming for downstream compatibility. For SmartSeq, each fastq file refers to exactly one cell.
    features <- data.frame(gene_id=rownames(counts))
  },
  stop(c("Don't find seqtype:", seqtype))   
)

# merge features with gtf
id_cols <- intersect(c("gene_id", "gene_name"), colnames(features)) # e.g. for SmartSeq we have the gene_id column only
features <- features |>
  dplyr::left_join(gtf,by=id_cols,multiple="any") |>
  dplyr::add_count(gene_name) |> # keep track of duplicate gene_names
  dplyr::mutate(gene_name=dplyr::case_when(is.na(gene_name) ~ gene_id, n > 1 ~ glue::glue("{gene_name}_{gene_id}"), .default = gene_name)) # replace missing values with gene_id; append gene_id to duplicate names


## create targets file and demultiplex barcodes if needed and merge with targets_pools
print(paste("Create cellwise targets file from cell barcodes and targets.txt"))
# load demux information if available
if(run_demux != "" & file.exists(demux_out) ) {
  print(paste("Sample multiplexing has been applied:", run_demux))
  if(!seqtype %in% c("tenX", "tenXmultiome")) {stop("Demultiplexing implemented only for 10X assays!")}
  
  switch(run_demux,
    demux_HTO= { # add HTO information if cell hashing applied
      
      print(paste("Reading demux files:", run_demux))
      targets <- readr::read_tsv(fs::dir_ls(path=demux_out, recurse=2, glob="*Seurat/demux.txt"),id="file") |> # obtain pruned file name from Souporcell dir
        dplyr::filter(HTO_classification.global=="Singlet") |> # retain only singlets
        dplyr::mutate(file=basename(stringr::str_remove(file,"/Seurat/demux.txt"))) |>
        dplyr::rename(barcode = cell_id) |>
        dplyr::mutate(GEMwell=as.character(match(file, gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", aggrcsv[,1])))) |>
        dplyr::inner_join(barcodes,by=c("barcode", "GEMwell")) |> # barcode alone is not unique, it may occur in multiple GEMwells.
        dplyr::mutate(name_HTO=basename(stringr::str_remove(HTO_maxID,"-.*$"))) |>
        dplyr::left_join(targets_pools, by=c("file", "name_HTO")) |>
        dplyr::rename(assignment = name_HTO) |> # rename for compatibility with demux_gt
        dplyr::select(-barcode, -orig.ident, -hash.ID, -starts_with("RNA_snn_res")) |>
        dplyr::select(cell_id, sample, file, group, replicate, GEMwell, assignment, samplename_suffix_HTO, file_HTO, everything())
      
     },     
    demux_GT_noAssignment= ,  # if sample de-multiplexing by genetic variance is applied (demux_GT_noAssignment or demux_GT)
    demux_GT= {
      
      if(run_demux == "demux_GT" & file.exists(demuxcluster_out)) { # if assignSouporcellCluster module used
        
        print("The cell cluster from each file refer to the same individuals. Harmonize the cluster assignment with assignSouporcellCluster output.")
        ## Prepare map of cluster re-assignments provided by assignSouporcellCluster module.
        # Note that even after harmonization of clusters (individuals) across files, this harmonized_assignment id
        # does still not match the 'replicate' column in targets.txt (which is arbitrary). The information to relate the cell cluster 
        # back to 'replicate' column in targets.txt is lost due to mixing of samples 
        # We need further external SNP data of those individuals for comparison, to relate the individuals back to targets.txt. 
        # Souporcell offers a script to do that if the SNP data is available: see 'Correlating Cluster to Donor Reference SNP Genotypes (optional)'
        # at: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Souporcell.html
        
        # select one (pruned) file name to use as reference. The clusters of the other files will be aligned to this one.
        global_reference <- targets_pools$file[1]
        
        # get assignSouporcellCluster results
        clusterassignFilenames <- list.files(demuxcluster_out, pattern = "\\.txt$", full.names=T)
    
        assignCl_map <- purrr::map_dfr(clusterassignFilenames, function(filepath) {
          filename <- basename(filepath)
          files <- stringr::str_remove(filename, "\\.txt$") |> stringr::str_split("_vs_", simplify = TRUE)
          file1 <- files[1]
          file2 <- files[2]
          
          # Skip if neither side is the global reference
          if (!(global_reference %in% files)) return(NULL)
          
          # Determine reference and target
          if (file1 == global_reference) {
            ref_file <- file1
            other_file <- file2
            from_col <- "experiment2_cluster"
            to_col <- "experiment1_cluster"
          } else {
            ref_file <- file2
            other_file <- file1
            from_col <- "experiment1_cluster"
            to_col <- "experiment2_cluster"
          }
          
          # Read clustercount and mapping
          clustercount <- read.delim(filepath, header = FALSE, skip = 3, nrows = 2) |>
            dplyr::pull(V1) |>
            stringr::str_extract("\\d+$") |>
            as.numeric() |>
            min()
          
          mapping <- read.delim(filepath, sep = "\t", skip = 5, nrows = clustercount) |>
            dplyr::select(from = !!from_col, to = !!to_col) |>
            dplyr::mutate(
              file = other_file, 
              assignment = as.character(from),
              harmonized_assignment = as.character(to) 
            )|>
            dplyr::select(file, assignment, harmonized_assignment)
          
          # Add identity map if we haven't already seen the global reference
          identity <- dplyr::tibble(
            file = ref_file,
            assignment = unique(as.character(mapping$harmonized_assignment)),
            harmonized_assignment = assignment
          )
          
          dplyr::bind_rows(mapping, identity)
        }) |>
          dplyr::distinct() # remove duplicated lines of global reference, if necessary
      }
        
      # Generate a harmonized_assignment column in targets_pools (from targets.txt) as sample index per file. 
      # In targets_pools this is arbitrary, but it is used for merging with targets to link the cell barcodes with information from targets.txt.
      # I.e. means that during merging, all cells per file with harmonized_assignment '0' in targets (from demuxfiles) get 
      # assigned the target_pools info with harmonized_assignment '0' (from targets.txt).
      # Remark: we don't transform 'harmonized_assignment' in assignCl_map to 'replicate' and then use this column for
      # merging with target_pools. The user defined 'replicate' may be named differently than the cluster names '0', '1', etc.
      # Also, multiple files may refer to the same group meaning in total there are more replicates per group than samples in a file.
      # We may need the replicate column later if external SNP data of individuals is available to explicitly assign the 
      # souporcell cell clusters to specific replicates defined by the user.
      targets_pools <- targets_pools |> 
        dplyr::arrange(file, replicate, sample) |> 
        dplyr::group_by(file) |>
        dplyr::mutate(harmonized_assignment = as.character(0:(dplyr::n_distinct(sample)-1))) |> # enumerate samples per file (start with 0)
        dplyr::ungroup()
        
      ## read cluster files from demux_gt module and combine with assigned clusters and targets_pools (targets.txt info)
      # need aggrcsv to match the file name to the respective GEMwell id.
      print(paste("Reading demux files:", run_demux))
      targets <- readr::read_tsv(fs::dir_ls(path=demux_out, recurse=1, glob="*clusters.tsv"),id="file") |> # obtain pruned file name from Souporcell dir
        dplyr::filter(status=="singlet") |> # retain only singlets
        dplyr::mutate(file=basename(stringr::str_remove(file,"/clusters.tsv"))) |>
        dplyr::mutate(barcode=stringr::str_remove(barcode,"-.*$")) |>
        dplyr::mutate(GEMwell=as.character(match(file, gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", aggrcsv[,1])))) |>
        dplyr::inner_join(barcodes,by=c("barcode","GEMwell")) # barcode alone is not unique, it may occur in multiple GEMwells.
        
        if (exists("assignCl_map")) {
          targets <- targets |>
            dplyr::left_join(assignCl_map,by=c("file", "assignment")) |> # add harmonized assignment if available
            dplyr::left_join(targets_pools,by=c("file", "harmonized_assignment")) 
        } else { # if demux_gt ran without subsequent cluster assignment, merge cluster arbitrarily to samples via assignment column
          targets <- targets |>
            dplyr::left_join(targets_pools |> dplyr::rename(assignment = harmonized_assignment), by=c("file", "assignment")) 
        }
      
      targets <- targets |>
        dplyr::select(-barcode, -starts_with("cluster")) |>
        dplyr::select(cell_id, sample, file, group, replicate, any_of(c("assignment", "harmonized_assignment")), GEMwell, everything())
    },
    stop(c("Don't find demux method:", run_demux))   
  )

} else { # no demux

  print(paste("No sample multiplexing has been applied."))
  switch(seqtype, # prepare targets when no demultiplexing applied
         tenX = {
           targets <- barcodes |>
            dplyr::mutate(
              file = aggrcsv[as.numeric(GEMwell), 1],
              file = gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", file)
              ) |>
              dplyr::left_join(targets_pools,by=c("file")) |> # without demux, 'file' in targets.txt should be unique!
              dplyr::select(cell_id, sample, file, group, replicate, GEMwell, everything())
         },
         tenXmultiome = {
           targets <- barcodes |>
             dplyr::mutate(
               file = aggrcsv[as.numeric(GEMwell), 1],
               file = gsub("(_S\\d{1,3}$)|(_S\\d{1,3}_L\\d{3}_R\\d_\\d{3}$)", "", file)
             ) |>
             dplyr::left_join(targets_pools,by=c("file")) |> # without demux, 'file' in targets.txt should be unique
             dplyr::select(cell_id, sample, file, group, replicate, GEMwell, everything())
         },
         ParseBio = {
           targets_pools <- targets_pools |> # summarize multiple rows per sample induced by distribution of samples across fastq files
            dplyr::group_by(sample) |>
            dplyr::summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) 
          targets <- barcodes |> # merge targets with targets_pools  
            dplyr::left_join(targets_pools, by="sample") |>
            dplyr::relocate(cell_id, colnames(targets_pools)) |>
            dplyr::relocate(ends_with("_wind"), ends_with("_well"), file, .after = last_col()) 
         },
         ScaleBio = {
          targets <- barcodes |> # merge targets with targets_pools  
             dplyr::left_join(targets_pools, by="file")
         },
         SmartSeq = {
           targets <- barcodes |> # merge targets with targets_pools  
             dplyr::left_join(targets_pools, by="file")
         }
  )
}


# write targets to disk
readr::write_tsv(targets, file=file.path(outdir, "targets_cellwise.txt"))
print(paste("dim barcodes:", paste0(dim(barcodes), collapse=", ")))
print(paste("dim targets:", paste0(dim(targets), collapse=", ")))

# prepare counts matrix
colnames(counts) <- barcodes$cell_id
rownames(counts) <- features$gene_id

# subset cell barcodes for targets and features for gex type (if given)
counts <- counts[,targets$cell_id]
if("feature_type" %in% colnames(features)) {
  counts <- counts[features$feature_type=="Gene Expression",]
  features <- features[features$feature_type=="Gene Expression",] # filter features as well for adding to rowData(sce) below
} 

# create SingleCellExperiment object
print("Prepare sce object.")
sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(counts)),colData=targets)

# add gene symbols from annotation file to SingleCellExperiment object
SummarizedExperiment::rowData(sce)$feature_id      <- features$gene_id
SummarizedExperiment::rowData(sce)$feature_symbol  <- features$gene_name
SummarizedExperiment::rowData(sce)$feature_chrom   <- features$seqnames

# spill to disk
print("Spill to disk.")
readr::write_rds(sce, file=file.path(resultsdir, "sce.RDS"))
readr::write_rds(gtf, file=file.path(resultsdir, "gtf.rds"))
save(pipeline_root, outdir, resultsdir, seqtype, gtf.file, ftargets, aggr_data_dir, run_demux, demux_out, demuxcluster_out, 
     targets, targets_pools,
     file=paste0(outdir,"/sc_bioc_readAggrData.RData"))
