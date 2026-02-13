# Single cell RNA-Seq analysis pipelines

Here we forge the tools to analyze single cell RNA-Seq experiments. The analysis workflow for the scRNA-Seq pipeline in *scRNA.pipeline.groovy* is based on the Bioconductor workflows by Amezquita RA, Lun ATL et al. [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA).

## Implemented scRNAseq protocols
The pipeline covers scRNA-Seq assay formats of several providers. Specify the name of the respective sequencing type given below as value for the ESSENTIAL_SEQTYPE parameter in *essential.vars.groovy*.

- tenX: droplet-based scRNA-Seq by 10X Genomics using Chromium Single Cell 3' Reagent Kits. For this pipeline, input fastq-files must follow the bcl2fastq naming convention (*[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz*, with *Read_Type* being one of R1, R2, I1, I2).
- ParseBio: scRNA-Seq based on combinatorial barcoding using Parse Biosciences assays (Evercode WT). For this pipeline, input fastq-files must be names as follows: *[Sample Name].[Read Type].fastq.gz*, with *Read_Type* being one of R1, R2. If you have fastq files from multiple lanes they must be concatenated. Since cells of a sample can be spread across multiple sub-libraries (fastq files), the targets.txt file needs to be adjusted accordingly (see below).
- ScaleBio: scRNA-Seq based on combinatorial barcoding using Scale Biosciences Single Cell RNA Sequencing kits. The initial analysis including mapping and read counting must be performed with ScaleBio's proprietary ScaleRNA software as a separated step outside of the NGSpipe2go pipeline. The result folder path name of this analysis run is then used as input parameter for the pipeline call here. I.e. instead of specifying fastq file names, you just provide the path of the ScaleRNA output folder as input when running the pipeline.
- SmartSeq: Plate-based scRNA-Seq with sorted cells. Each well contains a single cell. Control wells contain either zero ("0c") or ten cells ("10c"). Libraries are generated using the [Smart-seq2 kit](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2639.html). The sequencing reads are demultiplexed by cell barcode resulting in one (pair of) fastq file(s) per cell. For this pipeline, input fastq-files must be names as follows: *[Sample Name].[Read Type].fastq.gz*, with *Read_Type* being one of R1, R2.


## Pipeline Workflow
All analysis steps of the scRNA-Seq pipeline are illustrated in the following [flowchart](../../resources/NGSpipe2go_scRNAseq.drawio.svg). Specify desired analysis details for your data in the respective *essential.vars.groovy* file (see below) and run the scRNA-Seq pipeline as described [here](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/README.md). The pipeline allows further parameter fine-tuning within the header files of the respective analysis modules. Go through the parameters and modify the default settings where appropriate. When the pipeline run is finished, a *sc.report.Rmd* file has been generated in the *reports* sub-directory which can be rendered to a final html report in R (>v4.3). 


### The pipelines includes
- FastQC and FastQScreen modules for rawdata quality control
- optional adapter trimming with Cutadapt module (not necessary if done by manufacturer software)
- mapping reads to the reference genome and read count quantification using the respective manufacturer's software modules. Refer to the header files of the respective modules (e.g. *cellranger_count*, *splitpipe_all*) to adjust for user-specific settings.
- QC on alignment level including gene body coverage RNA class representation and mapping stats.
- generation of bigWig tracks for visualization of alignment
- optional sample de-multiplexing if samples were multiplexed per GEM well in 10X experiments (see additional remarks below).
- QC on cell level including library size, number of expressed genes, proportion of reads mapping to mitochondrial genes and predicted doublet score.
- optional QC filtering for these categories: modify the default filtering thresholds in the header file of the *filtCells_bioc* module. You can choose between applying absolute or relative thresholds. Also provide a threshold for filtering out low abundance genes to reduce complexity of the dataset.
- normalization using the *norm_bioc* module to produce log-transformed normalized expression values for downstream analysis. This module also identifies the top highly variable genes (HVGs) used to calculate the reduced dimension vectors for further data exploration. The header file of this module allows to adjust normalization settings according to your preferences (e.g. spike-in normalization, HVG thresholds, PCA and UMAP settings). Optionally, you can check whether cell-cycle assignment or any additional sample parameter provided in *targets.txt* (like batches, experimental factors) contribute substantially to the heterogeneity of gene expression. If so, the factor may be specified as blocking variable for HVG and marker detection and in differential expression analysis.
- clustering (kmeans / igraph / hclust): you can choose between k-means clustering, graph-based clustering and hierarchical clustering approaches as implemented in the modules *kmeans_bioc*, *igraph_bioc* and *hclust_bioc*. Settings for these algorithms can be customized for your data-set in the respective module header files. You can specify multiple values per parameter to explore multiple clustering settings. The results for each approach and combination of settings is stored in a separate clustering sub-directory in the results. Already existing sub-directories are not overwritten in case you re-run a clustering module with the same parameter settings. Mind that for huge data-sets the clustering can be quite time-consuming. Especially for hierarchical clustering, it may make sense to try a two-step approach to perform hierarchical clustering on the representative centroids obtained by an initial k-means vector quantization step (see settings in respective header file). The cluster modules also provide doublet detection approaches per cluster, helping you to decide whether a certain cluster is made up by a different cell type or rather by doublets from neighboring cluster. If you have decided for one or more clustering settings you want to use for downstream analysis, specify them in the header file of the *collectCl_bioc* module. 
- Gene expression pots: if you are interested in the expression of certain genes across all cells and clusters, you can specify them (symbols or ensemble IDs) in the *collectCl_bioc* as well to obtain gene expression plots.
- marker gene detection per cluster: the *findmarkers_bioc* module identifies the genes that drive separation between clusters via comparing each pair of clusters. The identified top marker genes per cluster are then used for enrichment analysis to detect GO-terms that are relatively active in that cluster. 
- cell-type annotation: marker-based cell-type annotation per cluster is implemented in the *scType_bioc* module. It uses positive and negative marker genes to annotate a cluster to a specific cell-type. The respective marker gene sets to use can be specified in the header file of that module (see default setting for table format). Mind that any cell-type annotation can be only as good or accurate as the provided marker gene sets it is based on! As for the clustering modules, you can specify multiple input settings to be run. If you have decided for one or more annotations you want to use for downstream analysis, specify them in the header file of the *collectCT_bioc* module.
- differential expression analysis: The *de_bioc* module uses a pseudo-bulking approach to test for changes in gene expression per cell cluster or cell-type for each specified contrast. The identified differentially expressed genes per cluster or cell-type are then used for enrichment analysis to detect GO-terms that are over-represented for that contrast. 


## Recommendation for running
Considering that some of the modules like clustering may be very time-consuming for large data-sets, it can be advisable to have the initial parameter settings (QC, normalization) optimized beforehand. For this, you can use bpipe's *-u* testing functionality, e.g.:  

    bpipe run -u filtCells_bioc scRNAseq.pipeline.groovy rawdata/*.fastq.gz
    
The command above will run the pipeline until (excluding) the *filtCells_bioc* module, which removes cells from the dataset based on the specified QC filter thresholds. You may inspect the QC plots and statistics provided so far by the *qc_bioc* module before applying the filter thresholds for downstream analysis. After adjustment of QC filter settings specific to your dataset in the header file of the *filtCells_bioc* module you can re-run the pipeline without *-u* flag. The same may be done for the following step of normalization by the *norm_bioc* module, in case you wish to make use of the setting adjustment options given in its header file. After the normalization step, the data matrices of the single cell object are not modified anymore. Any downstream module will produce separate result files stored in the respective results sub-directory. 


During exploration your data-set, you may want to try out different settings for clustering or other modules downstream of the normalization module. You can do so by deleting the *.RData* file in the respective results sub-directory and rerun the pipeline. Bpipe will then pick-up the analysis from that module. Existing results are not overwritten (the respective routine will be skipped if the specific result sub-directory exists already). For the clustering and the cell-type annotation modules you can specify which results shall be passed to downstream modules by specifying them in the *collectCl_bioc* or *collectCT_bioc* module, respectively. When creating the report, all available result sub-directories will be collected and included in the report.


### Pipeline parameter settings
- essential.vars.groovy: essential parameter describing the experiment 
  - ESSENTIAL_PROJECT: your project folder name.
  - ESSENTIAL_SEQTYPE: sequencing type, one of "tenX", "ParseBio", "ScaleBio", "SmartSeq".
  - ESSENTIAL_SAMPLE_PREFIX: common sample prefix to be removed from output plots.
  - ESSENTIAL_GENOME_REFERENCE: path to reference genome as requested for the respective assay.
  - ESSENTIAL_ORG: UCSC organism name.
  - ESSENTIAL_DB: UCSC assembly version
  - ESSENTIAL_GENESGTF: path to gtf file containing genome annotation.
  - ESSENTIAL_MTGENES: filename with list of gene_ids of mitochondrial genes (give path within ESSENTIAL_PROJECT). If not given (empty string), mitochondrial genes are identified automatically by default.
  - ESSENTIAL_PAIRED: either paired end ("yes") or single read ("no") design ("no" if paired but one read contains UMI and barcodes only).
  - ESSENTIAL_THREADS: number of threads for parallel tasks.
  - RUN_DEMUX: (10X assays only) specify de-multiplexing method if applied. Either "demux_GT" for demultiplexing by genetic variance, "demux_GT_noAssignment" for demultiplexing by genetic variance without subsequent assignment of corresponding individuals across files, "demux_HTO" for cell hashing or empty string for no demultiplexing.
  - ESSENTIAL_EXPECTED_CELLS: number of expected cells in experiment (this info can help in detecting valid cells).
  - ESSENTIAL_NUCLEI: (10X only) set to TRUE if assay runs with nuclei instead of cells.
  - ESSENTIAL_USE_AGGR_DATA: if true use sample data aggregated by manufacturer's software. Otherwise, load individual sample data and aggregate in Seurat (default: true. Other option not implemented yet)
  - RUN_BATCHCORRECT: (not implemented yet) whether to do batch correction. If *true*, specify the batch-variable as a column name in targets file.
- additional (more specialized) parameter can be given in the header files of the individual pipeline modules (see module header files linked in the flowchart for default parameter).

<!-- -->

- targets.txt: tab-separated txt-file giving information about the analyzed samples. The following columns can be provided:
  - sample: (required) shortened sample identifier to be used in plots.
  - file: (required) file name identifier. Must be a unique substring of the corresponding input sample file name, e.g. full sample name with filetype suffixes removed (the sample name itself must not contain dots). If cells of multiple samples are mixed in one raw data file (as e.g. in ParseBio sub-libraries), you need to add lines for each sample contained in this file. 
  - group: (required) default variable for cell grouping (e.g. by condition). You can provide additional grouping variables optionally.
  - replicate: (required) define sample replicates within groups.
  - wells: (for combinatorial barcoding approach only): specify the wells each sample was loaded per sub-library (fastq file pair).
    - ParseBio: If one sample is loaded to multiple sub-libraries, add additional lines for each sub-library specifying cells of this sample. Wells are specified row-wise in blocks, ranges, or individually like this: *A1:C6* specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6. *A1-B6* specifies a range as [start]-[end]; A1-A12, B1-6. Multiple selections are joined by commas (no spaces), e.g. *A1-A6,B1:D3,C4*.
    - ScaleBio: wells column not needed here but had to be specified when running ScaleBio's ScaleRNA software (in the respective samples.csv file). Mind that for ScaleBio ranges are read in column-wise order, e.g. 1A-2C, refers to 1A-1H (all of column 1) plus 2A-2C.
  - plate: (for SmartSeq only) plate ID (number).
  - row: (for SmartSeq only) plate row (letter).
  - col: (for SmartSeq only) plate column (number).
  - cells: (for SmartSeq only) number of cells per well (one of "0c", "1c", "10c" with "0c" and "10c" being control wells containing 0 or 10 cells, respectively).
  
<!-- -->

- contrasts.txt: tab-separated txt-file with information about the sample groups to be compared in pair-wise differential expression analysis per identified cluster. The following columns are required
  - contrast.name: an identifier/name for the groups you want to compare
  - contrast: indicate intended group comparisons for differential expression analysis; the group names must be defined in the targets.txt file. E.g. *(KO-WT)* if targets.txt contains the groups *KO* and *WT*.
  - mmatrix: the design formula which needs to be used to perform the differential analysis. Different formulas may be given for each contrast. E.g. ~group is sufficient if you do not need to correct for additional factors. The values you use in the design formula must be defined in the targets.txt file.


### Additional remarks
- Genome reference: mind that you need to provide genome references specific for each assay in *essential.vars.groovy*. E.g. for ParseBio, you have to create an indexed reference genome using [split-pipe](support.parsebiosciences.com) in mkref mode. Since split-pipe modifies the chromosome names by adding a genome name prefix (e.g. 'mm39_chr1'), mind to provide a compatible gtf file for the pipeline as well. Otherwise, you would need to remove the alignment QC modules (qualimap, geneBodyCov2, subread2rnatypes) from *scRNA.pipeline.groovy* to prevent them from failing. You can modify the chromosome names in your original gtf file like this (here: add *'mm39_'* prefix to chromosome names: *awk 'BEGIN{FS=OFS="\t"} $0 ~ /^#/ {print; next} {$1="mm39_"$1; print}' input.gtf > output.gtf*
- Sample multiplexing: Optional sample multiplexing per GEM well in 10X experiments can be done either by cell hashing or by utilizing the natural genetic variance between individuals (the latter for human samples only). 
  - Genetic variance: Samples coming from different individuals are pooled and sequenced as one library in one GEM well. We use the [Souporcell tool](https://github.com/wheaton5/souporcell) for subsequent de-multiplexing using the natural genetic variance between these individuals. In this case, *targets.txt* must contain a unique sample entry for each pooled sub-sample per input sample file name (i.e. as many rows per fastq file-pair name as there are sub-samples mixed in this pool). This is crucial because Souporcell will divide the cells of one file into that many clusters representing the mixed samples (it will not auto-detect the number of samples). REMARK: mind that we cannot refer those sub-samples back to the unique sample names given in *targets.txt*, because this information is lost when samples are mixed. For this you would need externally generated genotypes for each sub-sample and then compare them to those of the identified Souporcell clusters. However, if multiplexed samples belong to individuals used in all GEM wells (the same individual but with different treatment), we can apply another Souporcell process to align the identified cell clusters (representing the de-multiplexed sub-samples) across all GEM wells, meaning that replicate 1 from one GEM well refers to replicate 1 in the other GEM wells.
  - [Cell hashing](https://pubmed.ncbi.nlm.nih.gov/30567574/): cells belonging to one sample are labeled with barcoded antibodies before pooling. These barcodes are sequenced as separate library and are then used to demultiplex the pooled samples. This is done by using the [CITE-seq-Count](https://hoohm.github.io/CITE-seq-Count/Running-the-script/) tool followed by the HTODemux function of the [Seurat R-package](https://satijalab.org/seurat/articles/hashing_vignette.html). The *targets.txt* must contain additional columns specifying the corresponding hashtag oligo (HTO) libraries. 
    - file_HTO: file name identifier of the HTO library corresponding to *file*
    - seq_HTO: sequence of hashtag oligo (HTO)
    - name_HTO: name of hashtag oligo (HTO)



### Programs required
- FastQC
- FastQScreen
- STAR
- Cellranger (for tenX)
- spipe (for ParseBio)
- ScaleRNA (for ScaleBio)
- Souporcell (for de-multiplexing by natural genetic variance)
- CITE-seq-Count (for cell hashing)
- Samtools
- Bedtools
- Subread
- Picard
- UCSC utilities
- RSeQC
- UMI-tools
- R

## Resources
- [Orchestrating Single-Cell Analysis with Bioconductor](http://bioconductor.org/books/release/OSCA/)
- A [tutorial](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) from Hemberg lab
- Luecken and Theis 2019 [Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/10.15252/msb.20188746)
- How to estimate required [cell numbers](https://satijalab.org/howmanycells)
- Trajectory analysis (pseudotime): the [monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html) package.
- [Seurat vignette](https://satijalab.org/seurat/) 
- [Signac vignette](https://stuartlab.org/signac/articles/overview.html)
