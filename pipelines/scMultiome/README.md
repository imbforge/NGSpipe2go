# Single cell multiome analysis pipelines

Here we forge the tools to analyze single cell multiome experiments. The analysis workflow for the single cell multiome pipelines are based on the [Seurat](https://satijalab.org/seurat/) and [Signac](https://stuartlab.org/signac/articles/overview.html) package vignettes. 

## Implemented protocols
Provide the name of the respective sequencing type given below as value for the ESSENTIAL_SEQTYPE parameter in essential.vars.groovy.
- tenXmultiome: droplet-based multiome pipeline combining snRNA-Seq and snATAC-Seq data by 10X Genomics assays. For this pipeline, input fastq-files must follow the naming convention (*[Sample Name]_[Library Type]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz*, with *Read_Type* being one of R1, R2, I1, I2 and *Library Type* being one of gex or atac).


## Pipeline Workflow
Specify desired analysis details for your data in the respective *essential.vars.groovy* file (see below) and run the selected pipeline as described [here](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/README.md). The analysis allows further parameter fine-tuning within the header files of the respective analysis modules. Go through the parameters and modify the default settings where appropriate. 

Recommendation: For some of the time consuming modules like clustering it would make sense to have the initial parameter settings (QC, normalization) optimized beforehand. For this, you can use bpipe's *-u* testing functionality: 

    bpipe run -u sc_filter tenXmultiome.pipeline.groovy rawdata/*.fastq.gz
    
The command above will run the pipeline until (excluding) the sc_filter module, which removes cells from the dataset based on the specified QC filter thresholds. You may inspect the QC plots and statistics provided so far before applying the filter thresholds for downstream analysis. After adjustment of QC filter settings specific to your dataset in the respective header file you can re-run the pipeline without *-u* flag. Subsequently, the *scmultiome.report.Rmd* file can be converted to a final html report using the *knitr* R-package.

### The pipelines includes:
- FastQC, FastQScreen and other tools for rawdata quality control
- optional adapter trimming with Cutadapt (not necessary if provided by manufacturer software)
- mapping reads to the reference genome
- generation of bigWig tracks for visualisation of alignment
- read quantification 
- optional sample de-multiplexing if samples were multiplexed per GEM well in 10X experiments
- QC on cell level
- normalization
- clustering
- marker gene detection per cluster
- GO term enrichment analysis on the marker genes defining each cluster
- cell type annotation
- differential expression analysis between sample groups per identified cell cluster using pseudo bulking
- differentially accessible peaks 
- motif enrichment analysis for differentially accessible peaks
- differential motif activity analysis 


### Pipeline parameter settings
- essential.vars.groovy: essential parameter describing the experiment 
  - ESSENTIAL_PROJECT: your project folder name.
  - ESSENTIAL_SEQTYPE: sequencing type "tenXmultiome"
  - ESSENTIAL_SAMPLE_PREFIX: common sample prefix to be removed from output plots.
  - ESSENTIAL_GENOME_REFERENCE: path to reference genome as requested for the respective assay.
  - ESSENTIAL_ORG: UCSC organism name.
  - ESSENTIAL_DB: UCSC assembly version
  - ESSENTIAL_GENESGTF: path to gtf file containing genome annotation.
  - ESSENTIAL_MTGENES: filename with list of gene_ids of mitochondrial genes (give path within ESSENTIAL_PROJECT). If not given (empty string), mitochondrial genes are identified automatically by default.
  - ESSENTIAL_PAIRED: either paired end ("yes") or single read ("no") design ("no" for MARS-Seq and 10X, because R2 in MARS_Seq and R1 in 10X contain UMI and barcodes only).
  - ESSENTIAL_READLENGTH: read length of library.
  - ESSENTIAL_THREADS: number of threads for parallel tasks.
  - RUN_DEMUX: (10X assays only) specify de-multiplexing method if applied. Either "demux_GT" for demultiplexing by genetic variance with subsequent assignment of corresponding individuals across files, "demux_HTO" for cell hashing or empty string for no demultiplexing.
  - ESSENTIAL_CELLTYPE_ANNO: celltype annotation method to use (one or more of "Seurat", "Marker"). The first in the list is used for downstream processing. For scRNAseq pipeline the specification is done in the respective module header files.
  - ESSENTIAL_EXPECTED_CELLS: number of expected cells in experiment (this info can help in detecting valid cells).
  - ESSENTIAL_NUCLEI: (10X only) set to TRUE if assay runs with nuclei instead of cells (TRUE for "tenXmultiome").
- additional (more specialized) parameter can be given in the header files of the individual pipeline modules (see module header files linked in the flowchart for default parameter). 


- targets.txt: tab-separated txt-file giving information about the analyzed samples. The following columns can be provided:
  - sample: (required) shortened sample identifier to be used in plots. 
  - file: (required) file name identifier. Must be a unique substring of the corresponding input sample file name, e.g. full sample name with filetype suffixes removed (the sample name itself must not contain dots). 
  - group: (required) default variable for cell grouping (e.g. by condition). You can provide additional grouping variables optionally.
  - replicate: (required) define sample replicates within groups.
  - wells: (for combinatorial barcoding approach only): specify the wells each sample was loaded per sub-library (fastq file pair). If one sample is loaded to multiple sub-libraries, add additional lines for each sub-library specifying cells of this sample. Wells are specified in blocks, ranges, or individually like this: *A1:C6* specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6. Multiple selections are joined by commas (no spaces), e.g. *A1-A6,B1:D3,C4*.   
  - file_HTO: (if cell hashing is used only) file name identifier of the HTO library corresponding to *file*
  - seq_HTO: (if cell hashing is used only) sequence of hashtag oligo (HTO)
  - name_HTO: (if cell hashing is used only) name of hashtag oligo (HTO)


- contrasts.txt: tab-separated txt-file with information about the sample groups to be compared in pair-wise differential expression analysis per identified cluster. The following columns are required
  - contrast.name: an identifier/name for the groups you want to compare
  - contrast: indicate intended group comparisons for differential expression analysis; the group names must be defined in the targets.txt file. E.g. *(KO-WT)* if targets.txt contains the groups *KO* and *WT*. 


### Remarks
Optional sample de-multiplexing for multiplexed sampels per GEM well in 10X experiments can be done either by cell hashing or by utilizing natural genetic variance. The latter option is for human samples only. In this case, *targets.txt* must contain one unique sample entry for each de-multiplexed sample. If cell hashing is applied, *targets.txt* must contain additional columns specifying the corresponding hashtag oligo (HTO) libraries. De-multiplexing is done using the CITE-seq-Count tool followed by HTODemux function in Seurat R-package. For de-multiplexing samples by natural genetic variance we use the Souporcell tool. If multiplexed sub-samples belong to replicates used in all GEM wells, we can apply another Souporcell process to align the identified cell clusters (representing the multiplexed sub-samples) across all GEM wells, meaning that replicate 1 from one GEM well refers to replicate 1 in the other GEM wells. REMARK: please mind that we cannot refer those sub-samples back to the sample names given in *targets.txt*. For this you would need externally generated genotypes for each sub-sample and then compare them to those of the identified Souporcell clusters.


### Programs required
- FastQC
- FastQScreen
- STAR
- Cellranger (for tenX)
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
- Trajectory analysis (pseudotime): the [monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html) package.
- A [tutorial](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) from Hemberg lab
- Luecken and Theis 2019 [Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/10.15252/msb.20188746)
- How to estimate required [cell numbers](https://satijalab.org/howmanycells)
- [Seurat vignette](https://satijalab.org/seurat/) 
- [Signac vignette](https://stuartlab.org/signac/articles/overview.html)
