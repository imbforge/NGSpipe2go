# smallRNA-Seq pipeline
Here we provide the tools to perform a single read smallRNA-Seq analysis including raw data quality control, adapter and unique molecular identifier (UMI) trimming, read filtering and differential expression (DE) analysis. A DE analysis is run on all detected genes as well as only based on genes of a selected type of smallRNA (see variable ESSENTIAL_SMALLRNA below). Optionally, it is run only using counts of mature miRNAs. The pipeline requires zipped fastq-files (.fastq.gz) as input.


## Pipeline Workflow
Specify the desired analysis details for your data in file *smallrnaseq.essential.vars.groovy* (see below) and run the pipeline *smallrnaseq.pipeline.groovy* similarly as described [here](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/README.md). A markdown file *smallRNAreport.Rmd* will be generated in the output reports folder after running the pipeline. Additionally, a markdown file *smallRNAreport.TYPE.Rmd* will be generated in the same folder. It reports on DE results of just your selected type of smallRNA (see variable ESSENTIAL_SMALLRNA below). In case you choose to analyze mature miRNAs, a markdown file *smallRNAreport.miRNAmature.Rmd will be generated, reporting DE results of mature miRNAs. Subsequently, all markdown files can be converted to HTML reports using the *knitr* R-package.


### The pipelines includes
- quality control of raw data with FastQC
- adapter trimming with Cutadapt
- quality filtering using the FASTX-Toolkit
- duplicate removal (optional)
- UMI trimming on both ends of the reads using seqtk
- quality control of trimmed reads with FastQC
- competitive mapping to rRNAs and other contaminants (optional)
- read mapping to the reference genome using Bowtie
- read quantification with featureCounts (Subread)
- count extraction of a selected type of smallRNA (see variable ESSENTIAL_SMALLRNA below)
- generation of bigWig tracks for visualization of alignment with deepTools
- RNA class representation
- illustration of sample relatedness with MDS plots and heatmaps
- Differential Expression Analysis for depicted group comparisons with DESeq2


### Pipeline parameter settings
- targets.txt: tab-separated TXT file giving information about the analyzed samples. Required columns are:
  - sample: sample identifier for use in plots and tables
  - file: read counts file name (a unique sub-string of the file name is sufficient, this sub-string is grepped against the count file names produced by the pipeline) 
  - group: variable for sample grouping (e.g. condition)
  - replicate: replicate number of samples belonging to the same group
- contrasts.txt: indicate intended group comparisons for differential expression analysis. Give one contrast per line. Required columns are:, 
  - contrast.name: user chosen name of the contrast/comparison (e.g. *KO.vs.WT*, if targets.txt contains groups *KO* and *WT*)
  - contrast: description of the contrast/comparison (e.g. *(KO-WT)*)
  - mmatrix: column names of *targets.txt* and their relation to be used for the comparison by DESeq2 (format as required by DESeq2, e.g. *~group*)
- smallrnaseq.essential.vars.groovy: essential parameters describing the experiment including: 
  - ESSENTIAL_PROJECT: your project path and folder name
  - ESSENTIAL_THREADS: number of threads for parallel tasks
  - ESSENTIAL_SAMPLE_PREFIX: common sample prefix, will be removed from sample names for plotting
  - ESSENTIAL_SMALLRNA: type of smallRNA to be analyzed (e.g. "miRNA")
  - ESSENTIAL_BOWTIE_REF: path to and name of Bowtie index
  - ESSENTIAL_GENESGTF: properly formatted GTF annotation file including *different* gene and transcript IDs (e.g. from Ensembl or Gencode)
  - ESSENTIAL_MIRNAGFF: GFF annotation file as provided by miRBase
  - ESSENTIAL_FEATURETYPE: attribute name describing the RNA type as given in ESSENTIAL_GENESGTF, Gencode uses gene_type; Ensembl uses gene_biotype
  - ESSENTIAL_PAIRED: are paired-end reads provided?, should be set to "no" in case of typical smallRNA datasets
  - ESSENTIAL_STRANDED: strandness of library (no|yes|reverse)
  - ESSENTIAL_ORG: organism name
  - ESSENTIAL_BAMCOVERAGE: parameters used for bigWig coverage track generation 
  - ESSENTIAL_READLENGTH: read length of library  (incl. insert, UMIs, adapter)
  - ESSENTIAL_MINADAPTEROVERLAP: minimal overlap of read and adapter (used for adapter trimming)
  - ESSENTIAL_MINREADLENGTH: minimal read length after adapter trimming (incl. UMIs) 
  - ESSENTIAL_UMI_LENGTH: total length of all UMIs
  - ESSENTIAL_ADAPTER_SEQUENCE: adapter sequence to be trimmed by Cutadapt
  - ESSENTIAL_MINIMAL_QUAL: minimally required base quality, reads with any base quality below this value will be removed during quality filtering
  - ESSENTIAL_DESEQ2_FDR: FDR significance cutoff used by DESeq2
  - ESSENTIAL_DESEQ2_FC: optional fold change cutoff used for the DESeq2 model
  - ESSENTIAL_FASTQSCREEN_PERC: contaminant filter, if a contaminant is consuming at least this percentage of reads in at least one sample, contaminant will be shown in report 
  - ESSENTIAL_FASTQSCREEN_GENOME: genome names and indexes used by FastQScreen for the competitive mapping
  - additional (more specialized) parameters can be given in the header files of individual pipeline modules 


## Programs required
- Bowtie
- Cutadapt (optional)
- deepTools
- FastQC
- FastQScreen (optional)
- FASTX-Toolkit
- R packages, e.g. DESeq2, GenomicRanges, ggplot2, rtracklayer
- Samtools
- seqtk
- Subread


## Details worth mentioning

- As e.g. miRNAs are highly similar to each other, reads often map to multiple highly similar or even identical locations and, thus, multi-mapping reads need to be used when analyzing smallRNA sequencing data. One of the best hits is randomly picked. When running featureCounts, parameter -M is set. This can be changed by setting variable *count_multimapping* to *false* in module header files *subread.header* and/or *subread_mirnamature.header*.
- Various types of smallRNA are annotated at overlapping locations, e.g. miRNAs can be found as part of snoRNAs or overlapping intronic regions of protein coding genes. In order to not miss any smallRNA, which (partially) shares its location with another smallRNA or gene, ambiguous reads are used for quantification. When running featureCounts on all genes or on the selected TYPE of smallRNA, parameter -O is set. This can be changed by setting variable *count_ambiguous* to *false* in module header file *subread.header*. When ambiguous reads are used to quantify genes, they are counted on all overlapping genes and, thus, are counted more than once. This leads to "over-counting" of these reads and the total quantification might be larger than the total read count. If this is unfavorable, parameter *--fraction* can be used when running featureCounts in order to only count *1/n* for each of *n* overlapping genes.


