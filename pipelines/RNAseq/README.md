# RNA-Seq pipeline
Here we provide the tools to perform paired end or single read RNA-Seq analysis including raw data quality control, differential expression (DE) analysis and functional annotation. As input files you may use either zipped fastq-files (.fastq.gz) or mapped read data (.bam files). In case of paired end reads, corresponding fastq files should be named using *.R1.fastq.gz* and *.R2.fastq.gz* suffixes.


## Pipeline Workflow
All analysis steps are illustrated in the pipeline [flowchart](https://viewer.diagrams.net/?tags=%7B%7D&highlight=0000ff&edit=_blank&layers=1&nav=1&title=NGSpipe2go_RNAseq_pipeline.html#Uhttps%3A%2F%2Fdrive.google.com%2Fuc%3Fid%3D1wEeBp1znDO3MP_vDTiNXzys8kfabNfw2%26export%3Ddownload). Specify hthe desired analysis details for your data in the *essential.vars.groovy* file (see below) and run the pipeline *rnaseq.pipeline.groovy* as described [here](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/README.md). A markdown file *DEreport.Rmd* will be generated in the output reports folder after running the pipeline. Subsequently, the *DEreport.Rmd* file can be converted to a final html report using the *knitr* R-package.


### The pipelines includes
- quality control of rawdata with FastQC and MultiQC
- adapter trimming with Cutadapt (optional)
- Competitive mapping to rRNAs and other contaminants (optional)
- Read mapping to the reference genome using STAR
- generation of bigWig tracks for visualisation of alignment with deeptools
- Characterization of insert size for paired-end libraries
- Read quantification with featureCounts (Subread) 
- Library complexity assessment with dupRadar
- RNA class representation
- Check for strand specificity
- Visualization of gene body coverage
- Illustration of sample relatedness with MDS plots and heatmaps
- Differential Expression Analysis for depicted group comparisons with DESeq2
- Differential splicing analysis for depicted group comparisons with rMATS
- Enrichment analysis for DE results with clusterProfiler and ReactomePA
- Additional DE analysis including multimapped reads


### Pipeline parameter settings
- targets.txt: tab-separated txt-file giving information about the analysed samples. The following columns are required, but any number of additional columns can be added, for e.g. time point, sex, treatment, clone.
  - sample: sample identifier for use in plots and and tables
  - file: read counts file name (a unique sub-string of the file name is sufficient, this sub-string is grepped against the count file names produced by the pipeline) 
  - group: variable for sample grouping (e.g. by condition)
  - replicate: replicate number of samples belonging to the same group
- contrasts.txt: tab-separated txt-file with information about the samples to be compared. The following columns are required
  - contrast.name: an identifier/name for the groups you want to compare
  - contrast: indicate intended group comparisons for differential expression analysis; the group names must be defined in the targets.txt file.  E.g. *(KO-WT)* if targets.txt contains the groups *KO* and *WT*.
  - mmatrix: the design formula which needs to be used to perform the differential analysis. More than one formala maybe given in multiple lines of this file. E.g. *~group* is sufficient if you do not need to correct for additional factors. If you need to additionally correct for, say, sex information *~group+sex* should be used. The values you use in the design formula must be defined in the targets.txt file; in this case *sex* information should be present in the targets.txt file.  
- essential.vars.groovy: essential parameter describing the experiment including: 
  - ESSENTIAL_PROJECT: your project folder name
  - ESSENTIAL_STAR_REF: path to STAR indexed reference genome
  - ESSENTIAL_GENESGTF: properly formatted GTF annotation file having *different* gene and trancript IDs (e.g. from Ensembl or Gencode)
  - ESSENTIAL_PAIRED: either paired end ("yes") or single read ("no") design
  - ESSENTIAL_STRANDED: strandness of library (no|yes|reverse)
  - ESSENTIAL_ORG: UCSC organism name
  - ESSENTIAL_READLENGTH: read length of library
  - ESSENTIAL_THREADS: number of threads for parallel tasks
  - ESSENTIAL_ADAPTER_SEQUENCE: adapter sequence to trim with Cutadapt (optional)
  - ESSENTIAL_BASEQUALCUTOFF: base quality threshold to trim low-quality ends from reads with Cutadapt. If *ESSENTIAL_NEXTSEQTRIM* is true, qualities of terminal G bases are ignored. To switch off base quality trimming in Cutadapt entirely, set *ESSENTIAL_BASEQUALCUTOFF=0* and *ESSENTIAL_NEXTSEQTRIM=false*.
  - ESSENTIAL_NEXTSEQTRIM: most Illumina instruments use a two-color chemistry like the NextSeq (exceptions: MiSeq, HiSeq). This option accounts for terminal high quality G bases incorporated by faulty dark cycles during base quality trimming with Cutadapt.

- additional (more specialized) parameter can be given in the header-files of the individual pipeline modules 


## Programs required
- Bedtools
- Cutadapt (optional)
- DEseq2
- deeptools
- dupRadar (provided by another project from imbforge)
- FastQC
- FastQScreen (optional)
- MultiQC
- Picard
- R packages DESeq2, clusterProfiler, ReactomePA
- RSeQC
- Samtools
- STAR
- rMATS
- Subread
- UCSC utilities
