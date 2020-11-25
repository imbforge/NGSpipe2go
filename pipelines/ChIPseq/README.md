# ChIP-Seq pipeline
Here we provide the tools to perform paired end or single read ChIP-Seq analysis including raw data quality control, optional adapter trimming, read mapping, peak calling, differential binding analysis and functional annotation. As input files you may use either zipped fastq-files (.fastq.gz) or mapped read data (.bam files). In case of paired end reads, corresponding fastq files should be named using *.R1.fastq.gz* and *.R2.fastq.gz* suffixes. The Pipeline performes 2 branches for bam file processing in parallel:
- filtered branch: remove multimapping reads (MapQ filtered) as habitually done in most ChIP-seq studies and (optionally) removes duplicates from BAM file (assuming that duplicated library fragments are unwanted PCR artifacts. Duplicate removal is not recommended for single end design!).
- unfiltered branch: alternative workflow using the unfiltered BAM files. May be preferable when studying some types of repetitive regions. Make sure to have set ESSENTIAL_DUP="auto" for MACS2 peak calling.


## Pipeline Workflow
All analysis steps are illustrated in the pipeline [flowchart](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1#G1sfhhIib0KGAMbqAvbqYkbM8wFCXXymwB). Specify the desired analysis details for your data in the *essential.vars.groovy* file (see below) and run the pipeline *chipseq.pipeline.groovy* as described [here](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/README.md). A markdown file *ChIPreport.Rmd* will be generated in the output reports folder after running the pipeline. Subsequently, the *ChIPreport.Rmd* file can be converted to a final html report using the *knitr* R-package.


### The pipelines includes
- raw data quality control with FastQC, BamQC and MultiQC.
- adapter trimming with Cutadapt (optional)
- mapping reads or read pairs to the reference genome using bowtie2 (default) or bowtie1.
- filter out multimapping reads from bowtie2 output with samtools (performed in parallel pipeline branch).
- identify and remove duplicate reads with Picard MarkDuplicates (performed in parallel pipeline branch). 
- generation of bigWig tracks for visualisation of alignment with deeptools bamCoverage. For single end design, reads are extended to the average fragment size.
- characterization of insert size using Picard CollectInsertSizeMetrics (for paired end libraries only).
- characterize library complexity by PCR Bottleneck Coefficient using the GenomicAlignments R-package (for single read libraries only). 
- characterize phantom peaks by cross correlation analysis using the spp R-package (for single read libraries only).
- peak calling of IP samples vs. corresponding input controls using MACS2.
- peak annotation using the ChIPseeker R-package (optional).
- differential binding analysis using the diffbind R-package (optional). For this, contrasts of interest must be given in *NGSpipe2go/pipelines/ChIPseq/contrasts_diffbind.txt* (see below).


### Pipeline-specific parameter settings
- targets.txt: tab-separated txt-file giving information about the analysed samples. The following columns are required: 
  - IP: full IP sample name without any filetype suffixes like ".R1.fastq.gz" or ".bam" (sample name itself must not contain dots).
  - IPname: (shortened) IP sample name to be used in plots and tables.
  - INPUT: full input control sample name without any filetype suffixes like ".R1.fastq.gz" or ".bam" (sample name itself must not contain dots).
  - INPUTname: (shortened) input sample name to be used in plots and tables. 
  - group: variable for sample grouping (e.g. by condition).
  - Replicate: number of replicate (needed if differential binding analysis is selected).
  - PeakCaller: name of peak caller (in this pipeline it is macs; needed if differential binding analysis is selected).


- essential.vars.groovy: essential parameter describing the experiment including: 
  - ESSENTIAL_PROJECT: your project folder name.
  - ESSENTIAL_BOWTIE_REF: full path to bowtie2 indexed reference genome (bowtie1 indexed reference genome if bowtie1 is selected as mapper).
  - ESSENTIAL_BOWTIE_GENOME: full path to the reference genome FASTA file.
  - ESSENTIAL_BSGENOME: Bioconductor genome sequence annotation package.
  - ESSENTIAL_TXDB: Bioconductor transcript-related annotation package.
  - ESSENTIAL_ANNODB: Bioconductor genome annotation package.
  - ESSENTIAL_BLACKLIST: path to bed file with problematic 'blacklist regions' to be excluded from analysis (optional).
  - ESSENTIAL_PAIRED: either paired end ("yes") or single read ("no") design.
  - ESSENTIAL_READLEN: read length of library.
  - ESSENTIAL_FRAGLEN: mean length of library inserts.
  - ESSENTIAL_MACS2_BROAD: use broad setting for broad peak calling in MACS2 (default false).
  - ESSENTIAL_PEAK_MINLENGTH: MACS2 minimum peak length. Can be set to the same value as ESSENTIAL_FRAGLEN (MACS2 default is fragment size to filter out smaller peaks). For paired end data it may be beneficial to use a value below the fragment length. Should be increased for filtering if broad option is used.
  - ESSENTIAL_DEDUPLICATION: remove duplicated reads in filtered branch (default false is strongly recommended for single end data!).  
  - ESSENTIAL_THREADS: number of threads for parallel tasks.
  - ESSENTIAL_USE_BOWTIE1: if true use bowtie1 for read mapping, otherwise bowtie2 by default.
  - ESSENTIAL_FULLLIBRARYSIZE: if true (default) use total number of reads in bam for normalization in DiffBind analysis (if false only reads overlapping peaks).
 
- additional (more specialized) parameter can be given in the header files of the individual pipeline modules.

If differential binding analysis is selected it is required additionally:

- contrasts_diffbind.txt: indicate intended group comparisions for differential binding analysis, e.g. *KOvsWT=(KO-WT)* if targets.txt contains the groups *KO* and *WT*. Give 1 contrast per line.


## Programs required
- Bedtools
- Bowtie2 (or Bowtie1)
- Cutadapt
- deepTools
- encodeChIPqc (provided by another project from imbforge)
- FastQC
- MACS2
- MultiQC
- Picard
- R with packages ChIPSeeker, diffbind, GenomicAlignments, spp and genome annotation packages
- Samtools
- UCSC utilities
