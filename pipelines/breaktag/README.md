# (s)BLISS and BreakTag pipeline
Here we provide the tools to perform paired end or single read BreakTag raw data processing. The pipeline is also valid for (s)BLISS data. As input files you may use either gzipped fastq-files (.fastq.gz) or mapped read data (.bam files). In case of paired end reads, corresponding fastq files should be named using *.R1.fastq.gz* and *.R2.fastq.gz* suffixes. 
Some steps of the pipeline are based on the [blissNP pipeline](https://github.com/BiCroLab/blissNP) developed at the BriCo lab for BLISS data.

## The pipelines includes the following steps
- raw data quality control with FastQC, BamQC and MultiQC.
- mapping reads or read pairs to the reference genome using bwa.
- identify and remove duplicate reads with a custom approach which identifies duplicated reads with "close" UMIs (as in the original [blissNP pipeline](https://github.com/BiCroLab/blissNP)).
- generation of bigWig tracks for visualisation of alignment with deeptools bamCoverage. For single end design, reads are extended to the average fragment size.

## Pipeline-specific parameter settings (files you need to setup in order to run the pipeline):
- `targets.txt`: tab-separated txt-file giving information about the analysed samples. The following columns are required
  - name: sample name. Experiment ID found in fastq filename: expID_R1.fastq.gz
  - pattern: UMI+barcode pattern file used in the linker

- essential.vars.groovy: essential parameters describing the experiment
  - ESSENTIAL_PROJECT: root folder of the analysis
  - ESSENTIAL_SAMPLE_PREFIX: sample name prefix to be trimmed in the results
  - ESSENTIAL_THREADS: number of threads for parallel tasks
  - ESSENTIAL_BWA_REF: path to bwa indexed reference genome
  - ESSENTIAL_PAIRED: either paired end ("yes") or single read ("no") design
  - ESSENTIAL_QUALITY: minimum mapping quality desired

## Programs required
- Bedtools
- BWA
- FastQC
- MultiQC
- Samtools
- Several unix standard tools (perl, awk, etc.)
