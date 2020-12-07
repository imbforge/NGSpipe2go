![IMB-logo](resources/IMB_logo.png)

# NGSpipe2go #

An opinionated framework for building pipelines. It comprises set of NGS data analysis tools and pipelines developed and utilised at the Institute of Molecular Biology gGmbH in Mainz (https://www.imb.de/).

![NGSpipe2go scheme](resources/NGSpipe2go_scheme.png)

## DNAampliconseq_MPS pipeline

This pipeline performs multiplexed protein stability (MPS) profiling of DNA amplicon-seq data as described in the [publication](https://www.sciencedirect.com/science/article/pii/S1097276518302363) of Kats I, Khmelinskii A, Kschonsak M, Huber F, Knieß RA, Bartosik A, Knop M (2018). *Mapping Degradation Signals and Pathways in a Eukaryotic N-terminome.* Mol Cell. 2018 May 3;70(3):488-501.e5. doi: 10.1016/j.molcel.2018.03.033. It is designed as alternative or addition to the [CombinatorialProfiler](https://github.com/ilia-kats/CombinatorialProfiler) tool allowing for flexible amplicon design and UMI-deduplication. UMI deduplication and extraction of the variable region of the amplicon is implemented using [UMI-tools](https://umi-tools.readthedocs.io/en/latest/index.html). The variable region of interest is referred to as cell barcode within the UMI-tools regular expression describing the amplicon structure. The downstream processing of the count data like background subtraction and calculation of protein stabilty indicies (PSI) is analogous to CombinatorialProfiler. All processing steps are illustrated in the pipeline [flowchart](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1#G1Z44xRviBaLQYeuWo4pYM1tH6W63oMt6h). 

### The pipeline includes

- FastQC for rawdata quality control.
- optional adapter trimming with Cutadapt.
- read pair assembly with PEAR.
- extract UMI-sequences and cell barcodes with UMI-tools extract. Barcodes and UMIs are attached to the read names. Definition of the required regular expressions is explained [here](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction).
- optionally, generate custom whitelist for filtering and correction of cell barcodes according to the observed barcode (and UMI) distribution. 
- deduplication of PCR duplicates by UMI-sequences using UMI-tools.
- optional sample demultiplexing by 2nd barcode sequence (as in "amplicon3" experimental design).
- barcode counting with UMI-tools count_tab (or awk depending on experiment design).
- remove background noise by fitting a Gaussian mixture model to each sample and subtract the mean of the background distribution.
- calculate protein stability indices (with and without background subtraction) for sample fractions coming from the different signal intensity bins.

### Implemented experiment designs

- "amplicon1": Paired end sequencing with overlapping reads for assembly. Amplicon sequence contains cell barcodes to count and UMIs for deduplication. Other elements are optional, but must fit the regular expression given in essential.vars.groovy.
- "amplicon2": as "amplicon1", but without UMIs. Remark: the regular expression still needs an UMI segment for reasons of compatibility with umi_tools, though of length zero: (?P<umi_1>.{0}).
- "amplicon3": As "amplicon2", but with an additional barcode in read2 for sample demultiplexing (i.e. 2 independent cell barcodes but no UMIs). The 2nd barcode is extracted by an additional umi_tools extract step to keep it separated from the 1st barcode. It is copied as sample name into the 2nd column of count file (if both barcodes were extracted in a single umi_tools extract step they are merged in the read name and not separated by "_"). The occurrences of all CB combinations are counted.
- "amplicon4": Paired end sequencing with non-overlapping reads. No pear assembly of read pairs before barcode extraction. Contains two independent cell barcodes (one in each read) but no UMIs. The occurrences of all CB combinations are counted.

### Programs required

- Cutadapt
- FastQC
- MultiQC
- PEAR
- UMI-tools
- R

## NGSpipe2go preparations ##

### Put NGSpipe2go into the project dir ###

NGS projects should be run in a consistant and reproducible way, hence NGSpipe2go asks you to copy all tools into the project folder, which will ensure that you always use the same program versions at a later time point. This can be done either from a local NGSpipe2go copy or using the most recent development version from the GitLab repository. After cloning, cd into the repository and switch to the MPSprofiling branch.

    git clone https://gitlab.rlp.net/imbforge/NGSpipe2go <project_dir>/NGSpipe2go
    cd  <project_dir>/NGSpipe2go
    git checkout MPSprofiling

### Create symlinks for the pipeline ###

Go to your <project_dir> and make symlinks for the pipeline in the main project dir. 

    ln -s NGSpipe2go/pipelines/DNAampliconseq/* .

### Customise NGSpipe2go to your needs ###

- *essential.vars.groovy*: essential parameter describing the experiment 
  - project folder name
  - experiment design (see below)
  - regular expressions for barcode and UMI extraction
  - whitelist options: you may provide user-prepared barcode whitelists to filter the extracted barcodes using the ESSENTIAL_WHITELIST option. If ESSENTIAL_CORRECT_CB is set true, non matching barcodes will be corrected to barcode alternatives given in the whitelist. It also possible to generate a whitelist with likely true barcodes from the data set (as described [here](https://umi-tools.readthedocs.io/en/latest/reference/whitelist.html)) by using the ESSENTIAL_EXTRACT_WHITELIST flag. This whitelist would also contain possible barcode alternatives for correcting if possible. For each barcode extraction the user may use either the ESSENTIAL_WHITELIST or the ESSENTIAL_EXTRACT_WHITELIST option (but not both).
- additional (more specialized) parameter can be given in the header-files of the individual pipeline modules (e.g. Hamming distance for correction of barcodes to whitelist barcode in *addumibarcodetofastq.header*). The flowchart given above links to the header files of each module to inspect the default parameter set.
- *targets.txt*: comma-separated txt-file giving information about the analysed samples. The following columns are required (additional columns ignored): 
  - experiment: experiment name
  - sub_experiment: summarizes those samples, which have been distributed to stability bins and belong together for PSI calculation
  - unique_sample_id: sample identifier. If no sample demultiplexing by a 2nd barcode is necessary, may be same as pruned_file_name
  - pruned_file_name: Unique short form of the input fastq file name (common prefixes and suffixes can be removed). These names are grepped against the count file names to merge the *targets.txt* information to the count data.
  - fraction: fraction of cells assigned to each bin
  - bin: index number of signal intensity bin
  - barcode_demultiplex: (only required if design is "amplicon3") 2nd barcode for demultiplexing samples from count file

## Run the pipeline ##

Copy the input FastQ files into the <project_dir>/rawdata folder.

You may use GNU Screen for persistence when running the pipeline. Load the bpipe module customised for the Slurm job manager, e.g.

    screen
    module load bpipe/0.9.9.3.slurm

Start running the pipeline

    bpipe run DNAampliconseq_MPS.pipeline.groovy rawdata/*.fastq.gz

## Compile a project report ##

The results of the pipeline modules will be saved in the ./results folder. The final report Rmd-file is stored in the ./reports folder and can be edited or customised using a text editor before converting into a HTML report using knitr.
    
    R usage:
    rmarkdown::render("reports/mps.report.Rmd")

## PSI result tables in report file ##

The mps.report.html file contains overview statistics for every step performed in the pipeline. Finally, 4 result tables are given with PSIs calculated with and without background subtraction in two versions:
- PSIs calculated per nucleotides: Protein stability indices are calculated per nucleotide sequence from normalized count data. The statistics given in the table refer to count data of all fractions per sub_experiment and sequence. Column ‘nfractions’ gives the number of fractions (bins) with counts detected, while ‘totalfractions’ contains the total number of fractions for this sub_experiment.
- PSIs calculated per amino acids: Median protein stability indices per amino acid sequence are calculated as median of PSIs from all nucleotide sequences translated into the same amino acid sequence. Additionally, pooled PSIs are calculated directly from normalized count data by pooling counts from all nucleotide sequences translated into the same amino acid sequence. The statistics given in the table refer to the pooled approach. Column ‘nsequences’ gives the number of nucleotide sequences translated into the same amino acid sequence, ‘nfractions’ gives the number of fractions with counts detected for these sequences, while ‘totalfractions’ contains the total number of fractions for this sub_experiment. 

