![IMB-logo](resources/IMB_logo.png)

# NGSpipe2go #

A set of NGS data analysis tools and pipelines developed and utilised at the Institute of Molecular Biology gGmbH in Mainz (https://www.imb.de/). Currently NGSpipe2go contains bpipe-based pipelines for QC, processing, analysis and visualisation of RNA-seq, ChIP-seq and DNA-seq data.

## Prerequisites ##
### RNA-seq pipeline ###
#### Programs required ####
- FastQC
- STAR
- Samtools
- Subread package
- Picard tools
- BED tools
- UCSC tool set
- RSeQC
- EdgeR
- DEseq2
- dupRadar (*)

(*) is provided via another project from imbforge

#### Files required ####
- targets.txt (+) (*)
- contrasts.txt (+) (*)
- raw reads or mapped data

(+) files are needed to run the EdgeR and DEseq2 modules.

(*) examples provided within this project

### ChIP-seq pipeline ###
#### Programs required ####
- FastQC
- Bowtie 1
- Samtools
- BED tools
- Picard tools
- UCSC tool set
- encodeChIPqc (*)
- MACS2

(*) is provided via another project from imbforge

#### Files required ####
- targets.txt (*)
- raw reads or mapped data

(*) example provided within this project

## Preparations to run ##

### Get NGSpipe2go to your project ###
NGS projects should be run in a consistant and reproducible way, hence NGSpipe2go asks you to copy all tools into the project folder, which will ensure that you always use the same program versions at a later time point.
This can be done either from a local copy you created before:

    cp -r github/NGSpipe2go/ project_folder/

a freshly loaded version from the GitHub releases

    https://github.com/imbforge/NGSpipe2go/releases

or the most recent development version from the GitHub repository

    git clone https://github.com/imbforge/NGSpipe2go.git project_folder/NGSpipe2go

### Select the NGSpipe2go you run ###

Select a pipeline to run and make it available in the main project.

    ln -s project_folder/NGSpipe2go/pipelines/RNAseq/* project_folder/
    ln -s project_folder/NGSpipe2go/modules/RNAseq/essential.vars.groovy project_folder/
    ln -s project_folder/NGSpipe2go/modules/RNAseq/tool.locations.groovy project_folder/

or 

    ln -s project_folder/NGSpipe2go/pipelines/ChIPseq/* project_folder/
    ln -s project_folder/NGSpipe2go/modules/ChIPseq/essential.vars.groovy project_folder/
    ln -s project_folder/NGSpipe2go/modules/ChIPseq/tool.locations.groovy project_folder/


In case of RNAseq experiments you need to define which files belong to which groups (targets.txt) and which comparisons (contrasts.txt) to perform.

    ln -s project_folder/NGSpipe2go/tools/DE_DESeq2/targets.txt project_folder/
    ln -s project_folder/NGSpipe2go/tools/DE_DESeq2/contrasts.txt project_folder/


### Modify NGSpipe2go to your needs ###

Adjust the information found in the following files:

- Project_folder needs to be entered in essential.vars.groovy (+)
- Folder information in the recently linked rnaseq- or chipseq pipeline file, e.g. rnaseq_v1.2.txt (+)
- Folder information in project_folder/NGSpipe2go/modules/RNAseq/tool.locations or project_folder/NGSpipe2go/modules/ChIPseq/tool.locations (*)
- Compute requirements for the used queueing system according to the project data (Whole genome data might need more compute resources than a low coverage ChIPseq experiment) (*)
- File and sample descriptions in targets.txt (+)
- Comparison descriptions in contrasts.txt (+) 

(*) these steps need to be done once for the setup of NGSpipe2go

(+) these steps need to be repeated for each project to be run

## Run ##

We suggest to put the input files to a folder, e.g. project_folder/rawdata.

To start running the pipeline (tested for bpipe-0.9.8.7)

    bpipe run rnaseq_v1.2.txt rawdata/*.fastq.gz
or

    bpipe run chipseq_v1.2.txt rawdata/*.fastq.gz

## Report ##

The final result of the provided pipelines will be saved in project/reports.
The Rmd file can be adjusted to your needs using your favorite text editor (well except MS word).
The Rmd files can be transformed to easily readable documents using knitr, e.g.:
    
    R usage:
    rmarkdown::render("DEreport.Rmd")

