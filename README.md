![IMB-logo](resources/IMB_logo.png)

# NGSpipe2go #

NGSpipe2go is a software framework to facilitate the development, customization and deployment of NGS analysis pipelines developed and utilised at the Bioinformatics Core of IMB in Mainz (https://www.imb.de/). It uses Bpipe as a workflow manager, provides access to commonly used NGS tools and generate comprehensive reports based on R Markdown.

![NGSpipe2go scheme](resources/NGSpipe2go_scheme.png)

## Available Pipelines ## 

- [ChIP-Seq](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/pipelines/ChIPseq/README.md)
- [DNA-Seq](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/pipelines/DNAseq/README.md)
- [RNA-Seq](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/pipelines/RNAseq/README.md)
- RNA-Seq for variant calling
- smallRNA-Seq
- [single cell RNA-Seq](https://gitlab.rlp.net/imbforge/NGSpipe2go/-/blob/master/pipelines/scRNAseq/README.md)

## NGSpipe2go preparations ##

### Copy NGSpipe2go into the project dir ###

NGS projects should be run in a consistent and reproducible way, hence NGSpipe2go asks you to copy all tools into the project folder, which will ensure that you always use the same program versions at a later time point. This can be done either from a local NGSpipe2go copy or by using the most recent version from the Gitlab repository

    git clone https://gitlab.rlp.net/imbforge/NGSpipe2go <project_dir>/NGSpipe2go

### Choose one of the pipelines ###

Select a pipeline to run and make symlinks in the main project dir, e.g. for RNA-seq projects

    ln -s NGSpipe2go/pipelines/RNAseq/* .

or for ChIP-seq projects

    ln -s NGSpipe2go/pipelines/ChIPseq/* .

### Customise NGSpipe2go to your needs ###

Adjust the project-specific information in the pipeline dependent files (see pipeline specific README files for detailed information):

- *essential.vars.groovy* specifies the main project variables like project dir and reference genome
- *xxx.vars.groovy* additional software parameters can be customised in the vars-file accompanying each bpipe module.
- *xxx.pipeline.groovy* describes the steps of the selected pipeline and the location of the respective modules
- *targets.txt* and *contrasts.txt* contain the sample names and the differential group comparisons

Adjust general pipeline settings defined in the NGSpipe2go ***config*** folder:

- *bpipe.config.groovy*: define workload manager resources (default workload manager is "slurm", if not needed set executor="local"). Important: external users need to modify the default Slurm queues ("bcfshort" and "bcflong") to e.g. "short" and "long" according to their Slurm environment. 
- *preambles.groovy*: define module preambles if needed (or stay with default preambles)
- *tools.groovy*: define default versions and running environments for all installed pipeline tools, modify accordingly if new tools or tool versions are installed on your system. If you want to use a different tool version for a certain project you can overwrite the default value in the pipeline-specific file *NGSpipe2go/pipelines/<pipeline>/tools.groovy*.

## Run a pipeline ##

Copy the input FastQ files in the <project_dir>/rawdata folder.

Load the bpipe module customised for the Slurm job manager (we recommend to use GNU Screen for persistence), e.g.

    screen
    module load bpipe/0.9.9.8.slurm

Start running the pipeline of choice, e.g.

    bpipe run rnaseq.pipeline.groovy rawdata/*.fastq.gz

or

    bpipe run chipseq.pipeline.groovy rawdata/*.fastq.gz    

## Compile a project report ##

The results of the pipeline modules will be saved in the *./results* folder. The final report Rmd-file is stored in the *./reports* folder and can be edited or customised using a text editor before converting into a HTML report using knitr
    
    R usage:
    rmarkdown::render("DEreport.Rmd")
    or
    rmarkdown::render("ChIPreport.Rmd")
