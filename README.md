![IMB-logo](resources/IMB_logo.png)

# NGSpipe2go #

NGSpipe2go is a software framework to facilitate the development, customization and deployment of NGS analysis pipelines developed and utilised at the Bioinformatics Core of IMB in Mainz (https://www.imb.de/). It uses Bpipe as a workflow manager, provides access to commonly used NGS tools and generate comprehensive reports based on R Markdown.

![NGSpipe2go scheme](resources/NGSpipe2go_scheme.png)

## tenxSTARRseq Pipeline ## 


### Subpipelines ###

- tenxSTARRseq: TenX single cell/nucleus RNA-seq processing with STARR mRNAs separated out and processed separately.
- tenxSTARRampliconseq: STARR mRNAs amplified out of TenX single nucleus RNA-seq library and processed on their own.

### Pipeline stages ###



## NGSpipe2go preparations ##

### Copy NGSpipe2go into the project dir ###

NGS projects should be run in a consistant and reproducible way, hence NGSpipe2go asks you to copy all tools into the project folder, which will ensure that you always use the same program versions at a later time point. This can be done either from a local NGSpipe2go copy or using the most recent development version from the GitLab repository. After cloning, cd into the repository and switch to the `STARRseq` branch.

    git clone https://gitlab.rlp.net/imbforge/NGSpipe2go <project_dir>/NGSpipe2go
    cd  <project_dir>/NGSpipe2go
    git checkout tenxSTARRseq


### Create symlinks for the pipeline ###

Go to your <project_dir>, and make symlinks for the pipeline in the main project dir with:

    ln -s NGSpipe2go/pipelines/tenxSTARRseq/* .

### Customise NGSpipe2go to your needs ###

There are two pipeline files:

- *tenxSTARRseq.pipeline.groovy:* For processing the TenX sn/scRNA-seq data (uses cellranger)
- *tenxSTARRampliconseq.pipeline.groovy:* For CapSTARR-seq (STARR-seq with specific target regions enriched with a DNA capture approach)

Adjust the project-specific information in the pipeline dependent files (see pipeline specific README files for detailed information):

- *essential.vars.groovy* specifies the main project variables like project dir and reference genome
- *xxx.header* additional software parameters can be customised in the vars-file accompanying each bpipe module.
- *xxx.pipeline.groovy* describes the steps of the selected pipeline and the location of the respective modules
- *targets.txt* and *contrasts.txt* contain the sample names and the differential group comparisons

Optionally adjust some general pipeline settings defined in the NGSpipe2go ***config*** folder:

- *bpipe.config.groovy*: define workload manager resources (default workload manager is "slurm", if not needed set executor="local")
- *preambles.groovy*: define module preambles if needed (or stay with default preambles)
- *tools.groovy*: define default versions and running environments for all installed pipeline tools, modify accordingly if new tools or tool versions are installed on your system. If you want to use a different tool version for a certain project you can overwrite the default value in the pipeline-specific file *NGSpipe2go/pipelines/<pipeline>/tools.groovy*. Currently, we are using lmod, conda and/or singularity containers as running environments. Conda itself is loaded as lmod module if conda tools are specified. If you want to use conda but haven't installed it as lmod module make otherwise sure that it is available in the path.

## Run the pipeline ##

Copy the input FastQ files in the <project_dir>/rawdata folder.

Load the bpipe module customised for the Slurm job manager (we recommend to use GNU Screen for persistence), e.g.

    screen
    module load bpipe/0.9.9.8.slurm

Start running the relevant pipeline. For the main 10X run:

    bpipe run tenxSTARRseq.pipeline.groovy rawdata/*.fastq.gz

For the STARR mRNA amplification pipeline:

    bpipe run tenxSTARRampliconseq.pipeline.groovy rawdata/*.fastq.gz

## Compile a project report ##

The results of the pipeline modules will be saved in the *./results* folder. The final report Rmd-file is stored in the *./reports* folder and can be edited or customised using a text editor before converting into a HTML report using knitr.
    
R usage: main 10X pipeline

    rmarkdown::render("sc.report.Rmd")

The STARR mRNA amplification pipeline does not have a separate report script.

