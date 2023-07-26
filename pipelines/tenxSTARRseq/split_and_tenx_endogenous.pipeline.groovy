PIPELINE="tenx_endogenous"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/tenxSTARRseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/tenxSTARRseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/scRNAseq/cellranger_count.header"
load PIPELINE_ROOT + "/modules/scRNAseq/cellranger_aggr.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.header"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.header"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.header"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/tenxSTARRseq/splitmrna.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"


dontrun = { println "didn't run $module" }

// 10X endogenous mRNA processing pipeline
// Starts by splitting FASTQs into endogenous & STARR mRNAs, then runs Cellranger on endogenous mRNAs
// STARR mRNAs to be processed later with the separate "starr_mrnas_and_merge.pipeline.groovy" pipeline
Bpipe.run { 
    "%.fastq.gz" * [ FastQC + FastqScreen +
      (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) ] + 
      "%_R*_001.fastq.gz" * [ SplitmRNA ] +
      "%_S*_L*_R*_001.fastq.gz" * [              // Cell Ranger branch for endogenous mRNA
          cellranger_count + [
             bamCoverage,
             inferexperiment,
             qualimap,
             subread2rnatypes,
             dupRadar,
             geneBodyCov2
          ]
      ] + cellranger_aggr +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + collectBpipeLogs + MultiQC + shinyReports
}
