PIPELINE="tenXatac"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

load PIPELINE_ROOT + "/modules/scRNAseq/cellrangeratac_count.header"
load PIPELINE_ROOT + "/modules/scRNAseq/cellrangeratac_aggr.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/markdups2.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"


dontrun = { println "didn't run $module" }

Bpipe.run { 
    "%.fastq.gz" * [ FastQC + FastqScreen +
      (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) ] + 
      "%_S*_L*_R*_001.fastq.gz" * [
         cellrangeratac_count + [
            bamCoverage,
            InsertSize
         ]
    ] + 
    cellrangeratac_aggr +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + MultiQC + shinyReports
}

