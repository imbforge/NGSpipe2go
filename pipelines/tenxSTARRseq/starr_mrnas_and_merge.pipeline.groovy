PIPELINE="tenxSTARRampliconseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/tenxSTARRseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/tenxSTARRseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/RNAseq/star.header"
load PIPELINE_ROOT + "/modules/scRNAseq/addumibarcodetofastq.header"
load PIPELINE_ROOT + "/modules/scRNAseq/subread.header"
load PIPELINE_ROOT + "/modules/scRNAseq/umicount.header"
load PIPELINE_ROOT + "/modules/tenxSTARRseq/mergeampliconcounts.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"


dontrun = { println "didn't run $module" }

// Amplicon sequencing pipeline: Like the 10X STARR-seq pipeline, but skips the 10X part
Bpipe.run { 
    "%.fastq.gz" * [
          AddUMIBarcodeToFastq +
          STAR + BAMindexer +
          subread_count + BAMindexer + umicount
    ] + MergeAmpliconCounts +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + collectBpipeLogs + MultiQC + shinyReports
}

