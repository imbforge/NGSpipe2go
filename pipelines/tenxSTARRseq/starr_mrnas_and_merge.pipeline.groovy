PIPELINE="starr_mrnas"
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

// STARR mRNA processing pipeline
// Can be used for STARR mRNAs from original 10X sequencing, or for PCR-amplified STARR data
// Run after the 10X endogenous mRNA pipeline ("split_and_tenx_endogenous.pipeline.groovy")
// Merges its counts with the 10X pipeline's, and stores in SingleCellExperiment .RDS file
Bpipe.run { 
    "%.fastq.gz" * [
          AddUMIBarcodeToFastq +
          STAR + BAMindexer +
          subread_count + BAMindexer + umicount
    ] + MergeAmpliconCounts +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + collectBpipeLogs + MultiQC + shinyReports
}

