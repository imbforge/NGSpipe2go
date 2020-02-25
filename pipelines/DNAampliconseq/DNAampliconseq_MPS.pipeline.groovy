PIPELINE="DNAampliconseq_MPS"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./ngspipe2go/"  // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/DNAampliconseq/essential.vars.groovy"

load PIPELINE_ROOT + "/pipelines/DNAampliconseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"


load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/DNAampliconseq/addumibarcodetofastq.module.groovy"
load PIPELINE_ROOT + "/modules/DNAampliconseq/pear.module.groovy"
load PIPELINE_ROOT + "/modules/DNAampliconseq/barcode_count.module.groovy"
load PIPELINE_ROOT + "/modules/DNAampliconseq/MPSprofiling.module.groovy"

load PIPELINE_ROOT + "/modules/NGS/multiqc.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"

// load PIPELINE_ROOT + "/modules/DNAampliconseq/shinyreports.module.groovy"

// Main pipeline task
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%.fastq.gz" * [ FastQC ] +
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
        (RUN_PEAR ? pear : dontrun.using(module:"pear")) +
               AddUMIBarcodeToFastq + barcode_count
    ] +
            (RUN_MPSprofiling ? MPSprofiling : dontrun.using(module:"MPSprofiling")) +
            MultiQC + collectToolVersions + collectBpipeLogs
}
