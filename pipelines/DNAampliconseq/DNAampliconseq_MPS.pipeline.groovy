PIPELINE="DNAampliconseq_MPS"
PIPELINE_VERSION="1.2"
PIPELINE_ROOT="./NGSpipe2go/"  // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/DNAampliconseq/essential.vars.groovy"

load PIPELINE_ROOT + "/pipelines/DNAampliconseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"

load PIPELINE_ROOT + "/modules/DNAampliconseq/pear.header"
load PIPELINE_ROOT + "/modules/DNAampliconseq/addumibarcodetofastq.header"
load PIPELINE_ROOT + "/modules/DNAampliconseq/barcode_count.header"
load PIPELINE_ROOT + "/modules/DNAampliconseq/MPSprofiling.header"

load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/DNAampliconseq/shinyreports.header"

// Main pipeline task
dontrun = { println "didn't run $module" }

Bpipe.run {
   (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
        FastQC + FastqScreen +
        (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) +
        (RUN_PEAR ? pear : dontrun.using(module:"pear")) +
               AddUMIBarcodeToFastq + barcode_count
    ] +
            (RUN_MPSprofiling ? MPSprofiling : dontrun.using(module:"MPSprofiling")) +
            collectToolVersions + MultiQC + shinyReports
}
