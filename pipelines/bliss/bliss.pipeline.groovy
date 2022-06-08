PIPELINE="bliss"
PIPELINE_VERSION="1.0.0"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/bliss/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/bliss/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/bliss/pattern_filtering.header"
load PIPELINE_ROOT + "/modules/bliss/bwa.header"
load PIPELINE_ROOT + "/modules/bliss/count_breaks.header"
load PIPELINE_ROOT + "/modules/bliss/count_breaks_strandless.header"
load PIPELINE_ROOT + "/modules/bliss/collect_stats.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"


//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }
collect_bams = { forward inputs.bam }

Bpipe.run {
	(RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
    FastQC,
    pattern_filtering +
    bwa +
    count_breaks +
    count_breaks_strandless +
    collect_stats
  ] +
  collectToolVersions +
  MultiQC
}
