PIPELINE="breaktag"
PIPELINE_VERSION="1.0.0"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/breaktag/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/breaktag/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/breaktag/pattern_filtering.header"
load PIPELINE_ROOT + "/modules/breaktag/bwa.header"
load PIPELINE_ROOT + "/modules/breaktag/count_breaks.header"
load PIPELINE_ROOT + "/modules/breaktag/collect_stats.header"
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
    collect_stats
  ] +
  collectToolVersions +
  MultiQC
}
