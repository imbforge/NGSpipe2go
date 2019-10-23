PIPELINE="ChIPseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/ChIPseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/ChIPseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"

load PIPELINE_ROOT + "/modules/ChIPseq/GREAT.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/blacklist_filter.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/diffbind.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/ipstrength.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/macs2.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/pbc.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/peak_annotation.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/phantompeak.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/extend.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/markdups.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/trackhub.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/multiqc.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/shinyreports.module.groovy"

//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
  "%.fastq.gz" * [ FastQC , bowtie_se + BAMindexer + MarkDups + BAMindexer + [ extend + bamCoverage, BamQC , phantompeak , pbc , ipstrength , macs2 ] ] +
  (RUN_PEAK_ANNOTATION ? peak_annotation : dontrun.using(module:"peak_annotation")) +
  (RUN_DIFFBIND ? diffbind : dontrun.using(module:"diffbind")) +
  (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
  MultiQC + collectBpipeLogs + shinyReports
}

