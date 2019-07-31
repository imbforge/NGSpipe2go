MODULE_FOLDER="./NGSpipe2go/modules/"

load MODULE_FOLDER + "ChIPseq/essential.vars.groovy"
load MODULE_FOLDER + "ChIPseq/tool.locations.groovy"
load MODULE_FOLDER + "ChIPseq/tool.versions.groovy"

load MODULE_FOLDER + "ChIPseq/GREAT.module.groovy"
load MODULE_FOLDER + "ChIPseq/blacklist_filter.module.groovy"
load MODULE_FOLDER + "ChIPseq/bowtie1.module.groovy"
load MODULE_FOLDER + "ChIPseq/diffbind.module.groovy"
load MODULE_FOLDER + "ChIPseq/ipstrength.module.groovy"
load MODULE_FOLDER + "ChIPseq/macs2.module.groovy"
load MODULE_FOLDER + "ChIPseq/pbc.module.groovy"
load MODULE_FOLDER + "ChIPseq/peak_annotation.module.groovy"
load MODULE_FOLDER + "ChIPseq/phantompeak.module.groovy"
load MODULE_FOLDER + "ChIPseq/shinyreports.module.groovy"
load MODULE_FOLDER + "NGS/bamcoverage.module.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"
load MODULE_FOLDER + "NGS/bamqc.module.groovy"
load MODULE_FOLDER + "NGS/extend.module.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"
load MODULE_FOLDER + "NGS/markdups.module.groovy"
load MODULE_FOLDER + "NGS/rmdups.module.groovy"
load MODULE_FOLDER + "NGS/trackhub.module.groovy"
load MODULE_FOLDER + "NGS/trackhub_config.module.groovy"
load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

//MAIN PIPELINE TASK
nothing = segment { }
Bpipe.run {
  "%.fastq.gz" * [ FastQC , bowtie_se + BAMindexer + MarkDups + BAMindexer + [ extend + bamCoverage, BamQC , phantompeak , pbc , ipstrength , macs2 ] ] +
  (RUN_PEAK_ANNOTATION ? peak_annotation : nothing) +
  (RUN_DIFFBIND ? diffbind : nothing) +
  (RUN_TRACKHUB ? trackhub_config + trackhub : nothing) +
  collectBpipeLogs + shinyReports
}

