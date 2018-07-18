MODULE_FOLDER="./NGSpipe2go/modules/"

load MODULE_FOLDER + "ChIPseq/essential.vars.groovy"
load MODULE_FOLDER + "ChIPseq/tool.locations.groovy"
load MODULE_FOLDER + "ChIPseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/fastqc.vars.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"

load MODULE_FOLDER + "ChIPseq/bowtie2pe.vars.groovy"
load MODULE_FOLDER + "ChIPseq/bowtie2pe.module.groovy"

load MODULE_FOLDER + "ChIPseq/filbowtie2unique.vars.groovy"
load MODULE_FOLDER + "ChIPseq/filbowtie2unique.module.groovy"

load MODULE_FOLDER + "NGS/markdups.vars.groovy"
load MODULE_FOLDER + "NGS/markdups.module.groovy"

load MODULE_FOLDER + "NGS/bamindexer.vars.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"

load MODULE_FOLDER + "NGS/insertsize.vars.groovy"
load MODULE_FOLDER + "NGS/insertsize.module.groovy"

load MODULE_FOLDER + "NGS/bamcoverage.vars.groovy"
load MODULE_FOLDER + "NGS/bamcoverage.module.groovy"

load MODULE_FOLDER + "NGS/trackhub_config.vars.groovy"
load MODULE_FOLDER + "NGS/trackhub_config.module.groovy"

load MODULE_FOLDER + "NGS/trackhub.vars.groovy"
load MODULE_FOLDER + "NGS/trackhub.module.groovy"

load MODULE_FOLDER + "ChIPseq/ipstrength.vars.groovy"
load MODULE_FOLDER + "ChIPseq/ipstrength.module.groovy"

load MODULE_FOLDER + "ChIPseq/pbc.module.groovy"

load MODULE_FOLDER + "ChIPseq/macs2.vars.groovy"
load MODULE_FOLDER + "ChIPseq/macs2.module.groovy"

load MODULE_FOLDER + "ChIPseq/blacklist_filter.vars.groovy"
load MODULE_FOLDER + "ChIPseq/blacklist_filter.module.groovy"

load MODULE_FOLDER + "ChIPseq/peak_annotation.vars.groovy"
load MODULE_FOLDER + "ChIPseq/peak_annotation.module.groovy"  

load MODULE_FOLDER + "ChIPseq/GREAT.vars.groovy"
load MODULE_FOLDER + "ChIPseq/GREAT.module.groovy"             

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

load MODULE_FOLDER + "ChIPseq/shinyreports_pe.vars.groovy"
load MODULE_FOLDER + "ChIPseq/shinyreports_pe.module.groovy"

//MAIN PIPELINE TASK
run {
    "%.fastq.gz" * [ FastQC ] + "%.R*.fastq.gz" * 
    [ bowtie2_pe + BAMindexer + filbowtie2unique + BAMindexer + 
    [ MarkDups + BAMindexer, bamCoverage, InsertSize, pbc, ipstrength, macs2 ] ] + 
    peak_annotation +
//    trackhub_config + trackhub + 
    collectBpipeLogs + shinyReports
}

