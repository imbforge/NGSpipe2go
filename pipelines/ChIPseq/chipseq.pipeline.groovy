MODULE_FOLDER="./NGSpipe2go/modules/"

load MODULE_FOLDER + "ChIPseq/essential.vars.groovy"
load MODULE_FOLDER + "ChIPseq/tool.locations.groovy"
load MODULE_FOLDER + "ChIPseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/fastqc.vars.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"

load MODULE_FOLDER + "ChIPseq/bowtie1.vars.groovy"
load MODULE_FOLDER + "ChIPseq/bowtie1.module.groovy"

load MODULE_FOLDER + "NGS/bamindexer.vars.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"

load MODULE_FOLDER + "NGS/rmdups.vars.groovy"
load MODULE_FOLDER + "NGS/rmdups.module.groovy"

load MODULE_FOLDER + "NGS/extend.vars.groovy"
load MODULE_FOLDER + "NGS/extend.module.groovy"

load MODULE_FOLDER + "NGS/bam2bw.vars.groovy"
load MODULE_FOLDER + "NGS/bam2bw.module.groovy"

load MODULE_FOLDER + "ChIPseq/ipstrength.vars.groovy"
load MODULE_FOLDER + "ChIPseq/ipstrength.module.groovy"

load MODULE_FOLDER + "ChIPseq/phantompeak.vars.groovy"
load MODULE_FOLDER + "ChIPseq/phantompeak.module.groovy"

load MODULE_FOLDER + "ChIPseq/pbc.module.groovy"

load MODULE_FOLDER + "ChIPseq/macs2.vars.groovy"
load MODULE_FOLDER + "ChIPseq/macs2.module.groovy"

load MODULE_FOLDER + "ChIPseq/peak_annotation.vars.groovy"                                                                                                                    
load MODULE_FOLDER + "ChIPseq/peak_annotation.module.groovy"               

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

load MODULE_FOLDER + "ChIPseq/shinyreports.vars.groovy"
load MODULE_FOLDER + "ChIPseq/shinyreports.module.groovy"

//MAIN PIPELINE TASK
run {
	"%.fastq.gz" * 
	[ FastQC , bowtie_se + BAMindexer + [ extend + bam2bw , phantompeak , pbc , ipstrength , macs2 ] ] + 
	peak_annotation + collectBpipeLogs + shinyReports
}
