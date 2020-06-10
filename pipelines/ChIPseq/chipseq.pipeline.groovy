PIPELINE="ChIPseq"
PIPELINE_VERSION="1.1"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/ChIPseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/ChIPseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/ChIPseq/GREAT.header"
load PIPELINE_ROOT + "/modules/ChIPseq/blacklist_filter.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/diffbind.header"
load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.header"
load PIPELINE_ROOT + "/modules/ChIPseq/ipstrength.header"
load PIPELINE_ROOT + "/modules/ChIPseq/macs2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/pbc.header"
load PIPELINE_ROOT + "/modules/ChIPseq/peak_annotation.header"
load PIPELINE_ROOT + "/modules/ChIPseq/phantompeak.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/extend.header"
load PIPELINE_ROOT + "/modules/NGS/bamqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/markdups.header"
load PIPELINE_ROOT + "/modules/NGS/rmdups.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/shinyreports.header"

// quality control modules specific to paired end (pe) or single end (se) design
qc_paired = segment {
	[
		bamCoverage,
		InsertSize
	]
}

qc_single = segment {
	[
		extend + bamCoverage, 
		phantompeak
	]
}

//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%.fastq.gz" * [ FastQC ] + 
	(RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
		(ESSENTIAL_USE_BOWTIE1 ? bowtie1 : bowtie2) + BAMindexer + BamQC +            
                (RUN_USING_UNFILTERED_BAM ? (MarkDups + BAMindexer + pbc) : (filbowtie2unique + BAMindexer + RmDups + BAMindexer)) + 
		[
	    		(RUN_IN_PAIRED_END_MODE ? qc_paired : qc_single), 
                	 ipstrength, 
		         macs2
		] 
	] + 
        (RUN_PEAK_ANNOTATION ? peak_annotation : dontrun.using(module:"peak_annotation")) + 
        (RUN_DIFFBIND ? diffbind : dontrun.using(module:"diffbind")) +
        (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
        MultiQC + collectToolVersions + collectBpipeLogs + shinyReports
}

