PIPELINE="ChIPseq_pe"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/ChIPseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/ChIPseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/ChIPseq/GREAT.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/blacklist_filter.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2pe.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/diffbind.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/ipstrength.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/macs2.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/pbc.module.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/peak_annotation.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/insertsize.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/markdups.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/rmdups.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/trackhub.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/multiqc.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/ChIPseq/shinyreports_pe.module.groovy"

// standard PE workflow using deduplicated BAM files of MapQ filtered (~ unique) reads
// assumes that duplicated library fragments are unwanted PCR artifacts
// discards multimapping reads as habitually done in most ChIP-seq studies
//
filter_bam = segment {
  [ bowtie2_pe + BAMindexer + BamQC + filbowtie2unique + BAMindexer + RmDups + BAMindexer + [ bamCoverage, InsertSize, ipstrength, macs2 ] ]
}

// alternative PE workflow using the unfiltered BAM files
// may be preferable when studying some types of repetetive regions
// make sure to have ESSENTIAL_DUP="auto" for MACS2 peak calling
dont_filter_bam = segment {
  [ bowtie2_pe + BAMindexer + BamQC + [ MarkDups + BAMindexer, bamCoverage, InsertSize, ipstrength, macs2 ] ]
}

//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%.fastq.gz" * [ FastQC ] + "%.R*.fastq.gz" *
    (RUN_USING_UNFILTERED_BAM ? dont_filter_bam : filter_bam) +
    (RUN_DIFFBIND ? diffbind : dontrun.using(module:"diffbind")) +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    (RUN_PEAK_ANNOTATION ? peak_annotation : dontrun.using(module:"peak_annotation")) +
    MultiQC + collectBpipeLogs + shinyReports
}

