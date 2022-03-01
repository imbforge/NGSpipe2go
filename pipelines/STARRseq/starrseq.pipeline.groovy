PIPELINE="STARRseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/STARRseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/STARRseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/ChIPseq/GREAT.header"
load PIPELINE_ROOT + "/modules/ChIPseq/blacklist_filter.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/diffbind3.header" 
load PIPELINE_ROOT + "/modules/ChIPseq/diffbind2.header" 
load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.header"
load PIPELINE_ROOT + "/modules/ChIPseq/ipstrength.header"
load PIPELINE_ROOT + "/modules/ChIPseq/macs2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/pbc.header"
load PIPELINE_ROOT + "/modules/ChIPseq/peak_annotation.header"
load PIPELINE_ROOT + "/modules/ChIPseq/phantompeak.header"
load PIPELINE_ROOT + "/modules/ChIPseq/upsetPlot.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/extend.header"
load PIPELINE_ROOT + "/modules/NGS/bamqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/markdups.header"
load PIPELINE_ROOT + "/modules/NGS/rmdups.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/STARRseq/shinyreports.header"


//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }
collect_bams = { forward inputs.bam }

Bpipe.run {
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
        FastQC, 
        (RUN_FASTQSCREEN ? FastqScreen : dontrun.using(module: "FastqScreen")), 
        (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) +
        bowtie2 + BAMindexer + BamQC + filbowtie2unique + BAMindexer ] + collect_bams +   

     [ // parallel branches with and without deduplication

        "%.bam" * [MarkDups + BAMindexer + pbc] + collect_bams + "%.bam" * // branch withduplicates 
            [   // QC specific to paired end (pe) or single end (se) design
                (RUN_IN_PAIRED_END_MODE ? [bamCoverage.using(subdir:"withduplicates"), 
                                           InsertSize.using(subdir:"withduplicates")] : 
                                          [bamCoverage.using(subdir:"withduplicates"), 
                                           phantompeak.using(subdir:"withduplicates")]), 
                ipstrength.using(subdir:"withduplicates"), 
                macs2.using(subdir:"withduplicates")
            ], 
        "%.bam" * [RmDups + BAMindexer] + collect_bams + "%.bam" *  // branch deduplicated
            [   // QC specific to paired end (pe) or single end (se) design
                (RUN_IN_PAIRED_END_MODE ? [bamCoverage.using(subdir:"deduplicated"), 
                                           InsertSize.using(subdir:"deduplicated")] : 
                                          [bamCoverage.using(subdir:"deduplicated"), 
                                           phantompeak.using(subdir:"deduplicated")]), 
                ipstrength.using(subdir:"deduplicated"), 
                macs2.using(subdir:"deduplicated")  
            ]

    ] + // end parallel branches

    [ blacklist_filter.using(subdir:"withduplicates") +
     (RUN_PEAK_ANNOTATION ? peak_annotation.using(subdir:"withduplicates") : dontrun.using(module:"peak_annotation")) +
     (RUN_UPSETPLOT ? upsetPlot.using(subdir:"withduplicates") : dontrun.using(module:"upsetPlot")) + 
     (RUN_DIFFBIND ? (ESSENTIAL_DIFFBIND_VERSION >= 3 ? diffbind3.using(subdir:"withduplicates") : diffbind2.using(subdir:"withduplicates")) : dontrun.using(module:"diffbind")) +
     (RUN_ENRICHMENT ? GREAT.using(subdir:"withduplicates") : dontrun.using(module:"GREAT")),

      blacklist_filter.using(subdir:"deduplicated") +
     (RUN_PEAK_ANNOTATION ? peak_annotation.using(subdir:"deduplicated") : dontrun.using(module:"peak_annotation")) +
     (RUN_UPSETPLOT ? upsetPlot.using(subdir:"deduplicated") : dontrun.using(module:"upsetPlot")) +
     (RUN_DIFFBIND ? (ESSENTIAL_DIFFBIND_VERSION >= 3 ? diffbind3.using(subdir:"deduplicated") : diffbind2.using(subdir:"deduplicated")) : dontrun.using(module:"diffbind")) +
     (RUN_ENRICHMENT ? GREAT.using(subdir:"deduplicated") : dontrun.using(module:"GREAT"))
    ] +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + MultiQC + 
    shinyReports
}

