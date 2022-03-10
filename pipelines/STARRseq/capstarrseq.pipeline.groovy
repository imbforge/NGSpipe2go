PIPELINE="CapSTARRseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"

load PIPELINE_ROOT + "/pipelines/STARRseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/STARRseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bam2bw.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/bamqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/filterchromosomes.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/markdups.header"
load PIPELINE_ROOT + "/modules/NGS/rmdups.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.header"
load PIPELINE_ROOT + "/modules/RNAseq/deseq2.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread.header"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.header"
load PIPELINE_ROOT + "/modules/RNAseq/tpm.header"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.header"
load PIPELINE_ROOT + "/modules/STARRseq/capstarrseqfoldchange.header"
load PIPELINE_ROOT + "/modules/STARRseq/starrpeaker_procbam.header"
load PIPELINE_ROOT + "/modules/STARRseq/starrpeaker_callpeak.header"
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
    
        "%.bam" * [MarkDups + BAMindexer + dupRadar] + collect_bams + "%.bam" * // branch withduplicates 
            [    (RUN_IN_PAIRED_END_MODE ? InsertSize.using(subdir:"withduplicates") : dontrun.using(module: "InsertSize")),
                 bamCoverage.using(subdir:"withduplicates"),
                 subread_count.using(subdir:"withduplicates") + filter2htseq.using(subdir:"withduplicates") +
                 tpm.using(subdir:"withduplicates"),
                 (RUN_STARRPEAKER ? STARRPeaker_procBam.using(subdir:"withduplicates") +
                                    STARRPeaker_callPeak.using(subdir:"withduplicates") : dontrun.using(module:"starrpeaker"))
            ], 
        "%.bam" * [RmDups + BAMindexer] + collect_bams + "%.bam" *  // branch deduplicated
            [    (RUN_IN_PAIRED_END_MODE ? InsertSize.using(subdir:"deduplicated") : dontrun.using(module: "InsertSize")),
                 bamCoverage.using(subdir:"deduplicated"),
                 subread_count.using(subdir:"deduplicated") + filter2htseq.using(subdir:"deduplicated") +
                 tpm.using(subdir:"deduplicated"),
                 (RUN_STARRPEAKER ? STARRPeaker_procBam.using(subdir:"deduplicated") +
                                    STARRPeaker_callPeak.using(subdir:"deduplicated") : dontrun.using(module:"starrpeaker"))
            ]
        
    ] + // end parallel branches
    
    [
        (CAPSTARRSEQ_DIFFEXP ? DE_DESeq2.using(subdir:"withduplicates") : CapSTARRseq_FoldChange.using(subdir:"withduplicates")),
        (CAPSTARRSEQ_DIFFEXP ? DE_DESeq2.using(subdir:"deduplicated") : CapSTARRseq_FoldChange.using(subdir:"deduplicated"))
    ] +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module: "trackhub")) +
    collectToolVersions + MultiQC + shinyReports
}


