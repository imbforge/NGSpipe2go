PIPELINE="smallRNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/smallrnaseq.essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bam2bw.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/trim_umis.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter_stats.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastqscreen.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/bowtie1.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/cutadapt.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/cutadapt_stats.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup_stats.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mapping_stats.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/subread.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/filter2htseq.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/combined_stats.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/shinyreports.header"

//MAIN PIPELINE TASK
run {
    "%.fastq.gz" * [
        FastQC ,
        Cutadapt + FastQQualityFilter + FilterDuplicates + TrimUMIs + [
            FastQC,
            FastQScreen,
            Bowtie_se + BAMindexer + [
                SubreadCount + Filter2HTSeq,
                bam2bw,
                subread2rnatypes
            ]
        ]
    ] +
    [
        CutadaptStats, FastQQualityFilterStats, DedupStats, MappingStats, CombinedStats
    ] +
    collectToolVersions + collectBpipeLogs + shinyReports
}
