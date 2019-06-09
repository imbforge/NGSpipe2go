MODULE_FOLDER="NGSpipe2go/modules/" // adjust to your projects needs

load MODULE_FOLDER + "SmallRNAseq/essential.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/tool.versions.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/tool.locations.groovy"

load MODULE_FOLDER + "SmallRNAseq/fastqc.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/fastqc.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/cutadapt.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/cutadapt.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/fastq_quality_filter.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/fastq_quality_filter.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/dedup.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/dedup.module.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/dedup_stats.vars.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/dedup_stats.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/trim_umis.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/trim_umis.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/bowtie1.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/bowtie1.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/select_uniq_mappers.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/select_uniq_mappers.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/bamindexer.module.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/mapping_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/mapping_stats.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/read_count.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/read_count.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/read_count_summary.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/read_count_summary.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/count_mapped_reads.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/count_mapped_reads.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/aggregate_mapped_counts.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/aggregate_mapped_counts.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/split_read_strands.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/split_read_strands.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/bam2bw.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/bam2bw.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/nucleotide_signature.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/nucleotide_signature.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/ping_pong_signal.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/ping_pong_signal.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/ping_pong_pro.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/ping_pong_pro.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/collect_plots.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/collect_plots.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/repenrich.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/repenrich.vars.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/fastqscreen.vars.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/fastqscreen.module.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/subread.vars.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/subread.module.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/filter2htseq.vars.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/filter2htseq.module.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/cutadapt_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/cutadapt_stats.vars.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/fastq_quality_filter_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/fastq_quality_filter_stats.vars.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/combined_stats.vars.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/combined_stats.module.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/subread2rnatypes.vars.groovy"
load MODULE_FOLDER + "RNAseq/subread2rnatypes.module.groovy"

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

load MODULE_FOLDER + "smallRNAseq_BCF/shinyreports.vars.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/shinyreports.module.groovy"

//MAIN PIPELINE TASK
run {
    "%.fastq.gz" * 
    [ FastQC , Cutadapt + FastQQualityFilter + FilterDuplicates + TrimUMIs ] +
    "%.trimmed.fastq.gz" * 
    [ FastQC, FastQScreen, RepEnrich, Bowtie_se + [ BAMindexer, SelectUniqMappers + [ NucleotideSignature, PingPongSignal, PingPongPro] ] ] +
    "%.bam" *
    [ SubreadCount + Filter2HTSeq, CountReads, CountMappedReads, Bam2bw, subread2rnatypes, SplitReadStrands +
        "%sense.bam" * [Bam2bw] ] +
    [ CutadaptStats, FastQQualityFilterStats, DedupStats, MappingStats, CombinedStats, AggregateMappedCounts, CountReadsSummary ] +
    [ collectBpipeLogs + shinyReports ]
}
