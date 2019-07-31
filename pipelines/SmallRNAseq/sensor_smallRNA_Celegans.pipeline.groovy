MODULE_FOLDER="NGSpipe2go/modules/" // adjust to your projects needs

load MODULE_FOLDER + "SmallRNAseq/essential.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/tool.versions.groovy"

load MODULE_FOLDER + "RNAseq/subread2rnatypes.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/bam2bw.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/bamindexer.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/bowtie1.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/collect_bams.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/count_non_struct_reads.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/count_read_lengths.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/cutadapt.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/dedup.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/fastq_quality_filter.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/fastqc.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_classes.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/htseqcount.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/maping_stats.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/nucleotide_signature.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/ping_pong_pro.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/ping_pong_signal.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/plot_read_lengths.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/plot_sensor_coverage.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/read_count.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/read_count_summary.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/repenrich.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/select_uniq_mappers.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/sensor_coverage.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/split_read_strands.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/subread.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/trim_umis.module.groovy"
load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/combined_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/cutadapt_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/dedup_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/fastq_quality_filter_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/fastqscreen.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/filter2htseq.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/mapping_stats.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/shinyreports.module.groovy"
load MODULE_FOLDER + "smallRNAseq_BCF/subread.module.groovy"

//MAIN PIPELINE TASK
run {
   "%.bam" * [ HTseqCount ] + CountNonStrutReads +
   "%.bam" * [ FilterRNAClasses ] +
   "%.bam" * [ Bam2bw, SensorCoverage ] + [ PlotSensorCoverage ]
}

/*//MAIN PIPELINE TASK
 run { "%.bam" * [ FilterRNAClasses, HTseqCount ] +
    [ CountNonStrutReads , "%.22G.bam" * [ Bam2bw, SensorCoverage + PlotSensorCoverage ] ] }*/
