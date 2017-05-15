MODULE_FOLDER="NGSpipe2go/modules/" // adjust to your projects needs

load MODULE_FOLDER + "SmallRNAseq/essential.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/tool.locations.groovy"

load MODULE_FOLDER + "SmallRNAseq/fastqc.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/fastqc.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/cutadapt.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/cutadapt.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/fastq_quality_filter.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/fastq_quality_filter.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/dedup.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/dedup.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/dedup_stats.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/dedup_stats.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/trim_umis.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/trim_umis.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/count_read_lengths.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/count_read_lengths.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/plot_read_lengths.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/plot_read_lengths.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/bowtie1.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/bowtie1.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/select_uniq_mappers.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/select_uniq_mappers.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_classes.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_classes.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/subread.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/subread.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/htseqcount.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/htseqcount.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/bamindexer.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/maping_stats.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/maping_stats.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/read_count.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/read_count.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/count_non_struct_reads.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/count_non_struct_reads.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/read_count_summary.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/read_count_summary.vars.groovy"

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

load MODULE_FOLDER + "SmallRNAseq/repenrich.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/repenrich.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/collect_plots.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/collect_plots.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/collect_bams.module.groovy"

//MAIN PIPELINE TASK
run {
	"%.fastq.gz" *
      [ FastQC, Cutadapt + FastQQualityFilter + FilterDuplicates + TrimUMIs ] +
   "%.deduped_barcoded.trimmed.fastq.gz" *
      [ FastQC, CountReadLengths, RepEnrich, Bowtie_se + [ BAMindexer ] ] +
   "%.bam" * 
   		[ FilterRNAClasses +
   		"%.bam" * 
   			[ HTseqCount ], SplitReadStrands + 
   		"%sense.bam" * [Bam2bw] 
   		] +
   SubReadCount + CountNonStrutReads + PlotReadLengths
}

