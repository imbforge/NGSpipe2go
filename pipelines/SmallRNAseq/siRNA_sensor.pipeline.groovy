PIPELINE="siRNA sensor"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/SmallRNAseq/siRNA_sensor.essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/SmallRNAseq/tools.groovy"

MODULE_FOLDER="NGSpipe2go/modules/" // adjust to your projects needs
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

load MODULE_FOLDER + "SmallRNAseq/filter_non_structural.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/filter_non_structural.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_classes.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_classes.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_Sensor.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/filter_smRNA_Sensor.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/htseqcount.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/htseqcount.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/bamindexer.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/maping_stats.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/maping_stats.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/bam2bw.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/bam2bw.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/sensor_coverage.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/sensor_coverage.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/plot_sensor_coverage.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/plot_sensor_coverage.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/nucleotide_signature.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/nucleotide_signature.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/collect_plots.module.groovy"
load MODULE_FOLDER + "SmallRNAseq/collect_plots.vars.groovy"

load MODULE_FOLDER + "SmallRNAseq/collect_bams.module.groovy"

//MAIN PIPELINE TASK
run {
   "%.fastq.gz" *
      [ FastQC, Cutadapt + FastQQualityFilter + FilterDuplicates + TrimUMIs ] +
   "%.cutadapt.highQ.deduped.trimmed.fastq.gz" *
      [ FastQC, FastQScreen, RepEnrich, Bowtie_se + [ BAMindexer ] ] +
   "%.bam" *
         [ FilterNonStructuralReads ] +
   "%.non_structural.bam" *
         [ FilterRNAClasses, HTseqCount, FilterSensorClasses ] +
      "%.non_structural.bam" * [ Bam2bw ] +
    "%.sensor.22G.bam" * [ SensorCoverage ] +   
      [  MappingStatsPlot, PlotReadLengths ]           
}
