PIPELINE="smallRNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"

load PIPELINE_ROOT + "/modules/NGS/bam2bw.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/bowtie1.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/combined_stats.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/cutadapt.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/cutadapt_stats.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup_stats.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter_stats.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastqscreen.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/filter2htseq.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mapping_stats.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/subread.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/trim_umis.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/shinyreports.module.groovy"

//MAIN PIPELINE TASK
run {
    "%.fastq.gz" * [ FastQC , Cutadapt ] +
    "%.cutadapt.fastq.gz" * [ FastQQualityFilter + FilterDuplicates + TrimUMIs + [ FastQC , Bowtie_se + BAMindexer , FastQScreen ] ] +
    "%.bam" * [ SubreadCount + Filter2HTSeq, bam2bw, rnatypes ] + [ CutadaptStats, FastQQualityFilterStats, DedupStats, MappingStats, CombinedStats ] +
    collectBpipeLogs + shinyReports
}
