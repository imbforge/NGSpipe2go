MODULE_FOLDER="NGSpipe2go/modules/" // adjust to your projects needs

load MODULE_FOLDER + "SmallRNAseq/essential.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/fastqc.vars.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"

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

load MODULE_FOLDER + "SmallRNAseq/mirDeep2_mapper.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/mirDeep2_mapper.module.groovy"

load MODULE_FOLDER + "SmallRNAseq/mirDeep2.vars.groovy"
load MODULE_FOLDER + "SmallRNAseq/mirDeep2.module.groovy"


//MAIN PIPELINE TASK
run {
	"%.fastq.gz" * [ FastQC , Cutadapt + FastQQualityFilter + FilterDuplicates + TrimUMIs ] + "%.deduped_barcoded.trimmed.fastq" * [ FastQC + MirDeep2Mapper + MirDeep2] + [ DedupStats ]
}
