PIPELINE="miRNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/miRNA.essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/cutadapt.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup_stats.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mirDeep2.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mirDeep2_mapper.module.groovy"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/trim_umis.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.module.groovy"

//MAIN PIPELINE TASK
run {
    "%.fastq.gz" * [
        FastQC,
        Cutadapt + [
            FastQQualityFilter + FilterDuplicates + TrimUMIs + [
                FastQC,
                miRDeep2Mapper + miRDeep2
            ]
        ]
    ] +
    collectToolVersions
}
