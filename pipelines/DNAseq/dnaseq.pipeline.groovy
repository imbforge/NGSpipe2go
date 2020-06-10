PIPELINE="DNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"  // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/DNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/DNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/DNAseq/bwa.header"
load PIPELINE_ROOT + "/modules/DNAseq/realignment.header"
load PIPELINE_ROOT + "/modules/DNAseq/recalibration.header"
load PIPELINE_ROOT + "/modules/DNAseq/variantcallHC.header"
load PIPELINE_ROOT + "/modules/DNAseq/variantcallUG.header"
load PIPELINE_ROOT + "/modules/DNAseq/varianteval.header"
load PIPELINE_ROOT + "/modules/DNAseq/variantfuseHC.header"
load PIPELINE_ROOT + "/modules/DNAseq/variant_score_recalibration.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/rmdups.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/DNAseq/shinyreports.header"

// Main pipeline task
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%.fastq.gz" * [ FastQC ] +
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" * [ BWA_pe ] : "%.fastq.gz" * [ BWA_se ] ) +
    "%.bam" * [
        RmDups + BAMindexer + IndelRealignment + BaseRecalibration + [
            VariantCallHC + [
                VariantEval,
                VariantScoreRecalibration + VariantEval
            ],
            VariantCallUG + [
                VariantEval,
                VariantScoreRecalibration + VariantEval
            ]
        ]
    ] +
    collectToolVersions + collectBpipeLogs + shinyReports
}
