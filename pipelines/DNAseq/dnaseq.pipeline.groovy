PIPELINE_ROOT="./NGSpipe2go/"  // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/DNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/DNAseq/tools.groovy"

load PIPELINE_ROOT + "/modules/DNAseq/bwa.module.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/realignment.module.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/recalibration.module.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/variantcallHC.module.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/variantcallUG.module.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/varianteval.module.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/variantfuseHC.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/rmdups.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/DNAseq/shinyreports.module.groovy"

run {
    "%.fastq.gz" * [ FastQC ] + "%_R*.fastq.gz" * [ BWA_pe ] + "%.bam" * [ RmDups + BAMindexer + IndelRealignment + BaseRecalibration + [ VariantCallHC, VariantCallUG ] ] + "%.vcf.gz" * [ VariantEval ] + collectBpipeLogs + shinyReports
}
