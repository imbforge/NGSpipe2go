MODULE_FOLDER="./NGSpipe2go/modules/"  // adjust to your projects needs

load MODULE_FOLDER + "DNAseq/essential.vars.groovy"
load MODULE_FOLDER + "DNAseq/tool.locations.groovy"
load MODULE_FOLDER + "DNAseq/tool.versions.groovy"

load MODULE_FOLDER + "DNAseq/bwa.module.groovy"
load MODULE_FOLDER + "DNAseq/realignment.module.groovy"
load MODULE_FOLDER + "DNAseq/recalibration.module.groovy"
load MODULE_FOLDER + "DNAseq/shinyreports.module.groovy"
load MODULE_FOLDER + "DNAseq/variantcallHC.module.groovy"
load MODULE_FOLDER + "DNAseq/variantcallUG.module.groovy"
load MODULE_FOLDER + "DNAseq/varianteval.module.groovy"
load MODULE_FOLDER + "DNAseq/variantfuseHC.module.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"
load MODULE_FOLDER + "NGS/rmdups.module.groovy"
load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

run {
    "%.fastq.gz" * [ FastQC ] + "%_R*.fastq.gz" * [ BWA_pe ] + "%.bam" * [ RmDups + BAMindexer + IndelRealignment + BaseRecalibration + [ VariantCallHC, VariantCallUG ] ] + "%.vcf.gz" * [ VariantEval ] + collectBpipeLogs + shinyReports
}
