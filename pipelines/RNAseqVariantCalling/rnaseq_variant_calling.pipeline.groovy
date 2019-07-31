MODULE_FOLDER="./NGSpipe2go/modules/" // adjust to your projects needs

load MODULE_FOLDER + "RNAseqVariantCalling/essential.vars.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/tool.locations.groovy"

load MODULE_FOLDER + "NGS/bamindexer.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/add_read_group.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/base_recalibration.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/create_star_index_sjdb.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/mark_dups.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/merge_SJ_tab.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/splitNcigar.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/star1pass.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/star2pass.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/variantCall_HC.module.groovy"
load MODULE_FOLDER + "RNAseqVariantCalling/variant_filtration.module.groovy"

run {
    "%R*.fastq.gz" * [ STAR_pe ] +
    "*.SJ.out.tab" * [ FilterAndMergeSJtab + GenerateStarIndexFromSJ ] +
    "%R*.fastq.gz" * [ STAR_pe_2nd ] +
    "%.bam" * [ AddRG + MarkDups + SplitNCigarReads + BaseRecalibration + VariantCallHC ] +
    "%.vcf.gz" * [ VariantFiltration ]
}
