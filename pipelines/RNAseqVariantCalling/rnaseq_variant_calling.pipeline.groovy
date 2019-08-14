PIPELINE_ROOT="./NGSpipe2go/" // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/RNAseqVariantCalling/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/RNAseqVariantCalling/tools.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/add_read_group.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/base_recalibration.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/create_star_index_sjdb.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/mark_dups.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/merge_SJ_tab.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/splitNcigar.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/star1pass.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/star2pass.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/variantCall_HC.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/variant_filtration.module.groovy"

//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%R*.fastq.gz" * [ STAR_pe ] +
    "*.SJ.out.tab" * [ FilterAndMergeSJtab + GenerateStarIndexFromSJ ] +
    "%R*.fastq.gz" * [ STAR_pe_2nd ] +
    "%.bam" * [ AddRG + MarkDups + SplitNCigarReads + BaseRecalibration + VariantCallHC ] +
    "%.vcf.gz" * [ VariantFiltration ]
}
