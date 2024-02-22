PIPELINE="RNAseqVariantCalling"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/" // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/RNAseqVariantCalling/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/RNAseqVariantCalling/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/star1pass.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/merge_SJ_tab.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/star2pass.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/add_read_group.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/mark_dups.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/splitNcigar.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/base_recalibration.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/create_star_index_sjdb.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/variantCall_HC.header"
load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/variant_filtration.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"

//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%R*.fastq.gz" * [ STAR_pe ] +
    "*.SJ.out.tab" * [ FilterAndMergeSJtab + GenerateStarIndexFromSJ ] +
    "%R*.fastq.gz" * [
        STAR_pe_2nd + AddRG + MarkDups + SplitNCigarReads + BaseRecalibration + VariantCallHC + VariantFiltration
    ] + collectToolVersions
}
