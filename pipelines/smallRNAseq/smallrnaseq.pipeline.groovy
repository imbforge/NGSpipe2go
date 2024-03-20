PIPELINE="smallRNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/smallRNAseq/smallrnaseq.essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/smallRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"

load PIPELINE_ROOT + "/modules/smallRNAseq/bowtie1.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/dedup.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/deseq2.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/deseq2_mirnamature.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/fastq_quality_filter.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/filter2htseq.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/filter_smallrna_counts.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/subread.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/subread_mirnamature.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/smallRNAseq/trim_umis.header"

//load PIPELINE_ROOT + "/modules/smallRNAseq/mirDeep2.header"
//load PIPELINE_ROOT + "/modules/smallRNAseq/mirDeep2_mapper.header"


//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

run {
    "%.fastq.gz" * 
    [
        FastQC.using(subdir:"raw") ,
        Cutadapt + FastQQualityFilter + (REMOVE_DUPLICATES ? FilterDuplicates : dontrun.using(module: "FilterDuplicates")) + (TRIM_UMIS ? TrimUMIs : dontrun.using(module: "TrimUMIs")) +
        [
            FastQC.using(subdir:"trimmed"),
            (RUN_FASTQSCREEN ? FastqScreen : dontrun.using(module: "FastqScreen")),
            bowtie1_sRNA + BAMindexer + 
            [
                subread_count.using(subdir:"all") + filter2htseq.using(subdir:"all") + filter_smallRNA_counts.using(subdir:ESSENTIAL_SMALLRNA),
         	(RUN_MATUREMIRNA_ANALYSIS ? subread_miRNAmature_count.using(subdir:"miRNAmature") + filter2htseq.using(subdir:"miRNAmature") 
                                          : dontrun.using(module: "subread_miRNAmature_count")),
                bamCoverage,
                subread2rnatypes
            ]
        ]
    ] + 
    [ 
      DE_DESeq2.using(subdir:"all"),
      DE_DESeq2.using(subdir:ESSENTIAL_SMALLRNA), 
      (RUN_MATUREMIRNA_ANALYSIS ? DE_DESeq2_miRNAmature.using(subdir:"miRNAmature") : dontrun.using(module: "DE_DESeq2_miRNAmature"))
    ] +
    collectToolVersions + shinyReports 
}


// not tested pipeline using miRDeep
//run {
//    "%.fastq.gz" * [
//        FastQC.using(subdir:"raw"),
//        Cutadapt + FastQQualityFilter + (REMOVE_DUPLICATES ? FilterDuplicates : dontrun.using(module: "FilterDuplicates")) + TrimUMIs +
//        [
//            FastQC.using(subdir:"trimmed"),
//            miRDeep2Mapper + miRDeep2
//        ]
//    ] +
//    collectToolVersions
//}

