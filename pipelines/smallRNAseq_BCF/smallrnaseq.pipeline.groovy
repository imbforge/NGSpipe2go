PIPELINE="smallRNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/smallrnaseq.essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/smallRNAseq_BCF/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/bowtie1.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/deseq2.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/deseq2_mirnamature.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/filter2htseq.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/filter_mirna_counts.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/subread.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/subread_mirnamature.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/shinyreports.header"
load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/trim_umis.header"


//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

run {
    "%.fastq.gz" * 
    [
        FastQC.using(subdir:"raw") ,
        Cutadapt + FastQQualityFilter + (REMOVE_DUPLICATES ? FilterDuplicates : dontrun.using(module: "FilterDuplicates")) + TrimUMIs +
        [
            FastQC.using(subdir:"trimmed"),
            (RUN_FASTQSCREEN ? FastqScreen : dontrun.using(module: "FastqScreen")),
            bowtie1_sRNA + BAMindexer + 
            [
                subread_count.using(subdir:"all") + filter2htseq.using(subdir:"all") + filter_miRNA_counts.using(subdir:"miRNA"),
         	subread_miRNAmature_count.using(subdir:"miRNAmature") + filter2htseq.using(subdir:"miRNAmature"),
                bamCoverage,
                subread2rnatypes
            ]
        ]
    ] + 
    [ 
      DE_DESeq2.using(subdir:"all"),
      DE_DESeq2.using(subdir:"miRNA"), 
      DE_DESeq2_miRNAmature.using(subdir:"miRNAmature")
    ] +
    collectToolVersions + shinyReports 
}
