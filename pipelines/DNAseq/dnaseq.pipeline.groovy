PIPELINE="DNAseq"
PIPELINE_VERSION="2.0"
PIPELINE_ROOT="./NGSpipe2go/"  // adjust to your projects needs

load PIPELINE_ROOT + "/pipelines/DNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/DNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/DNAseq/bwa.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/DNAseq/realignment.header"
load PIPELINE_ROOT + "/modules/DNAseq/recalibration.header"
load PIPELINE_ROOT + "/modules/DNAseq/gatherBQSRReports.header"
load PIPELINE_ROOT + "/modules/DNAseq/variantcallHC.header"
load PIPELINE_ROOT + "/modules/DNAseq/genomicsDBImport.header"
load PIPELINE_ROOT + "/modules/DNAseq/genotypeGVCFs.header"
load PIPELINE_ROOT + "/modules/DNAseq/collectVariantCallingMetrics.header"
load PIPELINE_ROOT + "/modules/DNAseq/validateVariants.header"
load PIPELINE_ROOT + "/modules/DNAseq/varianteval.header"
load PIPELINE_ROOT + "/modules/DNAseq/variant_score_recalibration.header"
load PIPELINE_ROOT + "/modules/DNAseq/variantFiltration.header"
load PIPELINE_ROOT + "/modules/DNAseq/snpEff.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/NGS/rmdups.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/DNAseq/shinyreports.header"

// Main pipeline task
dontrun = { println "didn't run $module" }
collect_bams = { forward inputs.bam }

Bpipe.run {
    "%.fastq.gz" * [ FastQC, (RUN_FASTQSCREEN ? FastqScreen : dontrun.using(module: "FastqScreen")) ] +
    (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) +
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" * [ BWA_pe + BAMindexer ] : "%.fastq.gz" * [ BWA_se + BAMindexer ] ) + 
    collect_bams + "%.bam" * [ 
        (RUN_RMDUPS ? RmDups + BAMindexer : dontrun.using(module:"RmDups")) + [
            bamCoverage,
            BaseRecalibration + VariantCallHC 
            ]
        ] + 
    GenomicsDBImport + GenotypeGVCFs + 
    ( ESSENTIAL_FILTERVARIANTS=="VQSR" ? VariantScoreRecalibration : (ESSENTIAL_FILTERVARIANTS=="hard-filter" ? VariantFiltration : dontrun.using(module: "variant filtering"))) +
    [ VariantEval, CollectVariantCallingMetrics, GatherBQSRReports ] + 
    (RUN_SNPEFF ? snpEff : dontrun.using(module:"snpEff")) +
    collectToolVersions + MultiQC + shinyReports
}

