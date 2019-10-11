PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamcoverage.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/markdups2.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/star.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/subread.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/cutadapt.module.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.module.groovy"


//
// Typical workflow for SmartSeq data:
//
dontrun = { println "didn't run $module" }

run { 
    "%.fastq.gz" * [ FastQC ] + 
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [ Cutadapt + FastQC + STAR + BAMindexer + 
        [ subread_count +  bamCoverage , filter2htseq, inferexperiment , subread2rnatypes , qualimap, MarkDups2 + BAMindexer +
        [ dupRadar , geneBodyCov2 ]]] +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectBpipeLogs + shinyReports
}


