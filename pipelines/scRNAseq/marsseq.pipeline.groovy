PIPELINE="scRNAseq_marsseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bamcoverage.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/markdups2.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/star.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/addumibarcodetofastq.module.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/cutadapt.module.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/subread.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.module.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/umicount.module.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/umidedup.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.module.groovy"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.module.groovy"

//
// Typical workflow for MARS-Seq data:
//
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%.fastq.gz" * [ FastQC ] +
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
        AddUMIBarcodeToFastq + Cutadapt + [
            FastQC,
            STAR + BAMindexer + [
                subread_count + BAMindexer + umicount,
                bamCoverage,
                inferexperiment,
                subread2rnatypes,
                qualimap,
                geneBodyCov2
            ]
        ]
    ] +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + collectBpipeLogs + shinyReports
}

