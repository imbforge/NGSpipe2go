PIPELINE="RNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/RNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/RNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bam2bw.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/fastqc.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/filterchromosomes.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/insertsize.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/markdups2.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/trackhub.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/star.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/deseq2.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/subread.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/deseq2_mm.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/GO_Enrichment.module.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.module.groovy"
load PIPELINE_ROOT + "/modules/NGS/multiqc.module.groovy"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.groovy"
load PIPELINE_ROOT + "/modules/RNAseq/shinyreports.module.groovy"

//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
    "%.fastq.gz" * [ FastQC ] +
    (RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
        STAR + BAMindexer + [
            subread_count + filter2htseq,
            bamCoverage,
            inferexperiment,
            subread2rnatypes,
            MarkDups2 + BAMindexer + [
                dupRadar,
                geneBodyCov2
            ],
            (RUN_IN_PAIRED_END_MODE ? InsertSize : dontrun.using(module: "InsertSize"))
        ]
    ] +
    [ DE_DESeq2_MM , DE_DESeq2 + GO_Enrichment ] +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module: "trackhub")) +
    MultiQC + collectBpipeLogs + shinyReports
}

