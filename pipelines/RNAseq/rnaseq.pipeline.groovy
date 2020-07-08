PIPELINE="RNAseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/RNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/RNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/NGS/bam2bw.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/filterchromosomes.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/markdups2.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.header"
load PIPELINE_ROOT + "/modules/RNAseq/star.header"
load PIPELINE_ROOT + "/modules/RNAseq/deseq2.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.header"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.header"
load PIPELINE_ROOT + "/modules/RNAseq/deseq2_mm.header"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.header"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.header"
load PIPELINE_ROOT + "/modules/RNAseq/GO_Enrichment.header"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/RNAseq/shinyreports.header"

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
    collectToolVersions + collectBpipeLogs + MultiQC + shinyReports
}

