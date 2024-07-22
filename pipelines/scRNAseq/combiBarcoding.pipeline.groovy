PIPELINE="ParseBio"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

load PIPELINE_ROOT + "/modules/scRNAseq/splitpipe_all.header"
load PIPELINE_ROOT + "/modules/scRNAseq/splitpipe_comb.header"
load PIPELINE_ROOT + "/modules/scRNAseq/demux_hto.header"
load PIPELINE_ROOT + "/modules/scRNAseq/demux_gt.header"
load PIPELINE_ROOT + "/modules/scRNAseq/assignSouporcellCluster.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/markdups2.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.header"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.header"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.header"
load PIPELINE_ROOT + "/modules/RNAseq/star.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread.header"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"


dontrun = { println "didn't run $module" }

Bpipe.run { 
    "%.R*.fastq.gz" * [ FastQC + FastqScreen +
         (PIPELINE == "ParseBio" ? splitpipe_all : (PIPELINE == "ScaleBio" ? dontrun.using(module:"ScaleBio not implemented yet!") : dontrun.using(module:"Pipeline type not found!") )) + [
            (RUN_DEMUX ? (RUN_DEMUX == "demux_GT" ? demux_gt : (RUN_DEMUX == "demux_HTO" ? demux_hto : dontrun.using(module:"demux") )) : dontrun.using(module:"demux")),
            bamCoverage,
            inferexperiment,
            qualimap,
            subread2rnatypes,
            geneBodyCov2
         ]
    ] + 
    (PIPELINE == "ParseBio" ? splitpipe_comb : (PIPELINE == "ScaleBio" ? dontrun.using(module:"ScaleBio not implemented yet!") : dontrun.using(module:"Pipeline type not found!") )) +
    (RUN_DEMUX == "demux_GT" ? assignSouporcellCluster : dontrun.using(module:"assignSouporcellCluster")) +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + MultiQC + 
    shinyReports
}



