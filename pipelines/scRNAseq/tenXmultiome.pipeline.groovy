PIPELINE="tenXmultiome"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

load PIPELINE_ROOT + "/modules/scRNAseq/cellrangerarc_count.header"
load PIPELINE_ROOT + "/modules/scRNAseq/cellrangerarc_aggr.header"
load PIPELINE_ROOT + "/modules/scRNAseq/CRmotifCounts.header"
load PIPELINE_ROOT + "/modules/scRNAseq/CTannoSeurat.header"
load PIPELINE_ROOT + "/modules/scRNAseq/CTannoMarker.header"
load PIPELINE_ROOT + "/modules/scRNAseq/demux_hto.header"
load PIPELINE_ROOT + "/modules/scRNAseq/demux_gt.header"
load PIPELINE_ROOT + "/modules/scRNAseq/assignSouporcellCluster.header"
load PIPELINE_ROOT + "/modules/scRNAseq/DNAaccess.header"
load PIPELINE_ROOT + "/modules/scRNAseq/diffPeaks.header"
load PIPELINE_ROOT + "/modules/scRNAseq/diffExprSeurat.header"
load PIPELINE_ROOT + "/modules/scRNAseq/grn.header"
load PIPELINE_ROOT + "/modules/scRNAseq/motifEnrich.header"
load PIPELINE_ROOT + "/modules/scRNAseq/motifActivity.header"
load PIPELINE_ROOT + "/modules/scRNAseq/motifFootprinting.header"
load PIPELINE_ROOT + "/modules/scRNAseq/peaks2genes.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_readAggrData.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_readIndivSamplesAndMerge.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_filter.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_qc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/SCTransform.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_integrateRNA.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_integrateATAC.header"
load PIPELINE_ROOT + "/modules/scRNAseq/wnn.header"
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
    "%.fastq.gz" * [ FastQC + FastqScreen +
      (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) ] + 
      "%_gex_S*_L*_R*_001.fastq.gz" * [
         cellrangerarc_count + [
            (RUN_DEMUX ? (RUN_DEMUX == "demux_GT" ? demux_gt : (RUN_DEMUX == "demux_HTO" ? demux_hto : dontrun.using(module:"demux") )) : dontrun.using(module:"demux")),
            bamCoverage,
            inferexperiment,
            qualimap,
            subread2rnatypes,
            geneBodyCov2
         ]
    ] + 
    (RUN_DEMUX == "demux_GT" ? assignSouporcellCluster : dontrun.using(module:"assignSouporcellCluster")) +
    (ESSENTIAL_USE_AGGR_DATA ? cellrangerarc_aggr + sc_readAggrData : sc_readIndivSamplesAndMerge ) +
    sc_qc + sc_filter +
    (ESSENTIAL_USE_AGGR_DATA ? CRmotifCounts : dontrun.using(module:"CRmotifCounts")) + 
    SCTransform + DNAaccess + 
    (RUN_BATCHCORRECT ? sc_integrateRNA + sc_integrateATAC : dontrun.using(module:"No_Batch_Correction")) +
    wnn + peaks2genes + 
    (ESSENTIAL_CELLTYPE_ANNO.contains("Seurat")? CTannoSeurat : dontrun.using(module:"CTannoSeurat")) + 
    (ESSENTIAL_CELLTYPE_ANNO.contains("Marker")? CTannoMarker : dontrun.using(module:"CTannoMarker")) + 
    [diffPeaks, diffExprSeurat] + motifActivity + motifEnrich + motifFootprinting + grn +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + MultiQC + 
    shinyReports
}
