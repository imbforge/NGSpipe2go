PIPELINE="tenxSTARRseq"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/scRNAseq/cellranger_count.header"
load PIPELINE_ROOT + "/modules/scRNAseq/cellranger_aggr.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/markdups2.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/RNAseq/dupradar.header"
load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.header"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.header"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.header"
load PIPELINE_ROOT + "/modules/RNAseq/star.header"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/scRNAseq/addumibarcodetofastq.header"
load PIPELINE_ROOT + "/modules/scRNAseq/subread.header"
load PIPELINE_ROOT + "/modules/scRNAseq/umicount.header"
load PIPELINE_ROOT + "/modules/scRNAseq/umidedup.header"
load PIPELINE_ROOT + "/modules/tenxSTARRseq/splitmrna.header"
load PIPELINE_ROOT + "/modules/tenxSTARRseq/mergestarrcounts.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"


dontrun = { println "didn't run $module" }

Bpipe.run { 
    "%.fastq.gz" * [ FastQC + FastqScreen +
      (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) ] + 
      "%_R*_001.fastq.gz" * [ SplitmRNA ] +
      [         // parallel branches with endogenous and STARR mRNA processing
          "%_S*_L*_R*_001.endogenous.fastq.gz" * [    // Cell Ranger branch for endogenous mRNA
              cellranger_count + [
                 bamCoverage,
                 inferexperiment,
                 qualimap,
                 subread2rnatypes,
                 dupRadar,
                 geneBodyCov2
              ]
          ] + cellranger_aggr,
          "%_S*_L*_R*_001.starr.fastq.gz" * [        // STARR mRNA branch: map with STAR & count UMIs
              AddUMIBarcodeToFastq +
              STAR + BAMindexer + [
                  subread_count + BAMindexer + umicount
              ]
          ]
      ] + MergeSTARRCounts +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + collectBpipeLogs + MultiQC + shinyReports
}

