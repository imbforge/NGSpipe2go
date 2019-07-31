MODULE_FOLDER="./NGSpipe2go/modules/"    // may need adjustment for some projects

load MODULE_FOLDER + "scRNAseq/essential.vars.groovy"
load MODULE_FOLDER + "scRNAseq/tool.locations.groovy"
load MODULE_FOLDER + "scRNAseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/bamcoverage.module.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"
load MODULE_FOLDER + "NGS/markdups2.module.groovy"
load MODULE_FOLDER + "RNAseq/dupradar.module.groovy"
load MODULE_FOLDER + "RNAseq/genebodycov2.module.groovy"
load MODULE_FOLDER + "RNAseq/inferexperiment.module.groovy"
load MODULE_FOLDER + "RNAseq/qualimap.module.groovy"
load MODULE_FOLDER + "RNAseq/star.module.groovy"
load MODULE_FOLDER + "RNAseq/subread2rnatypes.module.groovy"
load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"
load MODULE_FOLDER + "scRNAseq/addumibarcodetofastq.module.groovy"
load MODULE_FOLDER + "scRNAseq/cutadapt.module.groovy"
load MODULE_FOLDER + "scRNAseq/shinyreports.module.groovy"
load MODULE_FOLDER + "scRNAseq/subread.module.groovy"
load MODULE_FOLDER + "scRNAseq/umicount.module.groovy"
load MODULE_FOLDER + "scRNAseq/umidedup.module.groovy"

//
// Typical workflow for MARS-Seq data:
//
run {
    "%.fastq.gz" * [ FastQC ] +
    "%.R*.fastq.gz" * [ AddUMIBarcodeToFastq + Cutadapt + FastQC + STAR + BAMindexer +
                      [ subread_count + BAMindexer + umicount , bamCoverage , inferexperiment , subread2rnatypes , qualimap, geneBodyCov2 ]] +
    //trackhub_config + trackhub +
    collectBpipeLogs + shinyReports
}

