MODULE_FOLDER="./NGSpipe2go/modules/"    // may need adjustment for some projects

load MODULE_FOLDER + "MARSseq/essential.vars.groovy"
load MODULE_FOLDER + "MARSseq/tool.locations.groovy"
load MODULE_FOLDER + "MARSseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/fastqc.vars.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"

load MODULE_FOLDER + "RNAseq/star.vars.groovy"
load MODULE_FOLDER + "RNAseq/star.module.groovy"

load MODULE_FOLDER + "MARSseq/addumibarcodetofastq.vars.groovy"
load MODULE_FOLDER + "MARSseq/addumibarcodetofastq.module.groovy"

load MODULE_FOLDER + "MARSseq/cutadapt.vars.groovy"
load MODULE_FOLDER + "MARSseq/cutadapt.module.groovy"
load MODULE_FOLDER + "MARSseq/umidedup.vars.groovy"
load MODULE_FOLDER + "MARSseq/umidedup.module.groovy"
load MODULE_FOLDER + "MARSseq/umicount.vars.groovy"
load MODULE_FOLDER + "MARSseq/umicount.module.groovy"

load MODULE_FOLDER + "NGS/bamindexer.vars.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"

load MODULE_FOLDER + "MARSseq/subread.vars.groovy"
load MODULE_FOLDER + "MARSseq/subread.module.groovy"

load MODULE_FOLDER + "RNAseq/genebodycov2.vars.groovy"
load MODULE_FOLDER + "RNAseq/genebodycov2.module.groovy"

load MODULE_FOLDER + "RNAseq/inferexperiment.vars.groovy"
load MODULE_FOLDER + "RNAseq/inferexperiment.module.groovy"

load MODULE_FOLDER + "RNAseq/subread2rnatypes.vars.groovy"
load MODULE_FOLDER + "RNAseq/subread2rnatypes.module.groovy"

load MODULE_FOLDER + "NGS/bam2bw.vars.groovy"
load MODULE_FOLDER + "NGS/bam2bw.module.groovy"

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

load MODULE_FOLDER + "MARSseq/shinyreports.vars.groovy"
load MODULE_FOLDER + "MARSseq/shinyreports.module.groovy"


//
// Typical workflow for MARS-Seq data:
//
run { "%.fastq.gz" *  [ FastQC ] + "%.R*.fastq.gz" * [ AddUMIBarcodeToFastq + Cutadapt + FastQC + STAR + BAMindexer + 
[ subread_count + BAMindexer + umicount , bam2bw , inferexperiment , subread2rnatypes, geneBodyCov2 ]] +
//trackhub_config + trackhub +
collectBpipeLogs + shinyReports
}

