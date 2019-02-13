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

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

//
// Typical workflow for MARS-Seq data:
//
run { "%.fastq.gz" *  [ FastQC ] + "%.R*.fastq.gz" * [ AddUMIBarcodeToFastq + Cutadapt + STAR + BAMindexer +  subread_count + BAMindexer + umidedup +  BAMindexer + umicount 
] +
//trackhub_config + trackhub +
collectBpipeLogs
// will be replaced by the reports file for the single cell? + shinyReports
}

