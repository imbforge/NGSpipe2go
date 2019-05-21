MODULE_FOLDER="./NGSpipe2go/modules/"    // may need adjustment for some projects

load MODULE_FOLDER + "RNAseq/essential.vars.groovy"
load MODULE_FOLDER + "RNAseq/tool.locations.groovy"
load MODULE_FOLDER + "RNAseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/fastqc.vars.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"

load MODULE_FOLDER + "RNAseq/star.vars.groovy"
load MODULE_FOLDER + "RNAseq/star.module.groovy"

load MODULE_FOLDER + "NGS/filterchromosomes.vars.groovy"
load MODULE_FOLDER + "NGS/filterchromosomes.module.groovy"

load MODULE_FOLDER + "NGS/bamindexer.vars.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"

load MODULE_FOLDER + "NGS/insertsize.module.groovy"
load MODULE_FOLDER + "NGS/insertsize.vars.groovy"

load MODULE_FOLDER + "RNAseq/subread.vars.groovy"
load MODULE_FOLDER + "RNAseq/subread.module.groovy"

load MODULE_FOLDER + "RNAseq/filter2htseq.vars.groovy"
load MODULE_FOLDER + "RNAseq/filter2htseq.module.groovy"

load MODULE_FOLDER + "NGS/bam2bw.vars.groovy"
load MODULE_FOLDER + "NGS/bam2bw.module.groovy"

load MODULE_FOLDER + "NGS/bamcoverage.module.groovy"
load MODULE_FOLDER + "NGS/bamcoverage.vars.groovy"

load MODULE_FOLDER + "NGS/trackhub_config.vars.groovy"
load MODULE_FOLDER + "NGS/trackhub_config.module.groovy"

load MODULE_FOLDER + "NGS/trackhub.vars.groovy"
load MODULE_FOLDER + "NGS/trackhub.module.groovy"

load MODULE_FOLDER + "NGS/markdups2.vars.groovy"
load MODULE_FOLDER + "NGS/markdups2.module.groovy"

load MODULE_FOLDER + "RNAseq/dupradar.vars.groovy"
load MODULE_FOLDER + "RNAseq/dupradar.module.groovy"

load MODULE_FOLDER + "RNAseq/genebodycov2.vars.groovy"
load MODULE_FOLDER + "RNAseq/genebodycov2.module.groovy"

load MODULE_FOLDER + "RNAseq/inferexperiment.vars.groovy"
load MODULE_FOLDER + "RNAseq/inferexperiment.module.groovy"

load MODULE_FOLDER + "RNAseq/subread2rnatypes.vars.groovy"
load MODULE_FOLDER + "RNAseq/subread2rnatypes.module.groovy"

load MODULE_FOLDER + "RNAseq/deseq2.vars.groovy"
load MODULE_FOLDER + "RNAseq/deseq2.module.groovy"

load MODULE_FOLDER + "RNAseq/deseq2_mm.vars.groovy"
load MODULE_FOLDER + "RNAseq/deseq2_mm.module.groovy"

load MODULE_FOLDER + "RNAseq/GO_Enrichment.vars.groovy"
load MODULE_FOLDER + "RNAseq/GO_Enrichment.module.groovy"

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

load MODULE_FOLDER + "RNAseq/shinyreports.vars.groovy"
load MODULE_FOLDER + "RNAseq/shinyreports.module.groovy"

//
// Typical workflow for SR data (default):
//
run { "%.fastq.gz" * 
[ FastQC , STAR + BAMindexer + 
[ subread_count + filter2htseq , bam2bw , inferexperiment , subread2rnatypes , MarkDups2 + BAMindexer + 
[ dupRadar , geneBodyCov2 ] ] ] + 
[ DE_DESeq2_MM , DE_DESeq2 + GO_Enrichment] + 
//trackhub_config + trackhub +
collectBpipeLogs + shinyReports }

//
// Typical workflow for PE data (optional):
//
// run { "%.fastq.gz" * [ FastQC ] + "%.R*.fastq.gz" * [ STAR + BAMindexer + [ subread_count + filter2htseq , bamCoverage , InsertSize , inferexperiment , subread2rnatypes , MarkDups2 + BAMindexer + [ dupRadar , geneBodyCov2 ] ] ] + [ DE_DESeq2_MM , DE_DESeq2 + GO_Enrichment] + trackhub_config + trackhub + collectBpipeLogs + shinyReports }
