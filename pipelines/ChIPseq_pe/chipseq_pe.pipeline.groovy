MODULE_FOLDER="./NGSpipe2go/modules/"

load MODULE_FOLDER + "ChIPseq/essential.vars.groovy"
load MODULE_FOLDER + "ChIPseq/tool.locations.groovy"
load MODULE_FOLDER + "ChIPseq/tool.versions.groovy"

load MODULE_FOLDER + "NGS/fastqc.vars.groovy"
load MODULE_FOLDER + "NGS/fastqc.module.groovy"

load MODULE_FOLDER + "ChIPseq/bowtie2pe.vars.groovy"
load MODULE_FOLDER + "ChIPseq/bowtie2pe.module.groovy"

load MODULE_FOLDER + "ChIPseq/filbowtie2unique.vars.groovy"
load MODULE_FOLDER + "ChIPseq/filbowtie2unique.module.groovy"

load MODULE_FOLDER + "NGS/rmdups.vars.groovy"
load MODULE_FOLDER + "NGS/rmdups.module.groovy"

load MODULE_FOLDER + "NGS/markdups.vars.groovy"
load MODULE_FOLDER + "NGS/markdups.module.groovy"

load MODULE_FOLDER + "NGS/bamindexer.vars.groovy"
load MODULE_FOLDER + "NGS/bamindexer.module.groovy"

load MODULE_FOLDER + "NGS/bamqc.vars.groovy"
load MODULE_FOLDER + "NGS/bamqc.module.groovy"

load MODULE_FOLDER + "NGS/insertsize.vars.groovy"
load MODULE_FOLDER + "NGS/insertsize.module.groovy"

load MODULE_FOLDER + "NGS/bamcoverage.vars.groovy"
load MODULE_FOLDER + "NGS/bamcoverage.module.groovy"

load MODULE_FOLDER + "NGS/trackhub_config.vars.groovy"
load MODULE_FOLDER + "NGS/trackhub_config.module.groovy"

load MODULE_FOLDER + "NGS/trackhub.vars.groovy"
load MODULE_FOLDER + "NGS/trackhub.module.groovy"

load MODULE_FOLDER + "ChIPseq/ipstrength.vars.groovy"
load MODULE_FOLDER + "ChIPseq/ipstrength.module.groovy"

load MODULE_FOLDER + "ChIPseq/pbc.module.groovy"

load MODULE_FOLDER + "ChIPseq/macs2.vars.groovy"
load MODULE_FOLDER + "ChIPseq/macs2.module.groovy"

load MODULE_FOLDER + "ChIPseq/blacklist_filter.vars.groovy"
load MODULE_FOLDER + "ChIPseq/blacklist_filter.module.groovy"

load MODULE_FOLDER + "ChIPseq/peak_annotation.vars.groovy"
load MODULE_FOLDER + "ChIPseq/peak_annotation.module.groovy"  

load MODULE_FOLDER + "ChIPseq/diffbind.vars.groovy"
load MODULE_FOLDER + "ChIPseq/diffbind.module.groovy"

load MODULE_FOLDER + "ChIPseq/GREAT.vars.groovy"
load MODULE_FOLDER + "ChIPseq/GREAT.module.groovy"             

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"

load MODULE_FOLDER + "ChIPseq/shinyreports_pe.vars.groovy"
load MODULE_FOLDER + "ChIPseq/shinyreports_pe.module.groovy"


// standard PE workflow using deduplicated BAM files of MapQ filtered (~ unique) reads
// assumes that duplicated library fragments are unwanted PCR artifacts
// discards multimapping reads as habitually done in most ChIP-seq studies
//
run {
    "%.fastq.gz" * [ FastQC ] + "%.R*.fastq.gz" * 
    [ bowtie2_pe + BAMindexer + BamQC + filbowtie2unique + BAMindexer + RmDups + BAMindexer +
    [ bamCoverage, InsertSize, ipstrength, macs2 ] ]  + 
    // diffbind +
    // trackhub_config + trackhub +
    peak_annotation + collectBpipeLogs + shinyReports
}

// alternative PE workflow using the unfiltered BAM files
// may be preferable when studying some types of repetetive regions
// make sure to have ESSENTIAL_DUP="auto" for MACS2 peak calling
//
// run {
//    "%.fastq.gz" * [ FastQC ] + "%.R*.fastq.gz" * 
//    [ bowtie2_pe + BAMindexer + BamQC + 
//    [ MarkDups + BAMindexer, bamCoverage, InsertSize, ipstrength, macs2 ] ] + 
//    diffbind +
//    // trackhub_config + trackhub +
//    peak_annotation + collectBpipeLogs + shinyReports
//}

