PIPELINE="ChIPseq"
PIPELINE_VERSION="1.2.2"
PIPELINE_ROOT="./NGSpipe2go"

load PIPELINE_ROOT + "/pipelines/ChIPseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/ChIPseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/ChIPseq/GREAT.header"
load PIPELINE_ROOT + "/modules/ChIPseq/blacklist_filter.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1.header"
load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/diffbind.header"
load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.header"
load PIPELINE_ROOT + "/modules/ChIPseq/ipstrength.header"
load PIPELINE_ROOT + "/modules/ChIPseq/macs2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/pbc.header"
load PIPELINE_ROOT + "/modules/ChIPseq/peak_annotation.header"
load PIPELINE_ROOT + "/modules/ChIPseq/phantompeak.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/extend.header"
load PIPELINE_ROOT + "/modules/NGS/bamqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/markdups.header"
load PIPELINE_ROOT + "/modules/NGS/rmdups.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub.header"
load PIPELINE_ROOT + "/modules/NGS/trackhub_config.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collectbpipes.module.2.header"
load PIPELINE_ROOT + "/modules/ChIPseq/shinyreports.header"


//MAIN PIPELINE TASK
dontrun = { println "didn't run $module" }

Bpipe.run {
	(RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [
		FastQC +  
                (RUN_CUTADAPT ? Cutadapt + FastQC : dontrun.using(module:"Cutadapt")) +
                (ESSENTIAL_USE_BOWTIE1 ? bowtie1 : bowtie2) + BAMindexer + BamQC +   

         [ // parallel branches with and without multi mappers

            [ MarkDups + BAMindexer + pbc +
		[       // QC specific to paired end (pe) or single end (se) design
	    		(RUN_IN_PAIRED_END_MODE ? [bamCoverage.using(subdir:"unfiltered"), 
                                                   InsertSize.using(subdir:"unfiltered")] : 
                                                  [extend + bamCoverage.using(subdir:"unfiltered"), 
                                                   phantompeak.using(subdir:"unfiltered")]), 
                	ipstrength.using(subdir:"unfiltered"), 
		        macs2.using(subdir:"unfiltered") + blacklist_filter.using(subdir:"unfiltered")                  
		]
            ], 
            [ filbowtie2unique + BAMindexer +  
               (ESSENTIAL_DEDUPLICATION ? [RmDups + BAMindexer] : [MarkDups + BAMindexer]) + 
                [       // QC specific to paired end (pe) or single end (se) design
	    		(RUN_IN_PAIRED_END_MODE ? [bamCoverage.using(subdir:"filtered"), 
                                                   InsertSize.using(subdir:"filtered")] : 
                                                  [extend + bamCoverage.using(subdir:"filtered"), 
                                                   phantompeak.using(subdir:"filtered")]), 
                	ipstrength.using(subdir:"filtered"), 
		        macs2.using(subdir:"filtered") + blacklist_filter.using(subdir:"filtered")
		]  
            ]
        ] // end parallel branches

      ] + 
    [(RUN_PEAK_ANNOTATION ? peak_annotation.using(subdir:"unfiltered") : dontrun.using(module:"peak_annotation")) +
     (RUN_DIFFBIND ? diffbind.using(subdir:"unfiltered") : dontrun.using(module:"diffbind")) +
     (RUN_ENRICHMENT ? GREAT.using(subdir:"unfiltered") : dontrun.using(module:"GREAT")),
     (RUN_PEAK_ANNOTATION ? peak_annotation.using(subdir:"filtered") : dontrun.using(module:"peak_annotation")) +
     (RUN_DIFFBIND ? diffbind.using(subdir:"filtered") : dontrun.using(module:"diffbind")) +
     (RUN_ENRICHMENT ? GREAT.using(subdir:"filtered") : dontrun.using(module:"GREAT")),
    ] +
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + collectBpipeLogs + MultiQC + 
    shinyReports
}

