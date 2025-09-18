PIPELINE="scRNAseq"
PIPELINE_VERSION="2.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

load PIPELINE_ROOT + "/pipelines/scRNAseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/scRNAseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.header"
load PIPELINE_ROOT + "/modules/RNAseq/inferexperiment.header"
load PIPELINE_ROOT + "/modules/RNAseq/qualimap.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread.header"
load PIPELINE_ROOT + "/modules/RNAseq/filter2htseq.header"
load PIPELINE_ROOT + "/modules/RNAseq/subread2rnatypes.header"
load PIPELINE_ROOT + "/modules/scRNAseq/cellranger_count.header"
load PIPELINE_ROOT + "/modules/scRNAseq/cellranger_aggr.header"
load PIPELINE_ROOT + "/modules/scRNAseq/splitpipe_all.header"
load PIPELINE_ROOT + "/modules/scRNAseq/splitpipe_comb.header"
load PIPELINE_ROOT + "/modules/scRNAseq/fwd_ext_bams.header"
load PIPELINE_ROOT + "/modules/scRNAseq/demux_hto.header"
load PIPELINE_ROOT + "/modules/scRNAseq/demux_gt.header"
load PIPELINE_ROOT + "/modules/scRNAseq/assignSouporcellCluster.header"
load PIPELINE_ROOT + "/modules/scRNAseq/sc_star.header"
load PIPELINE_ROOT + "/modules/scRNAseq/readAggrData_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/qc_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/filtCells_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/norm_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/igraph_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/hclust_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/kmeans_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/collectCl_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/findmarkers_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/scType_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/collectCT_bioc.header"
load PIPELINE_ROOT + "/modules/scRNAseq/de_bioc.header"
load PIPELINE_ROOT + "/modules/NGS/bamcoverage.header"
load PIPELINE_ROOT + "/modules/NGS/bamindexer.header"
load PIPELINE_ROOT + "/modules/NGS/fastqc.header"
load PIPELINE_ROOT + "/modules/NGS/fastqscreen.header"
load PIPELINE_ROOT + "/modules/NGS/markdups2.header"
load PIPELINE_ROOT + "/modules/NGS/insertsize.header"
load PIPELINE_ROOT + "/modules/NGS/cutadapt.header"
load PIPELINE_ROOT + "/modules/miscellaneous/collect_tool_versions.header"
load PIPELINE_ROOT + "/modules/scRNAseq/shinyreports.header"
load PIPELINE_ROOT + "/modules/NGS/multiqc.header"


dontrun = { println "didn't run $module" }
collect_bams = { forward inputs.bam }

Bpipe.run { 
    (ESSENTIAL_SEQTYPE == "ScaleBio" ? fwd_ext_dir :  
        ("%.fastq.gz" * [ FastQC + 
        (RUN_FASTQSCREEN ? FastqScreen : dontrun.using(module: "FastqScreen")), 
        (RUN_CUTADAPT ? Cutadapt + FastQC.using(subdir:"trimmed") : dontrun.using(module:"Cutadapt")) ])
        ) + 
	  
    (ESSENTIAL_SEQTYPE == "tenX" ? "%_L*_R*_001.fastq.gz" * [ cellranger_count ] : 
		(ESSENTIAL_SEQTYPE == "ParseBio" ? "%.R*.fastq.gz" * [ splitpipe_all ] : 
			(ESSENTIAL_SEQTYPE == "ScaleBio" ? fwd_ext_bams : 
				(ESSENTIAL_SEQTYPE == "SmartSeq" ? ((RUN_IN_PAIRED_END_MODE ? "%.R*.fastq.gz" : "%.fastq.gz") * [ STAR + BAMindexer + subread_count + filter2htseq ]) : 
					dontrun.using(module:"Mapping module. Sequencing type not found!") )))) + 
				
	collect_bams + "%.bam" * [

            (RUN_DEMUX ? (["demux_GT","demux_GT_noAssignment"].contains(RUN_DEMUX) ? demux_gt : (RUN_DEMUX == "demux_HTO" ? demux_hto : dontrun.using(module:"demux") )) : dontrun.using(module:"demux")),
            bamCoverage,
            inferexperiment,
            qualimap,
            subread2rnatypes,
            geneBodyCov2
    ] + 
    (RUN_DEMUX == "demux_GT" ? assignSouporcellCluster : dontrun.using(module:"assignSouporcellCluster")) +

	(ESSENTIAL_USE_AGGR_DATA ? [   
	    (ESSENTIAL_SEQTYPE == "tenX" ? [ cellranger_aggr + readAggrData_bioc] : 
		    (ESSENTIAL_SEQTYPE == "ParseBio" ? [ splitpipe_comb + readAggrData_bioc] : 
			    (ESSENTIAL_SEQTYPE == "ScaleBio" ? [ readAggrData_bioc ] : 
    			    (ESSENTIAL_SEQTYPE == "SmartSeq" ? [ readAggrData_bioc ] : 
    				    dontrun.using(module:"Sample aggregation. Sequencing type not found!") )))) ] :
    dontrun.using(module:"module readIndivSamplesAndMerge_bioc not implemented yet!") ) +

    qc_bioc + filtCells_bioc + norm_bioc + 
    [igraph_bioc, hclust_bioc, kmeans_bioc] + collectCl_bioc + findmarkers_bioc + 
    [scType_bioc] + collectCT_bioc + de_bioc +
	
    (RUN_TRACKHUB ? trackhub_config + trackhub : dontrun.using(module:"trackhub")) +
    collectToolVersions + MultiQC + 
    shinyReports
}

