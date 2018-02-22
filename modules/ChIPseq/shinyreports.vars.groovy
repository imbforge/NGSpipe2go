//vars for task shinyReports from catalog ChIPseq, version 1
SHINYREPS_PROJECT=PROJECT	//project directory
SHINYREPS_LOG=LOGS			//where the logs lie
SHINYREPS_QC=QC				//where the QC lie
SHINYREPS_RES=RESULTS		//where the results lie
SHINYREPS_PREFIX="Sample_imb_[A-Za-z]+_\\\\d+_\\\\d+_"  //standard sample prefix
SHINYREPS_TARGETS="targets.txt"
SHINYREPS_BOWTIE_LOG=LOGS + "/bowtie_se"	//where the Bowtie logs lie
SHINYREPS_BAMINDEX_LOG=LOGS + "/BAMindexer"	//where the Samtools/BamIndexer logs lie
SHINYREPS_MARKDUPS_LOG=LOGS + "/MarkDups"	//where the MarkDups logs lie
SHINYREPS_EXTEND_LOG=LOGS + "/extend"	//where the extend/BedTools logs lie
SHINYREPS_FASTQC=FASTQC_OUTDIR		        //where the Fastqc logs lie
SHINYREPS_FASTQC_LOG=LOGS + "/FastQC"  //where the Fastqc logs lie
SHINYREPS_IPSTRENGTH=QC + "/ipstrength"	//where the IPstrength files lie
SHINYREPS_IPSTRENGTH_LOG=LOGS + "/ipstrength" //where the IPstrength/R logs lie
SHINYREPS_PBC=QC + "/pbc"               	//where the PBC files lie
SHINYREPS_PHANTOMPEAK=QC + "/phantompeak"	//where the PhantomPeak files lie
SHINYREPS_PHANTOM_LOG=LOGS + "/phantompeak" //where the PhantomPeak/R logs lie
SHINYREPS_BUSTARD=QC + "/DemultiplexedBustardSummary.xml"	//where the bustard xml file lies
SHINYREPS_MACS2=RESULTS + "/macs2"	//where the MACS2 results lie
SHINYREPS_MACS2_LOG=LOGS + "/macs2" //where the macs2 logs lie
SHINYREPS_BLACKLIST_FILTER=RESULTS + "/macs2/peaks_detected_table.csv" //where the BlackList Filter module file lie
SHINYREPS_PLOTS_COLUMN=4            //number of columns to splits the plot grids (ipstrength, phantompeaks...). Min=2
SHINYREPS_PEAK_ANNOTATION=RESULTS + "/Peak_Annotation" // where the peak annotation results lie
SHINYREPS_GREAT=RESULTS + "/GREAT_analysis" // where the GO enrichment results lie
SHINYREPS_TOOL_VERSIONS= PROJECT + "/NGSpipe2go/modules/ChIPseq/tool.versions.groovy" //where the tool versions listed
