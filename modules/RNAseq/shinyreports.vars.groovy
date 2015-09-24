//vars for task shinyReports from catalog miscellaneous, version 1
SHINYREPS_PROJECT=PROJECT	//project directory
SHINYREPS_ORG=ESSENTIAL_ORG	//UCSC organism
SHINYREPS_DB=ESSENTIAL_DB	//UCSC assembly version
SHINYREPS_LOG=LOGS			//where the logs lie
SHINYREPS_QC=QC				//where the QC lie
SHINYREPS_RES=RESULTS		//where the results lie
SHINYREPS_PREFIX=ESSENTIAL_SAMPLE_PREFIX	//standard sample prefix
SHINYREPS_STAR_LOG=LOGS + "/STAR_se"	//where the STAR logs lie
SHINYREPS_STAR_SUFFIX="Log.final.out"	//extension given to the STAR log files
SHINYREPS_STARparms_SUFFIX="Log.out"	//extension given to the STAR log files
SHINYREPS_FASTQC_LOG=FASTQC_OUTDIR		//where the Fastqc logs lie
SHINYREPS_DUPRADAR_LOG=DUPRADAR_OUTDIR	//where the dupRadar logs lie
SHINYREPS_RNATYPES_LOG=QC + "/RNAtypes"	//where the RNAtypes logs lie
SHINYREPS_GENEBODYCOV_LOG=GENEBODYCOV_OUTDIR //where the geneBodyCov logs lie
SHINYREPS_BUSTARD=QC + "/DemultiplexedBustardSummary.xml"	//where the bustard xml file lies
SHINYREPS_DE_EDGER=DE_edgeR_OUTDIR + "/DE_edgeR.RData"   //where the DE_edgeR output lies
SHINYREPS_SUBREAD=RESULTS + "/subread-count" // location of the subread counts
SHINYREPS_SUBREAD_SUFFIX=".raw_readcounts.tsv.summary" // the extension of the subread stats file
SHINYREPS_PLOTS_COLUMN=4    //number of columns to splits the plot grids (dupradar, genebodycov...). Min=2

