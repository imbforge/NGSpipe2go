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
SHINYREPS_FASTQC_OUT=FASTQC_OUTDIR		//where the Fastqc output lie
SHINYREPS_FASTQC_LOG=LOGS + "/FastQC"		//where the Fastqc logs lie
SHINYREPS_BAMINDEX_LOG=LOGS + "/BAMindexer"	//where the Samtools/BamIndexer logs lie
SHINYREPS_DUPRADAR_LOG=QC + "/dupRadar"	//where the dupRadar logs lie
SHINYREPS_RNATYPES_LOG=QC + "/RNAtypes"	//where the RNAtypes logs lie
SHINYREPS_GENEBODYCOV_LOG=QC + "/geneBodyCov"
SHINYREPS_BUSTARD=QC + "/DemultiplexedBustardSummary.xml"	//where the bustard xml file lies
SHINYREPS_DE_EDGER=""       //where the DE_edgeR output lies
SHINYREPS_DE_DESEQ=RESULTS + "/DE_DESeq2/DE_DESeq2.RData"   //where the DE_DESeq2 output lies
SHINYREPS_DE_DESEQ_MM=RESULTS + "/DE_DESeq2_MM/DE_DESeq2.RData"   //where the DE_DESeq2_MM output lies
SHINYREPS_SUBREAD=RESULTS + "/subread-count" // location of the subread counts
SHINYREPS_SUBREAD_SUFFIX=".raw_readcounts.tsv.summary" // the extension of the subread stats file
SHINYREPS_SUBREAD_LOG=LOGS + "/subread_count"	//where the Subread/FeatureCounts logs lie
SHINYREPS_BAM2BW_LOG=LOGS + "/bam2bw"        	//where the Bam2BW logs lie
SHINYREPS_MARKDUPS_LOG=LOGS + "/MarkDups"	//where the picard MarkDuplicates logs lie
SHINYREPS_EDGER_LOGS=""                  //where the DE_edgeR logs lie
SHINYREPS_DESEQ_LOGS=LOGS + "/DE_DESeq2" //where the DE_DESeq2 logs lie
SHINYREPS_PLOTS_COLUMN=4L    //number of columns to splits the plot grids (dupradar, genebodycov...). Min=2L. L=integer in R
SHINYREPS_INFEREXPERIMENT_LOGS=QC + "/inferexperiment" //where the inferexperiment logs lie
SHINYREPS_QUALIMAP_LOGS=QC + "/qualimap" //where the qualimap output files are

