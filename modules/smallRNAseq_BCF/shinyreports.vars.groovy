//vars for task shinyReports from catalog miscellaneous, version 1
SHINYREPS_PROJECT=PROJECT	//project directory
SHINYREPS_MINADAPTEROVERLAP=ESSENTIAL_MINADAPTEROVERLAP
SHINYREPS_MINREADLENGTH=ESSENTIAL_MINREADLENGTH
SHINYREPS_PREFIX=ESSENTIAL_SAMPLE_PREFIX	//standard sample prefix
SHINYREPS_BOWTIE_RES=MAPPED	//where the Bowtie logs lie
SHINYREPS_BOWTIE_SUFFIX=".bowtie.log"
SHINYREPS_BOWTIE_PLOT_DIR=MAPPED + "/plots"
SHINYREPS_BOWTIE_LOG=LOGS + "/Bowtie_se"
SHINYREPS_STAR_SUFFIX="Log.final.out"   //extension given to the STAR log files
SHINYREPS_STARparms_SUFFIX="Log.out"    //extension given to the STAR log files
SHINYREPS_STAR_LOG=LOGS + "/STAR"    //where the STAR logs lie
SHINYREPS_STAR_PLOT_DIR=MAPPED + "/plots"
SHINYREPS_CUTADAPT_PLOT_DIR=PROCESSED + "/plots"
SHINYREPS_CUTADAPT_LOG=LOGS + "/Cutadapt"
SHINYREPS_QUALITYFILTER_PLOT_DIR=PROCESSED + "/plots"
SHINYREPS_QUALITYFILTER_LOG=LOGS + "/FastQQualityFilter"
SHINYREPS_DEDUP_PLOT_DIR=PROCESSED + "/plots"
SHINYREPS_DEDUP_LOG=LOGS + "/FilterDuplicates"
SHINYREPS_RAWFILTERSUMMARY_PLOT_DIR=PROCESSED + "/plots"
SHINYREPS_FASTQC_OUT=FASTQC_OUTDIR		//where the Fastqc output lie
SHINYREPS_FASTQC_LOG=LOGS + "/FastQC"		//where the Fastqc logs lie
SHINYREPS_FASTQSCREEN_OUT=FASTQSCREEN_OUTDIR
SHINYREPS_FASTQSCREEN_LOG=LOGS + "/FastQScreen"
SHINYREPS_BAMINDEX_LOG=LOGS + "/BAMindexer"	//where the Samtools/BamIndexer logs lie
SHINYREPS_RNATYPES_LOG=QC + "/RNAtypes"	//where the RNAtypes logs lie
SHINYREPS_RNATYPES=QC + "/RNAtypes"	//where the RNAtypes count results lie
SHINYREPS_RNATYPES_SUFFIX=".readcounts.tsv" // the extension of the subread results files
SHINYREPS_RNATYPES_CUTOFF=0.005  //all types occurring in a larger fraction than SHINYREPS_RNATYPES_CUTOFF will be plotted individually
SHINYREPS_SUBREAD=RESULTS + "/subread-count" // location of the subread counts
SHINYREPS_SUBREAD_SUFFIX=".raw_readcounts.tsv.summary" // the extension of the subread stats file
SHINYREPS_SUBREAD_LOG=LOGS + "/SubreadCount"	//where the Subread/FeatureCounts logs lie
SHINYREPS_PLOTS_COLUMN=3L    //number of columns to splits the plot grids (dupradar, genebodycov...). Min=2L. L=integer in R
SHINYREPS_TOOL_VERSIONS= PROJECT + "/NGSpipe2go/modules/smallRNAseq_BCF/tool.versions.groovy" //where the tool versions listed

