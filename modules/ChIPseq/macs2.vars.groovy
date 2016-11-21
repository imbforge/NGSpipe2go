//vars for task macs2 from catalog ChIPseq, version 1
MACS2_TARGETS="targets.txt" // targets file describing the samples
MACS2_MFOLD="-m 5,50"	// range of enrichment ratio (default: 5,50)
MACS2_GSIZE="-g " + ESSENTIAL_MACS2_GSIZE // the mappable genome size
MACS2_BWIDTH="--bw " + Integer.toString(ESSENTIAL_FRAGLEN)	  // bandwidth use for model building
MACS2_MAPPED=MAPPED // where the bam files are stored
MACS2_EXTRA=""		// other parms sent to macs2

