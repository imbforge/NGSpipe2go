//vars for task macs2 from catalog ChIPseq, version 1
MACS2_TARGETS="targets.txt" // targets file describing the samples
MACS2_MFOLD="1 30"	// range of enrichment ratio (default: 10,30)
MACS2_GSIZE=ESSENTIAL_MACS2_GSIZE // the mappable genome size
MACS2_BWIDTH=ESSENTIAL_FRAGLEN	  // bandwidth use for model building
MACS2_OTHER=""		// other parms sent to macs2
MACS2_MAPPED=MAPPED // where the bam files are stored

