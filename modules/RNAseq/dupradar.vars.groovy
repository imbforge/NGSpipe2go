DUPRADAR_STRANDED="stranded=" + ESSENTIAL_STRANDED  // strandness
DUPRADAR_PAIRED="paired="     + ESSENTIAL_PAIRED    // is a paired end experiment
DUPRADAR_OUTDIR="outdir="     + QC + "/dupRadar"    //output dir. If you change it here, change it in the module file also
DUPRADAR_THREADS="threads="   + Integer.toString(ESSENTIAL_THREADS)	// number of threads to be used
DUPRADAR_GTF="gtf=" + ESSENTIAL_GENESGTF    // gene model
DUPRADAR_EXTRA=""   // extra parms sent to the tool
