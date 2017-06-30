//vars for task geneBodyCov from catalog RNAseq, version 1
INFEREXPERIMENT_BED=ESSENTIAL_GENESBED // gene model
INFEREXPERIMENT_OUTDIR=QC + "/inferexperiment"
INFEREXPERIMENT_BED="-r " + ESSENTIAL_GENESBED //this variable is essential for the module to run do not set it to the empty string! 
INFEREXPERIMENT_EXTRA="-s 4000000" //add other options here in this case it is the samples size (how many reads should be samples from the bam file) 
