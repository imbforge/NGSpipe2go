//vars for task geneBodyCov2
GENEBODYCOV_GTF="gtf=" + ESSENTIAL_GENESGTF  // the gencode annotation GTF (can be compressed)
GENEBODYCOV_PAIRED="paired=" + ESSENTIAL_PAIRED   // paired end yes|no
GENEBODYCOV_STRANDED="stranded=" + ESSENTIAL_STRANDED // strandness yes|no|reverse
GENEBODYCOV_OUTDIR="outdir=" + QC + "/geneBodyCov"
GENEBODYCOV_THREADS="threads=" + Integer.toString(ESSENTIAL_THREADS) // number of cores to use
