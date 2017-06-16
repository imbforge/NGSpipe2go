//vars for task geneBodyCov2
GENEBODYCOV2_GTF="gtf=" + ESSENTIAL_GENESGTF  // the gencode annotation GTF (can be compressed)
GENEBODYCOV2_PAIRED="paired=" + ESSENTIAL_PAIRED   // paired end yes|no
GENEBODYCOV2_STRANDED="stranded=" + ESSENTIAL_STRANDED // strandness yes|no|reverse
GENEBODYCOV2_OUTDIR="outdir=" + QC + "/geneBodyCov"
GENEBODYCOV2_THREADS="threads=" + Integer.toString(ESSENTIAL_THREADS) // number of cores to use
