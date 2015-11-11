BOWTIE_THREADS=ESSENTIAL_THREADS // threads to use
BOWTIE_REF=ESSENTIAL_BOWTIE_REF  // prefix of the bowtie reference genome
BOWTIE_INSERT=28			// seed size. Match with fragment size
BOWTIE_MM_SEED=2			// number of mismatches allowed in the seed
BOWTIE_MAQERR=70			// max sum of quals for -n mismatches
BOWTIE_MULTIMAP=1			// discard reads mapping to more than MULTIMAP positions
BOWTIE_BEST="--best --strata --chunkmbs 256"	// bowtie best mode
BOWTIE_QUALS="--phred33-quals"	// phred33-quals. Use --phred64-quals for old sequencing runs

