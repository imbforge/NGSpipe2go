BOWTIE2_THREADS=" -p " + Integer.toString(ESSENTIAL_THREADS) // threads to use
BOWTIE2_SAMTOOLS_THREADS=" -@ " + Integer.toString(ESSENTIAL_THREADS) // threads to use
BOWTIE2_REF=" -x " + ESSENTIAL_BOWTIE_REF // prefix of the bowtie reference genome
BOWTIE2_INSERT=""			// seed size. Match with fragment size, the default 20 is set by --very-sensitive
BOWTIE2_MM_SEED="--very-sensitive"	// for very polymorphic samples one can try the "super-sensitive" but >2x slower parameter combination "-D 20 -R 3 -N 1 -L 20 -i S,1,0.50"
BOWTIE2_QUALS="--phred33"	// phred33-quals. Use --phred64-quals for old sequencing runs
BOWTIE2_EXTRA="--fr --end-to-end --maxins 1000 --minins 0"  // mates align fw/rev, entire, read must align, maximum and minimum fragmnet length
BOWTIE2_MAPPED=MAPPED
