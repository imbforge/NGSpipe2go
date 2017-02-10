BOWTIE2_THREADS=" -p 12" //+ Integer.toString(ESSENTIAL_THREADS) // threads to use
BOWTIE2_SAMTOOLS_THREADS=Integer.toString(ESSENTIAL_THREADS) // threads to use
BOWTIE2_REF=" -x " + ESSENTIAL_BOWTIE2_REF // prefix of the bowtie reference genome
BOWTIE2_INSERT=""			// seed size. Match with fragment size, the default 20 is set by --very-sensitive
BOWTIE2_MM_SEED="-N 1"			// number of mismatches allowed in the seed
BOWTIE2_QUALS="--phred33"	// phred33-quals. Use --phred64-quals for old sequencing runs
BOWTIE2_EXTRA="--fr --very-sensitive --end-to-end --maxins 1000 --minins 0"
BOWTIE2_MAPPED=MAPPED
