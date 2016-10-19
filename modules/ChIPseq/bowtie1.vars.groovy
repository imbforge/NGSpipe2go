// bowtie parameters with suggested typical defaults
BOWTIE_THREADS=" -p" + Integer.toString(ESSENTIAL_THREADS) // threads to use
BOWTIE_REF=ESSENTIAL_BOWTIE_REF // prefix of the bowtie reference genome
BOWTIE_INSERT="-l 40"			// seed length, the optimum depends on the read length and quality
BOWTIE_MM_SEED="-n 2"			// maximum number of mismatches allowed in the seed sequence
BOWTIE_MAQERR="-e 70"			// maximum permitted total of quality values at all mismatched positions throughout the entire alignment
BOWTIE_MULTIMAP="-m 1"			// discard (-m 1) or keep one random alignment (-M 1) of all reads mapping to multiple locations
BOWTIE_BEST="--best --strata --tryhard --chunkmbs 256"	// bowtie best mapping mode
BOWTIE_QUALS="--phred33-quals"	// phred33-quals. Use --phred64-quals for old sequencing runs
BOWTIE_EXTRA=""                 // extra parms to be passed to bowtie, e.g. for trimming barcodes

// samtools parameters
SAMTOOLS_THREADS="-@" + Integer.toString(ESSENTIAL_THREADS)
