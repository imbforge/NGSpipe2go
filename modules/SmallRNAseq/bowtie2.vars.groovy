MULTIMAP_OUT_DIR=BOWTIE2_MAPPED
BOWTIE2_PATH=TOOL_BOWTIE2 // bowtie2 bin to use
BOWTIE2_THREADS=8			// threads to use
BOWTIE2_REF=ESSENTIAL_BOWTIE2_REF // prefix of the bowtie2 reference genome
BOWTIE2_SEED_MM=6             // number of mismatches allowed in seed
BOWTIE2_SEED_LENGTH=10             // length of mapping seed
BOWTIE2_SEED_EXT=20             // consecutive seed extension attempts can "fail" before Bowtie 2 moves on
BOWTIE2_RESEED=3             // maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds
BOWTIE2_GAPS=3             // Disallow gaps within <int> positions of the beginning or end of the read
BOWTIE2_RGAP_PEN="30,30"             // Sets the read gap open (<int1>) and extend (<int2>) penalties
BOWTIE2_FGAP_PEN="30,30"             // Sets the reference gap open (<int1>) and extend (<int2>) penalties
BOWTIE2_MM_PEN="4,4"             // Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers.
BOWTIE2_MIN_SCORE="L,-1,-1"             // Sets a function governing the minimum alignment score needed for an alignment to be considered "valid" (i.e. good enough to report)



