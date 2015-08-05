//vars for task STAR_se from catalog RNAseq, version 1
STAR_THREADS = 8
STAR_REF = ESSENTIAL_STAR_REF
STAR_MAXRAM = "31000000000"	// around 30Gb for mammals
STAR_BUFSIZE = "150000000"	// buffer size
STAR_MM = "2"				// number of mismatches allowed
STAR_MULTIMAP = "10"		// max multimap positions per read
STAR_MININTRO = "21"		// minimum intron size
STAR_FILTER_SEC="YES"		// filter out secondary alignments from the bam file?

