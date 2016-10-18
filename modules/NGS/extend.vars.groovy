//vars for task extend from catalog NGS, version 1
EXTEND_FRAGLEN=ESSENTIAL_FRAGLEN - ESSENTIAL_READLEN	//the average fragment length
EXTEND_SAMTOOLS_THREADS="-@ " + Integer.toString(ESSENTIAL_THREADS)
