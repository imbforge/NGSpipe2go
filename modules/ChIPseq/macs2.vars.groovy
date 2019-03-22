// vars for task macs2 from catalog ChIPseq
MACS2_TARGETS="targets.txt" // targets file describing the samples
MACS2_GSIZE="-g " + ESSENTIAL_MACS2_GSIZE // the mappable genome size
MACS2_BWIDTH="--bw " + Integer.toString(ESSENTIAL_FRAGLEN) // bandwidth for model building with SR data
MACS2_MINLEN="--min-length " + Integer.toString(ESSENTIAL_FRAGLEN) // minimum peak length
MACS2_MAPPED=MAPPED // where the bam files are stored
MACS2_EXTRA="--keep-dup " + ESSENTIAL_DUP		// other parameters sent to macs2
MACS2_PAIRED=ESSENTIAL_PAIRED // for PE data use fragments in peak calling
