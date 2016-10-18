//vars for task macs2 from catalog ChIPseq, version 1
NORMBIGWIG_TARGETS=ESSENTIAL_MACS_TARGETS // targets file describing the samples
NORMBIGWIG_OUTDIR=TRACKS + "/input_normalised_cov" 
EXTENSION_LENGTH=ESSENTIAL_FRAGLEN-ESSENTIAL_READLEN
//for paired end the pairs are automatically used the extension length is only used if there are singletons
NORMBIGWIG_OTHER="--scaleFactorsMethod readCount --ratio subtract --normalizeUsingRPKM " + "--extendReads " + Integer.toString(EXTENSION_LENGTH) + " --outFileFormat bedgraph"
NORMBIGWIG_MAPPED=MAPPED
NORMBIGWIG_THREADS=Integer.toString(ESSENTIAL_THREADS) 
