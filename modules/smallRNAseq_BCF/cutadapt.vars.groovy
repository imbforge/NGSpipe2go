//vars for task Cutadapt from catalog smallRNAseq_BCF, version 0.1

CUTADAPT_OUTDIR=PROCESSED
ADAPTER_SEQUENCE="-a TGGAATTCTCGGGTGCCAAGG"  // adapter to be trimmed off
MINIMUM_OVERLAP=ESSENTIAL_MINADAPTEROVERLAP  // minimal overlap of the read and the adapter
MINIMUM_LENGTH_KEEP=ESSENTIAL_MINREADLENGTH  // minimal length of reads to be kept
MAXIMUM_LENGTH_KEEP=Integer.toString(ESSENTIAL_READLENGTH - ESSENTIAL_MINADAPTEROVERLAP) // maximal length of reads to be kept for further analysis
MAXIMUM_LENGTH_KEEP_PLUS1=Integer.toString(ESSENTIAL_READLENGTH - ESSENTIAL_MINADAPTEROVERLAP + 1) // minimal length of discarded reads, but they are kept in a separate fastq.gz file
