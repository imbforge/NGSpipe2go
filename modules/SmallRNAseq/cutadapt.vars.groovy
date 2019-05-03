//vars for task Cutadapt from catalog SmallRNAseq, version 0.1

CUTADAPT_OUTDIR=RESULTS + "/processed_reads"
ADAPTER_SEQUENCE="-a " + ESSENTIAL_ADAPTER_SEQUENCE // adapter to be trimmed off
MINIMUM_OVERLAP=5
MINIMUM_LENGTH_KEEP=MIN_LENGTH
MAXIMUM_LENGTH_KEEP=MAX_LENGTH
