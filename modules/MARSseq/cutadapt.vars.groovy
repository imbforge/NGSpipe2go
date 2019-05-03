CUTADAPT_OUTDIR=TRIMMED
CUTADAPT_ADAPTER_SEQUENCE="--adapter " + ESSENTIAL_ADAPTER_SEQUENCE  // adapter to be trimmed off
CUTADAPT_POLYA="-a A{100}"
CUTADAPT_MINIMUM_OVERLAP="--overlap=" + ESSENTIAL_MINADAPTEROVERLAP  // minimal overlap of the read and the adapter
CUTADAPT_MINIMUM_LENGTH_KEEP="--minimum-length " + ESSENTIAL_MINREADLENGTH  // minimal length of reads to be kept
CUTADAPT_ERRORRATE="--error-rate 0.1" // default 0.1
CUTADAPT_EXTRA="--times 2"
