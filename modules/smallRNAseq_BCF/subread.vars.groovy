//vars for task SubreadCount from catalog smallRNAseq_BCF, version 0.1

SUBREAD_STRANDED=ESSENTIAL_STRANDED //whether the data is from a strand-specific assay (illumina SR: always reverse)
SUBREAD_PAIRED=ESSENTIAL_PAIRED	 //paired end design
SUBREAD_DUPLICATES="--ignoreDup" //set the parm to ignore duplicates
SUBREAD_GENESGTF="-a " + ESSENTIAL_GENESGTF
SUBREAD_CORES="-T " + Integer.toString(ESSENTIAL_THREADS)
SUBREAD_OUTDIR=RESULTS + "/subread-count"
SUBREAD_EXTRA="-M "    // extra parms to sent to the tool (-M also count multi-mapping reads)

