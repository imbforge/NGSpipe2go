// vars for task SplitNCigarReads step from catalog RNAseqVariantCalling, version 1
GATK_REF=ESSENTIAL_GENOME_REF
GATK_THREADS=ESSENTIAL_THREADS
SPLITNCIGAR_MAXMEM = "2400m"  // max mem in MB to reserve for the jvm
READ_FILTER_FLAG="ReassignOneMappingQuality"
MAP_Q_FROM_FLAG=255
MAP_Q_TO_FLAG=60
UNSAFE_FLAG="ALLOW_N_CIGAR_READS"
