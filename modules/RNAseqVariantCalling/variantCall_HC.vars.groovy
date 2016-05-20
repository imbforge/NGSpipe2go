// vars for task BaseRecalibration step from catalog RNAseqVariantCalling, version 1
OUTPUT_HC=RESULTS +  "/HC"
GATK_REF=ESSENTIAL_GENOME_REF
VCF_REF=ESSENTIAL_VCF_REF
GATK_THREADS=ESSENTIAL_THREADS
HC_MAXMEM = "2400m"  // max mem in MB to reserve for the jvm
MIN_SCORE_CALL=20
MIN_SCORE_EMIT=20

