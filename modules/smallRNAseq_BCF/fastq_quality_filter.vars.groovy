//vars for task FastQQualityFilter from catalog smallRNAseq_BCF, version 0.1

FASTQ_QUALITY_FILTER_OUTDIR=PROCESSED
MIN_QUAL=20      // minimal quality of bases in reads to be kept
MIN_PERCENT=100  // percentage of bases fulfilling the minimal quality requirement
QUAL_FORMAT=33   // format of the quality scores
FASTQ_QUALITY_FILTER_OTHER=" -v -z"

