FastQQualityFilter_vars=[
    outdir     : PROCESSED,
    logdir     : PROCESSED + "/logs",
    min_qual   : 20,  // minimal quality of bases in reads to be kept
    min_percent: 100, // percentage of bases fulfilling the minimal quality requirement
    qual_format: 33,  // format of the quality scores
    extra      : "-v -z"
]
