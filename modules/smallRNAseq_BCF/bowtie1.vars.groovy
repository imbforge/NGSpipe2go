Bowtie_se_vars=[
    mapped     : MAPPED,           // output dir
    threads    : Integer.toString(ESSENTIAL_THREADS), // threads to use
    samtools_threads: Integer.toString(ESSENTIAL_THREADS),
    ref        : ESSENTIAL_BOWTIE_REF, // prefix of the bowtie reference genome
    mm         : 1,                // number of mismatches allowed
    multireport: 1,                // if a read has more than <int> reportable alignments, one is reported at random.
    best       : true,             // bowtie best mode (implies --best --strata --tryhard). Doesn't apply to PE
    quals      : "--phred33-quals",// phred33-quals. Use --phred64-quals for old sequencing runs
    extra      : ""
]
