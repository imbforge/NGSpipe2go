bowtie_pe_vars=[
    mapped          : MAPPED,          // output dir
    threads         : Integer.toString(ESSENTIAL_THREADS), // threads to use
    samtools_threads: Integer.toString(ESSENTIAL_THREADS),
    ref             : ESSENTIAL_BOWTIE_REF, // prefix of the bowtie reference genome
    insert          : "40",            // seed size. Match with fragment size
    mm_seed         : "2",             // number of mismatches allowed in the seed
//    maqerr          : "70",            // max sum of quals for -n mismatches
    multimap_mode   : "discard",       // discard (-m) or keep a random alignment (-M) of reads mapping to multiple locations
    multimap        : "1",             // discard (-m 1) or keep one random alignment (-M 1) of all reads mapping to multiple locations
    best            : true,            // bowtie best mode (implies --best --strata --tryhard). Doesn't apply to PE
    quals           : "--phred33-quals",    // phred33-quals. Use --phred64-quals for old sequencing runs
    extra           : "--fr"           // extra parms to be passed to bowtie
]
