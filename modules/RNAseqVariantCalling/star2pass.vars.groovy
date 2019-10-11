STAR_pe_2nd_vars=[
    outdir    : MAPPED + "/2ndPass",
    logdir    : LOGS + "/STAR_2ndPass",
    threads   : Integer.toString(ESSENTIAL_THREADS),
    ref       : FilterAndMergeSJtab_vars.outdir,
    maxram    : "31000000000", // around 30Gb for mammals
    bufsize   : "150000000",   // buffer size
    mm        : "2",           // number of mismatches allowed
    multimap  : "10",          // max multimap positions per read
    minintro  : "21",          // minimum intron size
    filter_sec: true,          // filter out secondary alignments from the bam file?
    samtools_threads: Integer.toString(ESSENTIAL_THREADS),
    overhang  : Integer.toString(ESSENTIAL_READLENGTH - 1),
    extra     : ""
]
