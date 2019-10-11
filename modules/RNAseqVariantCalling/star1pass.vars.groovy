STAR_pe_vars=[
    outdir    : MAPPED + "/1stPass",
    logdir    : LOGS + "/STAR_1stPass",
    threads   : Integer.toString(ESSENTIAL_THREADS),
    ref       : ESSENTIAL_STAR_REF,
    maxram    : "31000000000", // around 30Gb for mammals
    bufsize   : "150000000",   // buffer size
    mm        : "2",           // number of mismatches allowed
    multimap  : "10",          // max multimap positions per read
    minintro  : "21",          // minimum intron size
    overhang  : Integer.toString(ESSENTIAL_READLENGTH - 1),
    extra     : ""
]
