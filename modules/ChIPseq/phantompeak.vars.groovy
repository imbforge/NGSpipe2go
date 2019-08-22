phantompeak_vars=[
    outdir  : QC + "/phantompeak",
    threads : Integer.toString(ESSENTIAL_THREADS), // number of threads to use
    minshift: Integer.toString(-500), // left 'x' coordinate in plot
    maxshift: Integer.toString(1500), // right 'x' coordinate in plot
    binsize : Integer.toString(5),    // stepsize for cc calculation
    readlen : Integer.toString(ESSENTIAL_READLEN), // read length
    extra   : ""                      // extra parms to pass to the tool
]
