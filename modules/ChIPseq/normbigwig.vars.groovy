normbigwig_vars=[
    outdir          : TRACKS + "/input_normalised_cov",
    targets         : "targets.txt", // targets file describing the samples
    extension_length: ESSENTIAL_FRAGLEN-ESSENTIAL_READLEN, //for paired end the pairs are automatically used. The extension length is only used if there are singletons
    mapped          : MAPPED,
    threads         : Integer.toString(ESSENTIAL_THREADS),
    extra           : "--scaleFactorsMethod readCount --operation subtract " + "--extendReads " + Integer.toString(EXTENSION_LENGTH) + " --outFileFormat bedgraph"
]
