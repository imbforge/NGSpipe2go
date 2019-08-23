STAR_vars=[
    outdir          : MAPPED,
    logdir          : LOGS + "/STAR",
    mm              : "0.04",          // max. permitted mismatches in the alignment, relative to the read length (4% as used by ENCODE)
    unmapped_bam    : "Within",        // report unmapped reads to bam file? (choices: Within, None)

    // settings imported from essential vars
    threads         : Integer.toString(ESSENTIAL_THREADS),
    ref             : ESSENTIAL_STAR_REF,
    overhang        : Integer.toString(ESSENTIAL_READLENGTH - 1),
    gtf             : ESSENTIAL_GENESGTF,   // gene model GTF file
    paired          : (ESSENTIAL_PAIRED == "yes"),

    // additional settings
    filter_sec      : true,            // remove secondary alignments in order to keep just 1 alignment of the multi-mapping reads ?
    samtools_threads: Integer.toString(ESSENTIAL_THREADS),
    extra           : ""               // extra parms to sent to the tool
]
