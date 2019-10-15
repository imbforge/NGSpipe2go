Cutadapt_vars=[
    outdir             : TRIMMED,
    logdir             : LOGS + "/Cutadapt",
    adapter_sequence   : ESSENTIAL_ADAPTER_SEQUENCE,   // adapter to be trimmed off
    polya              : "-a A{100}",
    minimum_overlap    : ESSENTIAL_MINADAPTEROVERLAP,  // minimal overlap of the read and the adapter
    minimum_length_keep: ESSENTIAL_MINREADLENGTH,      // minimal length of reads to be kept
    errorrate          : "0.1",                        // default 0.1
    extra              : "--cut 5 -a T{100} --times 3" // MARS-Seq: included 5 random bases on 5'sequence
]
