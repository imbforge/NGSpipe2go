Cutadapt_vars=[
    outdir             : TRIMMED,
    logdir             : LOGS + "/Cutadapt",
    paired             : RUN_IN_PAIRED_END_MODE,       // run in se or pe mode    
    adapter_sequence   : ESSENTIAL_ADAPTER_SEQUENCE,   // adapter to be trimmed off the 3' end (paired data: of the first read)
    polya              : " -a A{100} ",                // trim polyA stretches off the 3' end (paired data: of the first read)
    minimum_overlap    : ESSENTIAL_MINADAPTEROVERLAP,  // minimal overlap of the read and the adapter
    minimum_length_keep: ESSENTIAL_MINREADLENGTH,      // minimal length of reads to be kept
    errorrate          : "0.1",                        // default 0.1
    extra              : "-a T{100} --times 3" // MARS-Seq: included 5 random bases on 5'sequence --cut 5
]

// additional adapter settings which may be used in 'extra' (replace "ADAPTER" with the actual adapter sequence) 
// -a ADAPTER     // Sequence of an adapter ligated to the 3' end (additional to the one given in ESSENTIAL_ADAPTER_SEQUENCE)
// -g ADAPTER     // Sequence of an adapter ligated to the 5' end
// -A ADAPTER     // 3' adapter to be removed from second read in a pair
// -G ADAPTER     // 5' adapter to be removed from second read in a pair
