MPSprofiling_vars=[
    targets  : "targets.txt",              // comma-sep targets file. Needs columns "pruned_file_name", "unique_sample_id", "experiment", "sub_experiment", "bin", "fraction" and optionally "barcode_demultiplex"
    prefix   : "",                         // prefix to be removed from file names
    suffix   : "",                         // suffix to be removed from file names
    inputdir : RESULTS + "/barcode_count", // directory where the .tsv files with the barcode counts are located
    outdir   : RESULTS + "/MPSprofiling",  // output directory
    logdir   : LOGS + "/MPSprofiling",     // log directory
    expdesign : ESSENTIAL_EXPDESIGN,       // experimental design
    threshold_rel_countssum : 0.001,       // lower threshold for first qc filtering (samples with sample sum < threshold_rel_countssum * mean sample sum are removed). 0 means no exclusion.
    excludeSeqsNotInAllFractions : "FALSE",  // discard all sequences per sub_experiment, which were not detected in all fractions of this sub_experiment
    extra    : ""                          // extra parms to sent to the tool
]
