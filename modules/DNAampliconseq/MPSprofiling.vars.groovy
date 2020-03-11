MPSprofiling_vars=[
    targets  : "targets.txt",              //targets file. Check the bin directory for the format
    prefix   : "",                         //prefix to be removed from file names
    suffix   : "",                         //suffix to be removed from file names
    inputdir : RESULTS + "/barcode_count", //directory where the .tsv files with the barcode counts are located
    outdir   : RESULTS + "/MPSprofiling",  //output directory
    expdesign : ESSENTIAL_EXPDESIGN,        // experimental design
    extra    : ""                         //extra parms to sent to the tool
]
