DE_edgeR_vars=[
    targets  : "targets.txt",         //targets file. Check the bin directory for the format
    contrasts: "contrasts.txt",       //contrasts file. Check the bin directory for the format
    mmatrix  : "~0+group",            //glm formula with the contrasts to be tested. Check edgeR man
    filter   : true,                  //whether to filter invariant genes
    prefix   : "",                    //prefix to be removed from file names if a sample column is not provided in targets.txt
    suffix   : "",                    //suffix to be removed from file names if a sample column is not provided in targets.txt
    cwd      : "" + RESULTS + "/subread-count", //current working directory
    robust   : false,                 //robust estimation of the negative binomial dispersion
    outdir   : RESULTS + "/DE_edgeR", //output filename base pattern
    gtf      : ESSENTIAL_GENESGTF,
    extra    : ""                     // extra parms to sent to the tool
]
