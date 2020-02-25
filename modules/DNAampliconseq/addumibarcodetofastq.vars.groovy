AddUMIBarcodeToFastq_vars=[
    outdir     : PROJECT + "/rawdata_processed",
    logdir  : LOGS + "/AddUMIBarcodeToFastq",
    extractmethod: ESSENTIAL_EXTRACTMETHOD, // either "string" or "regex"
    bcpattern  : ESSENTIAL_BCPATTERN, // pattern of the umi and the barcode in the second read. The C are the barcode bases the Ns are the UMI bases
    bcpattern2 : ESSENTIAL_BCPATTERN_2, // in case a 2nd pattern is needed for paired end approach
    barcodelist: ESSENTIAL_WHITELIST, // list of valid barcodes
    extra      : ""
]
