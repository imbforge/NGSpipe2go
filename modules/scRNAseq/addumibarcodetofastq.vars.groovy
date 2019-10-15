AddUMIBarcodeToFastq_vars=[
    outdir     : PROJECT + "/rawdata_processed",
    bcpattern  : ESSENTIAL_BCPATTERN, // pattern of the umi and the barcode in the second read. The C are the barcode bases the Ns are the UMI bases
    barcodelist: ESSENTIAL_WHITELIST, // list of valid barcodes
    extra      : ""
]
