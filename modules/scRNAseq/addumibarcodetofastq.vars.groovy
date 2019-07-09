ADDUMIBARCODE_OUTDIR=PROJECT + "/rawdata_processed"
ADDUMIBARCODE_BCPATTERN="--bc-pattern=" + ESSENTIAL_BCPATTERN //pattern of the umi and the barcode in the second read. The C are the barcode bases the Ns are the UMI bases
ADDUMIBARCODE_BARCODELIST="--whitelist=" + ESSENTIAL_WHITELIST + " --filter-cell-barcode" // list of valid barcodes these options do only make sense together
ADDUMIBARCODE_EXTRA=""
