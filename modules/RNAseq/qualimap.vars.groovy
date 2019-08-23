qualimap_vars=[
    outdir  : QC + "/qualimap",
    stranded: ESSENTIAL_STRANDED, // options are no/yes/reverse 
    genesgtf: ESSENTIAL_GENESGTF,
    paired  : (ESSENTIAL_PAIRED == "yes"),
    extra   : "--java-mem-size=10G"
]
