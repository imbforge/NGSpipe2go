StringTie_vars=[
    outdir  : RESULTS + "/stringtie",
    gtf     : ESSENTIAL_GENESGTF, 
    stranded: ESSENTIAL_STRANDED,
    threads : Integer.toString(ESSENTIAL_THREADS),
    extra   : "-f 0.1 -B -e"
]
