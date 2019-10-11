GenerateStarIndexFromSJ_vars=[
    outdir    : MAPPED + "/sjdbStarIndex",
    sjdbfile  : MAPPED + "/sjdbStarIndex/SJ.out.tab.Pass1.sjdb",
    outdir_2nd_index: FilterAndMergeSJtab_vars.outdir,
    threads   : STAR_pe_vars.threads,
    genome_ref: ESSENTIAL_GENOME_REF,
    maxram    : STAR_pe_vars.maxram,
    bufsize   : STAR_pe_vars.bufsize,
    overhang  : STAR_pe_vars.overhang,
    extra     : ""
]
