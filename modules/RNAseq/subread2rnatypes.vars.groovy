subread2rnatypes_vars=[
    outdir    : QC + "/RNAtypes",
    stranded  : ESSENTIAL_STRANDED,              //whether the data is from a strand-specific assay (illumina SR: always reverse)
    paired    : (ESSENTIAL_PAIRED == "yes"),     //paired end design
    genesgtf  : ESSENTIAL_GENESGTF,
    feature   : "exon",                          // type of feature that is to be counted in
    accumulate: ESSENTIAL_FEATURETYPE,           // type of annotation counts should be accumulated on. Usually that would be gene_id, but in this case we choose gene_biotype
    threads   : Integer.toString(ESSENTIAL_THREADS),
    extra     : ""                               // extra parms to sent to the tool
]
