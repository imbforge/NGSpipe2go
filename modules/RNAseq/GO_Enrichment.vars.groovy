GO_Enrichment_vars=[
    rdata   : DE_DESeq2_vars.outdir,
    log2fold: "0",
    padj    : "0.01",
    org     : ESSENTIAL_ORG,
    univ    : "expressed",
    type    : "gene_name",
    category: "5",
    outdir  : RESULTS + "/GO_Analysis",
    cores   : Integer.toString(ESSENTIAL_THREADS),
    extra   : ""
]
