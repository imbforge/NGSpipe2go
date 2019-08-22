GREAT_vars=[
    outdir    : RESULTS + "/GREAT_analysis",
    files     : RESULTS + "/macs2",
    targets   : "targets.txt", // targets file describing the samples
    padj      : "0.01",
    nterms    : "5",
    db        : ESSENTIAL_DB,
    upstream  : "5",    // 5 kb upstream of the TSS 
    downstream: "1",    // 1 kb downstream of the TSS 
    extra     : ""
]
