peak_annotation_vars=[
    outdir         : RESULTS + "/Peak_Annotation",
    files          : RESULTS + "/macs2",
    transcript_type: "Bioconductor",
    transcript_db  : ESSENTIAL_TXDB,   // eg, TxDb.Mmusculus.UCSC.mm9.knownGene"
    orgdb          : ESSENTIAL_ANNODB, // eg, org.Mm.eg.db
    regiontss      : "3000",
    extra          : ""
]
