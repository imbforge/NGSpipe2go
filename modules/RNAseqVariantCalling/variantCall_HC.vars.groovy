VariantCallHC_vars=[
    outdir        : RESULTS + "/HC",
    java_flags    : "-Xmx2400m",
    gatk_ref      : ESSENTIAL_GENOME_REF,
    vcf_ref       : ESSENTIAL_VCF_REF,
    threads       : ESSENTIAL_THREADS,
    min_score_call: 20,
    min_score_emit: 20
]
