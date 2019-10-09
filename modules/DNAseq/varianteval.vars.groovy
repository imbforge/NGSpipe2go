VariantEval_vars=[
    outdir        : QC + "/GATK_varianteval",
    java_flags    : "-Xmx5000m",
    threads       : Integer.toString(ESSENTIAL_THREADS),
    bwa_ref       : ESSENTIAL_BWA_REF,
    known_variants: ESSENTIAL_KNOWN_VARIANTS
]
