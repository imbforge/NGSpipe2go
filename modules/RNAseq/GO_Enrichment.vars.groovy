GO_Enrichment_RDATA ="rData=" + DE_DESeq2_OUTDIR
GO_Enrichment_LOG2FOLD="log2Fold=0" 
GO_Enrichment_PADJ="padj=0.01"
GO_Enrichment_ORG="organism=" + ESSENTIAL_ORG
GO_Enrichment_UNIV="univ=expressed"
GO_Enrichment_TYPE="type=gene_name" 
GO_Enrichment_CATEGORY="plotCategory=5"
GO_Enrichment_OUTDIR="out=" + RESULTS + "/GO_Analysis"
GO_Enrichment_CORES="cores=" + Integer.toString(ESSENTIAL_THREADS) // number of cores to use
GO_Enrichment_EXTRA=""

