RNATYPES_STRANDED=ESSENTIAL_STRANDED//whether the data is from a strand-specific assay (illumina SR: always reverse)
RNATYPES_PAIRED=ESSENTIAL_PAIRED	 //paired end design
RNATYPES_GENESGTF="-a " + ESSENTIAL_GENESGTF
RNATYPES_FEATURE="-t exon" // type of feature that is to be counted in
RNATYPES_ACCUMULATE="-g " + ESSENTIAL_FEATURETYPE // type of annotation counts should be accumulated on. Usually that would be gene_id, but in this case we choose gene_biotype
RNATYPES_CORES="-T " + Integer.toString(ESSENTIAL_THREADS)
RNATYPES_OUTDIR=QC + "/rnatypes-count"

RNATYPES_EXTRA=""    // extra parms to sent to the tool

