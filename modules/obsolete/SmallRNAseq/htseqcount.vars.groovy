HTSEQCOUNT_STRANDED=ESSENTIAL_STRANDED//whether the data is from a strand-specific assay (illumina SR: always reverse)
HTSEQCOUNT_FILE="-f bam"   //paired end design
HTSEQCOUNT_MODE="-m intersection-nonempty" //set the parm to ignore duplicates
HTSEQCOUNT_GENESGTF=ESSENTIAL_GENESGTF
HTSEQCOUNT_TRANSPSONSGTF=ESSENTIAL_TRANSPOSONS
HTSEQCOUNT_OUTDIR=RESULTS + "/htseq-count"
HTSEQCOUNT_EXTRA=""    // extra parms to sent to the tool
