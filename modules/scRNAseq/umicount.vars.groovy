UMICOUNT_OUTDIR=RESULTS + "/umicount"
UMICOUNT_LOG="--verbose=1"
//this assumes that the labeling is done on the bam file e.g. by processing with featureCounts beforehand
//additionaly this is configured to fit marsseq paramers. It might be necessary to add --read-length if you want to ensure that not only the position + UMI is used to deduplicate but also the read-length. For marsseq this opition is not set because we expect reads with the same umi+starting position to be PCR duplicates event if they are of different length
UMICOUNT_PARAM="--per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts" 
UMICOUNT_EXTRA="--edit-distance-threshold=0 " //Spliced reads are treated different from unspliced
