filter2htseq = {
	doc title: "filter_featureCounts_to_htseq",
	desc: "filter featureCount output to fit HTSeq format, extract column 1 and 7 as well as skipping the header",
	constraints: "none.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"
	
	output.dir = SUBREAD_OUTDIR
	
	transform(".raw_readcounts.tsv") to (".readcounts.tsv") {
		exec """
		tail -n +3 $input | awk '{print \$1\"\\t\"\$7}' > $output
		"""
	}
	
}
