//rule for task Insertsize from catalog NGS
//desc: get the insert size from a bam paired end bam file
InsertSize = {
	doc title: "InsertSize",
		desc:  "Call picard tools create insert size values",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Nastasja Kreim"

	output.dir=INSERTSIZE_OUTDIR
	def INSERTSIZE_FLAGS = INSERTSIZE_OTHER
	def JAVA_FLAGS  = INSERTSIZE_MAXMEM 

	transform(".bam") to ("_insertsizemetrics.tsv") {
		exec """
        module load R/${R_VERSION} &&
		    module load picard/${PICARD_VERSION} && 
        java $JAVA_FLAGS -jar ${TOOL_PICARD}/picard.jar CollectInsertSizeMetrics $INSERTSIZE_FLAGS INPUT=$input OUTPUT=$output HISTOGRAM_FILE=${output.prefix}_hist.pdf
		""","InsertSize"
	}
}

