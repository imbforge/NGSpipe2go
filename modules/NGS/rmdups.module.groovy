//rule for task MarkDups from catalog NGS, version 1
//desc: Mark with/without removing duplicated reads from a bam file
RmDups = {
	doc title: "MarkDups",
		desc:  "Call picard tools to mark with/without removing duplicated reads from a bam file",
		constraints: "Picard tools version >= 1.141",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir=MAPPED
	def JAVA_FLAGS  = "-Xmx" + MARKDUPS_MAXMEM + "m"
	def MARKDUPS_FLAGS  = "REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE"

	transform(".bam") to (".duprm.bam") {
		exec """
			module load java/${JAVA_VERSION} &&
            module load picard/${PICARD_VERSION} &&
			
			java $JAVA_FLAGS -jar ${TOOL_PICARD}/picard.jar MarkDuplicates $MARKDUPS_FLAGS INPUT=$input OUTPUT=$output METRICS_FILE=${input.prefix}_dupmetrics.tsv
		""","RmDups"
	}
}

