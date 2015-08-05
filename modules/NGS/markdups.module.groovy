//rule for task MarkDups from catalog NGS, version 1
//desc: Mark with/without removing duplicated reads from a bam file
MarkDups = {
	doc title: "MarkDups",
		desc:  "Call picard tools to mark with/without removing duplicated reads from a bam file",
		constraints: "",
		author: "Sergi Sayols"

	output.dir=MAPPED
	def MARKDUPS_FLAGS  = "REMOVE_DUPLICATES=" + MARKDUPS_REMOVE + " ASSUME_SORTED=TRUE"

	transform(".bam") to ("_duprm.bam") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_PICARD}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			java -Xmx50000m -jar ${TOOL_PICARD}/MarkDuplicates.jar $MARKDUPS_FLAGS
				INPUT=$input
				OUTPUT=$output
				METRICS_FILE=${input.prefix}_dupmetrics.tsv
		""","MarkDups"
	}
}

