//rule for task BAMindexer from catalog NGS, version 1
//desc: Samtools index a bam file
BAMindexer = {
	doc title: "BAMindexer",
		desc:  "Call samtools to index a bam file",
		constraints: "Define a global SAMTOOLS var pointing to the bin file",
		author: "Sergi Sayols"

	output.dir = MAPPED

	transform(".bam") to(".bam.bai") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			${TOOL_SAMTOOLS} index $input
		""","BAMindexer"
	}

	forward input
}

