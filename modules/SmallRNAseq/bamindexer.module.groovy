//rule for task BAMindexer from catalog NGS, version 1
//desc: Samtools index a bam file
BAMindexer = {
	doc title: "BAMindexer",
		desc:  "Call samtools to index a bam file",
		author: "Sergi Sayols"

	output.dir = MAPPED

	transform(".bam") to(".bam.bai") {
		exec """

			echo 'VERSION INFO'  1>&2 &&
			samtools --version 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&

			samtools index $input
		""","BAMindexer"
	}

	forward input
}

