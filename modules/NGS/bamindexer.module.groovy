//rule for task BAMindexer from catalog NGS, version 1
//desc: Samtools index a bam file
BAMindexer = {
	doc title: "BAMindexer",
		desc:  "Call samtools to index a bam file",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = MAPPED

	transform(".bam") to(".bam.bai") {
		exec """
			module load samtools &&
			
			echo 'VERSION INFO'  1>&2 &&
			echo \$(samtools --version | grep samtools | cut -d' ' -f2) 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&
			
			samtools index $input
		""","BAMindexer"
	}

	forward input
}

