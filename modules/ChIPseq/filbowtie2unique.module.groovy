filbowtie2unique = {
	doc title: "Bowtie PE alignment",
		desc:  "Align paired end reads",
		constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Nastasja Kreim"

	output.dir = BOWTIE2_MAPPED
	
	transform(".bam") to (".unique.bam") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(samtools --version | grep bowtie | cut -d' ' -f3)   1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			${TOOL_SAMTOOLS} view  ${input} | grep -v XS | samtools view -bhSu -T $BOWTIE2_GENOME - |  ${TOOL_SAMTOOLS} sort -@ $BOWTIE2_SAMTOOLS_THREADS - ${output.prefix}
		""","filbowtie2unique"
	}
}

