filbowtie2unique = {
	doc title: "filter out multimapping reads from bowtie2 out",
		desc:  "filter out multimapping reads from bowtie2 out. output bam file",
		constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Nastasja Kreim"

	output.dir = BOWTIE2_MAPPED
	
	transform(".bam") to (".unique.bam") {
		exec """
			module load samtools/${SAMTOOLS_VERSION} &&

            if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
			
			view  ${input} | grep -v XS | samtools view -bhSu -T $BOWTIE2_GENOME - | sort -@ $BOWTIE2_SAMTOOLS_THREADS -T $TMPDIR/\$(basename $output.prefix) -o ${output}
		""","filbowtie2unique"
	}
}

