//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
bowtie_se = {
	doc title: "Bowtie SE alignment",
		desc:  "Align single end reads",
		constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
		author: "Sergi Sayols"

	output.dir = MAPPED

	def BOWTIE_FLAGS = " -q --sam"  +
	                   " "   + BOWTIE_QUALS    +
					   " "   + BOWTIE_BEST     +
	                   " -n" + BOWTIE_MM       +
					   " -l" + BOWTIE_INSERT   +
					   " -e" + BOWTIE_MAQERR   +
					   " -m" + BOWTIE_MULTIMAP +
					   " -p" + BOWTIE_THREADS

	transform(".fastq.gz") to (".bam") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_BOWTIE}/env.sh   &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&

			zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - | ${TOOL_SAMTOOLS} view -bhSu - | ${SAMTOOLS} sort -@ $BOWTIE_THREADS - $output.prefix
		""","bowtie_se"
	}
}

