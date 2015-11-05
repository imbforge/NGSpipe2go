//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
bowtie_se = {
	doc title: "Bowtie SE alignment",
		desc:  "Align single end reads",
		constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = MAPPED

	def BOWTIE_FLAGS = " -q --sam"  +
                       " "   + BOWTIE_QUALS    +
                       " "   + BOWTIE_BEST     +
                       " -n" + Integer.toString(BOWTIE_MM)       +
                       " -l" + Integer.toString(BOWTIE_INSERT)   +
                       " -e" + Integer.toString(BOWTIE_MAQERR)   +
                       " -m" + Integer.toString(BOWTIE_MULTIMAP) +
                       " -p" + Integer.toString(BOWTIE_THREADS)

	transform(".fastq.gz") to (".bam") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_BOWTIE}/env.sh   &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(bowtie --version | grep bowtie | cut -d' ' -f3)   1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			
			zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - | ${TOOL_SAMTOOLS} view -bhSu - | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_THREADS - $output.prefix
		""","bowtie_se"
	}
}

