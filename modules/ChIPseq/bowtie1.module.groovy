//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
bowtie_se = {
	doc title: "Bowtie SE alignment",
		desc:  "Align single end reads",
		constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = MAPPED

	def BOWTIE_FLAGS = "-q --sam "  +
                       BOWTIE_QUALS    + " " + 
                       BOWTIE_BEST     + " " + 
                       BOWTIE_MM_SEED  + " " + 
                       BOWTIE_INSERT   + " " + 
                       BOWTIE_MAQERR   + " " + 
                       BOWTIE_MULTIMAP + " " + 
                       BOWTIE_THREADS  + " " + 
                       BOWTIE_EXTRA

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
			
			zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - | ${TOOL_SAMTOOLS} view -bhSu - | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_SAMTOOLS_THREADS - $output.prefix
		""","bowtie_se"
	}
}

