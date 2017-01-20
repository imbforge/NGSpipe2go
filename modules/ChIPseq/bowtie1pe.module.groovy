//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
bowtie_pe = {
	doc title: "Bowtie PE alignment",
		desc:  "Align paired end reads",
		constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols modified for paired end by Nastasja Kreim"

	output.dir = BOWTIE_MAPPED
	def OUTPUTFILE = input1
	int path_index = OUTPUTFILE.lastIndexOf("/")
	OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
	println(OUTPUTFILE)
	OUTPUTFILE = (OUTPUTFILE =~ /_read1.fastq.gz/).replaceFirst("")

	def BOWTIE_FLAGS = "-q --sam "  +
                       BOWTIE_QUALS    + " " + 
                       BOWTIE_BEST     + " " + 
                       BOWTIE_MM_SEED  + " " + 
                       BOWTIE_INSERT   + " " + 
                       BOWTIE_MAQERR   + " " + 
                       BOWTIE_MULTIMAP + " " + 
                       BOWTIE_THREADS  + " " + 
                       BOWTIE_EXTRA

	produce(OUTPUTFILE + ".bam") {
		exec """
			module load bowtie &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
			echo $output;
			echo 'VERSION INFO'  1>&2 ;
			echo \$(bowtie --version | grep bowtie | cut -d' ' -f3)   1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			base1=`basename $input1`;
			base2=`basename $input2`;
			zcat $input1 > \$TMPDIR/\${base1%.gz};
			zcat $input2 > \$TMPDIR/\${base2%.gz};
			
			bowtie $BOWTIE_FLAGS $BOWTIE_REF -1 $TMPDIR/\${base1%.gz} -2 $TMPDIR/\${base2%.gz} | ${TOOL_SAMTOOLS} view -bhSu - | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_SAMTOOLS_THREADS - $output.prefix
		""","bowtie_pe"
	}
}

