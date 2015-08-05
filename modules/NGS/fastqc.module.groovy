FastQC = {
	doc title: "FastQC",
		desc:  "Quality control of input file",
		constraints: "Only supports compressed FASTQ files",
		author: "Sergi Sayols"

	output.dir   = FASTQC_OUTDIR
	def FASTQC_FLAGS = "--extract --quiet -t " + FASTQC_THREADS
	
	transform(".fastq.gz") to ("_fastqc.zip") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_FASTQC}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&
			
			fastqc $FASTQC_FLAGS -o $output.dir $input
		""","FastQC"
	}

	forward input
}
