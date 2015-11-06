FastQC = {
	doc title: "FastQC",
		desc:  "Quality control of input file",
		constraints: "Only supports compressed FASTQ files",
		author: "Sergi Sayols"

	output.dir   = FASTQC_OUTDIR
	def FASTQC_FLAGS = "--extract --quiet -t " + FASTQC_THREADS

	transform(".fastq.gz") to ("_fastqc.zip") {
		exec """
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			echo 'VERSION INFO'  1>&2 &&
			fastqc --version 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&

			fastqc $FASTQC_FLAGS -o $output.dir $input
		""","FastQC"
	}

	forward input
}
