Cutadapt = {
	doc title: "Cutadapt",
		desc:  "remove adapter from reads",
		constraints: "Only supports compressed FASTQ files",
		author: "Antonio Domingues"

	output.dir = CUTADAPT_OUTDIR
   MACS2_FLAGS= " -m "   + MACS2_MFOLD  +
                " -g "   + MACS2_GSIZE  +
                " --bw " + MACS2_BWIDTH +
                MACS2_OTHER

	transform(".fastq.gz") to (".cutadapt28to51nt.fastq") {
		exec """
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			echo 'VERSION INFO'  1>&2 &&
			cutadapt --version 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&

			cutadapt $ADAPTER_SEQUENCE -O $MINIMUM_OVERLAP -m $MINIMUM_LENGTH_KEEP -M $MAXIMUM_LENGTH_KEEP -o $output $input  2>&1 >> cutadapt.log
		""","Cutadapt"
	}
}