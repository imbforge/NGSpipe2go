FastxTrimmer = {
	doc title: "Read trimming",
		desc:  "Trims reads up to a given position. Useful when input files with the same read length are needed.",
		constraints: "Only supports uncompressed FASTQ files",
		author: "Antonio Domingues"

	output.dir = FASTQ_TRIMMER_OUTDIR
   FASTQ_TRIMMER_FLAGS= " -f "   + FIRST_BASE  +
                        " -l "   + LAST_BASE  +
                        FASTQ_TRIMMER_OTHER

   transform(".fastq.gz") to (".23bp.fastq.gz") {
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_FASTX}/fastx_trimmer -h 2>&1 | grep 'Toolkit') 1>&2 &&
         echo '/VERSION INFO'  1>&2 &&

         zcat $input | ${TOOL_FASTX}/fastx_trimmer $FASTQ_TRIMMER_FLAGS -o $output
      ""","FastxTrimmer"
   }
}
