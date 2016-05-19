FastQQualityFilter = {
	doc title: "Remove sequences",
		desc:  "filter reads containing low-quality (Phred score below 20) bases in order to facilitate the PCR duplicates removal.",
		constraints: "Only supports compressed FASTQ files",
		author: "Antonio Domingues"

	output.dir = FASTQ_QUALITY_FILTER_OUTDIR
   FASTQ_QUALITY_FILTER_FLAGS= " -q "   + MIN_QUAL  +
                " -p "   + MIN_PERCENT  +
             " -Q " + QUAL_FORMAT +
             FASTQ_QUALITY_FILTER_OTHER

   transform(".cutadapt28to51nt.fastq.gz") to (".highQ.fastq.gz") {
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_FASTX}/fastq_quality_filter -h 2>&1 | grep 'Toolkit') 1>&2 &&
         echo '/VERSION INFO'  1>&2 &&

         zcat $input | ${TOOL_FASTX}/fastq_quality_filter $FASTQ_QUALITY_FILTER_FLAGS -o $output
      ""","FastQQualityFilter"
   }
}
