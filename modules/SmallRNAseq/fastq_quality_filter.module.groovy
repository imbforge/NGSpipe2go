FastQQualityFilter = {
	doc title: "Remove sequences",
		desc:  "filter reads containing low-quality (Phred score below 20) bases in order to facilitate the PCR duplicates removal.",
		constraints: "Only supports uncompressed FASTQ files",
		author: "Antonio Domingues"

	output.dir = FASTQ_QUALITY_FILTER_OUTDIR
   FASTQ_QUALITY_FILTER_FLAGS= " -q "   + MIN_QUAL  +
                " -p "   + MIN_PERCENT  +
             " -Q " + QUAL_FORMAT +
             FASTQ_QUALITY_FILTER_OTHER

   transform(".cutadapt28to51nt.fastq") to (".highQ.fastq") {
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         fastq_quality_filter -h 1>&2 &&
         echo '/VERSION INFO' 1>&2 &&

         fastq_quality_filter $FASTQ_QUALITY_FILTER_FLAGS -i $input > $output
      ""","FastQQualityFilter"
   }
}
