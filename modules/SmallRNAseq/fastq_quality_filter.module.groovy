FastQQualityFilter = {
	doc title: "Remove sequences",
		desc:  "filter reads containing low-quality (Phred score below 20) bases in order to facilitate the PCR duplicates removal.",
		constraints: "Only supports compressed FASTQ files",
		author: "Antonio Domingues"

	output.dir = FASTQ_QUALITY_FILTER_OUTDIR

    FASTQ_QUALITY_FILTER_FLAGS = " -q "   + MIN_QUAL  +
                                 " -p "   + MIN_PERCENT  +
                                 " -Q " + QUAL_FORMAT +
                                 FASTQ_QUALITY_FILTER_OTHER

    transform(".cutadapt.fastq.gz") to (".highQ.fastq.gz") {

        def EXP = input.prefix.prefix.prefix

        exec """
            module load fastx_toolkit/${FASTX_VERSION} &&

            zcat $input | fastq_quality_filter $FASTQ_QUALITY_FILTER_FLAGS -o $output 2>&1 >> ${FASTQ_QUALITY_FILTER_OUTDIR}/${EXP}.fastq_quality_filter.log
            
      ""","FastQQualityFilter"
   }
}
