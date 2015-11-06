FilterDuplicates = {
	doc title: "Remove Duplicate sequences",
		desc:  """Identify sequences that share the same random barcodes, and removed them. These are most likely PCR duplicates. It takes two steps to accomplish this:
         (i) convert FastQ to tabular (comma separated) format and filter exact duplicates (NNNN-insert-NNNN) which are most likely PCR clones NOTE: FastQ-Sanger quality scores may have "," for Phred=11 in the raw FastQ files; however the "highQ" files don't contain "," which can be used as a field separator.
         (ii) convert the filtered data back to FastQ format. NOTE: the random barcodes are still present, will be removed during mapping.""",
		author: "Antonio Domingues"

	output.dir = REMOVE_DUP_OUTDIR

	transform(".highQ.fastq") to (".deduped_barcoded.fastq.gz") {

      def SAMPLE = input.prefix.prefix
      exec """
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			echo 'REMOVING DUPLICATES' 1>&2 &&

			cat $input | paste -d, - - - -  > ${SAMPLE}.fastq2table &&
         sort -u -t, -k2,2 ${SAMPLE}.fastq2table > ${SAMPLE}.fastq2table.unique &&
         tr ',' '\n' < ${SAMPLE}.fastq2table.unique > ${SAMPLE}.deduped_barcoded.fastq &&
         gzip ${SAMPLE}.deduped_barcoded.fastq

		""","FilterDuplicates"
	}
}
