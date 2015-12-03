FilterDuplicates = {
	doc title: "Remove Duplicate sequences",
		desc:  """Identify sequences that share the same random barcodes, and removed them. These are most likely PCR duplicates. It takes two steps to accomplish this:
         (i) convert FastQ to tabular (comma separated) format and filter exact duplicates (NNNN-insert-NNNN) which are most likely PCR clones NOTE: FastQ-Sanger quality scores may have "," for Phred=11 in the raw FastQ files; however the "highQ" files don't contain "," which can be used as a field separator.
         (ii) convert the filtered data back to FastQ format. NOTE: the random barcodes are still present, will be removed during mapping.""",
		author: "Antonio Domingues"

	output.dir = REMOVE_DUP_OUTDIR

	transform(".highQ.fastq.gz") to (".deduped_barcoded.fastq.gz") {

      def SAMPLE_NAME = input.prefix.prefix

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi                                          &&

         nreads=\$(zcat $input | echo \$((`wc -l`/4))) &&
         echo \$nreads \$input >> dedup.stats.txt

			zcat $input | paste -d, - - - - | sort -u -t, -k2,2 | tr ',' '\\n' | gzip > $output &&

         nreads=\$(zcat $output | echo \$((`wc -l`/4))) &&
         echo \$nreads \$output >> dedup.stats.txt

		""","FilterDuplicates"
	}
}
