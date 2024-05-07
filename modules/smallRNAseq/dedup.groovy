FilterDuplicates = {
    doc title: "Remove Duplicate sequences",
        desc: """
          Identify sequences that are identical and share the same random barcodes. Removed them.
            These are most likely PCR duplicates. It takes two steps to accomplish this:
              (i) convert FastQ to tabular (comma separated) format and filter exact duplicates
                  (NNNN-insert-NNNN) which are most likely PCR clones. NOTE: FastQ-Sanger quality
                  scores may have ',' for Phred=11 in the raw FastQ files; however the 'highQ' files
                  don't contain ',' which can be used as a field separator.
              (ii) convert the filtered data back to FastQ format. NOTE: the random barcodes are
                   still present, will be removed during mapping.
        """,
        author: "Antonio Domingues, Anke Busch"

    output.dir = FilterDuplicates_vars.outdir

    // create the log folder if it doesn't exists
    def REMOVE_DUP_LOGDIR = new File(FilterDuplicates_vars.logdir)
    if (!REMOVE_DUP_LOGDIR.exists()) {
        REMOVE_DUP_LOGDIR.mkdirs()
    }

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".fastq.gz") to (".dedup.fastq.gz") {
        def SAMPLENAME = input.prefix.prefix
         exec """
             ${PREAMBLE} &&

             EXP=\$(basename ${SAMPLENAME}) &&
             nreads=\$(zcat $input | echo \$((`wc -l`/4))) &&
             echo inclDup \$nreads > ${FilterDuplicates_vars.logdir}/\${EXP}.dedup.log &&

             zcat $input | paste -d, - - - - | sort -u -t, -k2,2 | tr ',' '\\n' | gzip > $output &&
             nreads=\$(zcat $output | echo \$((`wc -l`/4))) &&
             echo exclDup \$nreads >> ${FilterDuplicates_vars.logdir}/\${EXP}.dedup.log

        ""","FilterDuplicates"
    }
}
