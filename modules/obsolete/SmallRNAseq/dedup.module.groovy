load MODULE_FOLDER + "SmallRNAseq/dedup.vars.groovy"

FilterDuplicates = {
    doc title: "Remove Duplicate sequences",
        desc:  """Identify sequences that share the same random barcodes, and removed them. These are most likely PCR duplicates. It takes two steps to accomplish this:
         (i) convert FastQ to tabular (comma separated) format and filter exact duplicates (NNNN-insert-NNNN) which are most likely PCR clones NOTE: FastQ-Sanger quality scores may have "," for Phred=11 in the raw FastQ files; however the "highQ" files don't contain "," which can be used as a field separator.
         (ii) convert the filtered data back to FastQ format. NOTE: the random barcodes are still present, will be removed during mapping.""",
        author: "Antonio Domingues"

   output.dir = REMOVE_DUP_OUTDIR

   // create the log folder if it doesn't exists
   def REMOVE_DUP_LOGDIR = new File( REMOVE_DUP_OUTDIR + "/logs")
   if (!REMOVE_DUP_LOGDIR.exists()) {
          REMOVE_DUP_LOGDIR.mkdirs()
   }

   transform(".fastq.gz") to (".deduped.fastq.gz") {

      def SAMPLE_NAME = input.prefix.prefix

      exec """
         EXP=\$(basename ${SAMPLE_NAME})

         nreads=\$(zcat $input | echo \$((`wc -l`/4))) &&
         echo \$nreads \${EXP}.highQ > ${REMOVE_DUP_LOGDIR}/${EXP}.dedup_stats.txt &&

         zcat $input | paste -d, - - - - | sort -u -t, -k2,2 | tr ',' '\\n' | gzip > $output &&

         ureads=\$(zcat $output | echo \$((`wc -l`/4))) &&
         echo \$ureads ${EXP}.unique >> ${REMOVE_DUP_LOGDIR}/${EXP}.dedup_stats.txt
      ""","FilterDuplicates"
   }
}
