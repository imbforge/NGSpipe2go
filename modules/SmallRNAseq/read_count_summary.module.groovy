CountReadsSummary = {
   doc title: "Count Reads Summary",
      desc:  "Aggregates the read counts and plots results. Might be modified in the future.",
      constraints: "The bed files should contains an identifier, for instance, LINE or SINE, in col7.",
      author: "Antonio Domingues"


   output.dir = COUNT_SUMMARY_OUTDIR

   from ("*.counts") produce(
      COUNT_READS_OUTDIR + "/piRNA_quantification.RData"
      ) {


      exec """
         module load R/${R_VERSION} &&

         Rscript $AGREGATE_SCRIPT -i $inputs -o $output.dir -l ${MAPPED}/mappedReads.txt

      ""","CountReadsSummary"
   }
}
