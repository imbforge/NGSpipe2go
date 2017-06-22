CountReadsSummary = {
   doc title: "Count Reads Summary",
      desc:  "Aggregates the read counts and plots results. Might be modified in the future.",
      constraints: "The bed files should contains an identifier, for instance, LINE or SINE, in col7. Uses bedtools v2.23.0",
      author: "Antonio Domingues"


   output.dir = COUNT_READS_OUTDIR

   produce(
      COUNT_READS_OUTDIR + "/figure/PercentageOfFeature.pdf",
      COUNT_READS_OUTDIR + "/piRNA_quantification.RData"
      ) {


      exec """
         module load R/${R_VERSION} &&

         cd $COUNT_READS_OUTDIR &&
         Rscript $AGREGATE_SCRIPT

      ""","CountReadsSummary"
   }
}
