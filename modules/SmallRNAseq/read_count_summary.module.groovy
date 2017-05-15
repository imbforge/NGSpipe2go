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
         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_R}/bin/Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
         echo '/VERSION${TOOL_R}/bin/Rscript INFO'  1>&2 &&

         cd $COUNT_READS_OUTDIR &&
         ${TOOL_R}/bin/Rscript $AGREGATE_SCRIPT

      ""","CountReadsSummary"
   }
}
