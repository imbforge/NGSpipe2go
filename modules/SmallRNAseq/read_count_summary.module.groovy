CountReadsSummary = {
   doc title: "Count Reads Summary",
      desc:  "Agregates the reac counts and plots results. Mighr be modified in the future.",
      constraints: "The bed files should contains an identifier, for instance, LINE or SINE, in col7. Uses bedtools v2.23.0",
      author: "Antonio Domingues"


   output.dir = COUNT_READS_OUTDIR

   produce(COUNT_READS_OUTDIR + "/figure/PercentageOfFeature.pdf") {


      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_R}/bin/Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
         echo '/VERSION${TOOL_R}/bin/Rscript INFO'  1>&2 &&

         cd $COUNT_READS_OUTDIR &&
         ${TOOL_R}/bin/Rscript $AGREGATE_SCRIPT

      ""","CountReadsSummary"
   }
}
