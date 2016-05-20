CountReads = {
   doc title: "Intersects reads genomic locations",
      desc:  "The ultimate goal is aggregate how many reads (piRNAs) map to each class of genomic features.",
      constraints: "The bed files should contains an identifier, for instance, LINE or SINE, in col7. Uses bedtools v2.23.0",
      author: "Antonio Domingues"


   output.dir = COUNT_READS_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")

   produce(COUNT_READS_OUTDIR + "/" + SAMPLE_NAME + ".sense.counts",
           COUNT_READS_OUTDIR + "/" + SAMPLE_NAME + ".antisense.counts") {


      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_BEDTOOLS}/bedtools --version) 1>&2 ;
         echo '/VERSION INFO'  1>&2 &&

         ${TOOL_BEDTOOLS}/bamToBed -i $input | ${TOOL_BEDTOOLS}/bedtools intersect -a - -b $FEATURES_PATH -s -wa -wb -bed -f 1.0 -nonamecheck > $output1 &&
         ${TOOL_BEDTOOLS}/bamToBed -i $input | ${TOOL_BEDTOOLS}/bedtools intersect -a - -b $FEATURES_PATH -S -wa -wb -bed -f 1.0 -nonamecheck > $output2

      ""","CountReads"
   }
}
