//rule for task feature_count from Subread package, version 1
//desc: Counting reads in features with featureCounts
SubReadCount = {
   doc title: "subread_count_se",
   desc:  "Counting reads in features with feature-count out of the subread package",
   constraints: """Default: strand specific counting.""",
   bpipe_version: "tested with bpipe 0.9.8.7",
   author: "Oliver Drechsel"

   output.dir  = SUBREAD_OUTDIR
   def SUBREAD_FLAGS = SUBREAD_CORES    + " " +
                  SUBREAD_GENESGTF + " " +
                        SUBREAD_EXTRA    + " "


   if(SUBREAD_PAIRED == "yes") {
      SUBREAD_FLAGS = "-p " + SUBREAD_FLAGS
   }

   // no|yes|reverse
   if(SUBREAD_STRANDED == "no") {
      SUBREAD_FLAGS = "-s 0 " + SUBREAD_FLAGS
   }
   else if (SUBREAD_STRANDED == "yes") {
      SUBREAD_FLAGS = "-s 1 " + SUBREAD_FLAGS
   }
   else {
      SUBREAD_FLAGS = "-s 2 " + SUBREAD_FLAGS
   }

   // run the chunk
   from("*.bam") produce("featureCounts_summary.txt") {
      exec """

         module load subread/${SUBREAD_VERSION} &&

         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         featureCounts $SUBREAD_FLAGS -o $output $inputs 2> ${output.prefix}_subreadlog.stderr

      ""","SubReadCount"
   }
}
