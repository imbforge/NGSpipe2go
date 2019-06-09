//rule for task feature_count from htseq-count package, version 1
//desc: Counting reads in features with htseq-count
HTseqCount = {
   doc title: "subread_count_se",
   desc:  "Counting reads in features with htseq-count",
   constraints: """Default: unstranded counting, intersection-nonempty. Assumes BAM files""",
   bpipe_version: "tested with bpipe 0.9.9 beta_1",
   author: "Antonio Domingues"

   def EXP = HTSEQCOUNT_GENESGTF.split("/")[-1].replaceAll(".gtf", "")
   HTSEQCOUNT_OUT = HTSEQCOUNT_OUTDIR + "/" + EXP
   output.dir  = HTSEQCOUNT_OUT 
   
   def HTSEQCOUNT_FLAGS = HTSEQCOUNT_FILE    + " " +
                          HTSEQCOUNT_MODE + " " +
                          HTSEQCOUNT_EXTRA    + " "


   // no|yes|reverse
   if(HTSEQCOUNT_STRANDED == "no") {
      HTSEQCOUNT_FLAGS = "-s no " + HTSEQCOUNT_FLAGS
   }
   else if (HTSEQCOUNT_STRANDED == "yes") {
      HTSEQCOUNT_FLAGS = "-s yes " + HTSEQCOUNT_FLAGS
   }
   else {
      HTSEQCOUNT_FLAGS = "-s reverse " + HTSEQCOUNT_FLAGS
   }


   // run the chunk
   transform(".bam") to (".counts") {
      exec """

         htseq-count $HTSEQCOUNT_FLAGS $input $HTSEQCOUNT_GENESGTF > $output 2> ${output.prefix}_htseqcountlog.stderr

      ""","HTseqCount"
   }
}
