//rule for task feature_count from htseq-count package, version 1
//desc: Counting reads in features with htseq-count
HTseqCount = {
   doc title: "subread_count_se",
   desc:  "Counting reads in features with htseq-count",
   constraints: """Default: unstranded counting, intersection-nonempty. Assumes BAM files""",
   bpipe_version: "tested with bpipe 0.9.9 beta_1",
   author: "Antonio Domingues"

   def EXP = GTF.split("/")[-1].replaceAll(".gtf", "")
   output.dir  = HTSEQCOUNT_OUTDIR + "/" + EXP
   
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
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 ;
         echo \$(htseq-count 2>&1 | grep version) 1>&2 ;
         echo '/VERSION INFO' 1>&2 ;

         htseq-count $HTSEQCOUNT_FLAGS $input $GTF > $output 2> ${output.prefix}_htseqcountlog.stderr

      ""","HTseqCount"
   }
}
