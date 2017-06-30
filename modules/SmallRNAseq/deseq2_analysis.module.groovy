PerformDEGAnalaysis = {
   doc title: "PerformDEGAnalaysis",
   desc: "Creates DEG report",
   constraints: "needs a table with gene names and biotypes",
   author: "António Domingues"

   output.dir = HTSEQCOUNT_OUT

   produce("normlization_factors.txt"){
      exec """

         module load R/${R_VERSION} &&
         
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         cd $output.dir &&
         ln -s $DEG_SCRIPT . &&
         ln -s $DEG_SCRIPT_HELPERS . &&

         knit_to_html.sh $DEG_SCRIPT 
         
      ""","PerformDEGAnalaysis"
   }
}
