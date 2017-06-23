CountNonStrutReads = {
   doc title: "CountNonStrutReads",
   desc: "Summarizes the number of non-structural reads for each the libraries. It uses the output from featureCounts",
   constraints: "needs a table with gene names and biotypes",
   author: "António Domingues"

   output.dir = COUNT_NONSTRUCT_OUTDIR

   produce("normlization_factors.txt"){
      exec """

         module load R/${R_VERSION} &&
         
         cd $output.dir &&

         Rscript $COUNT_NONSTRUCT_TOOL_PATH $input $ESSENTIAL_BIOTYPES_TABLE
         
      ""","CountNonStrutReads"
   }
}
