PlotReadLengths = {
   doc title: "PlotReadLengths",
   desc: "Summarizes the number of non-structural reads for each the libraries. It uses the output from featureCounts",
   constraints: "needs a table with gene names and biotypes",
   author: "António Domingues"

   output.dir = CUTADAPT_OUTDIR

   produce("figure/PercentageReadsLengthDistribution.pdf"){
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         cd $output.dir

         Rscript $PLOT_NONSTRUCT_TOOL_PATH 
      ""","PlotReadLengths"
   }
}
