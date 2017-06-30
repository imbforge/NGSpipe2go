FilterRNAClasses = {
   doc title: "Filter small RNA classes",
         desc:  "Filter alignments to select C. elegans small RNA classes: 21U, 22G, 26G",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = FILTER_CLASSES_OUTDIR

   transform(".bam") to (".21U.bam", ".22G.bam", ".26G.bam"){

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         python $FILTER_CLASSES_TOOL_PATH -i $input -c 21U -o $output1 &&
         python $FILTER_CLASSES_TOOL_PATH -i $input -c 22G -o $output2 &&
         python $FILTER_CLASSES_TOOL_PATH -i $input -c 26G -o $output3

      ""","FilterRNAClasses"
   }
   forward input, output1, output2, output3
}
