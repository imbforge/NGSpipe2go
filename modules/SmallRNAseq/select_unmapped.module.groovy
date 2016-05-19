//rule for task SelectUnMapped from catalog NGS, version 1
SelectUnMapped = {
   doc title: "SelectUnMapped",
      desc:  "Call samtools to create a new BAM with only the reads that failed to map to the genome.",
      author: "Antonio Domingues"

   output.dir = UNIQUEMAP_OUT_DIR

   transform(".bam") to(".unmapped.bam") {
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         ${TOOL_SAMTOOLS} --version 1>&2 &&
         echo '/VERSION INFO' 1>&2 &&

         ${TOOL_SAMTOOLS} view -hb -f 4 $input | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_THREADS - -o $output &&
         ${TOOL_SAMTOOLS} index $output
      ""","SelectUnMapped"
   }
}
