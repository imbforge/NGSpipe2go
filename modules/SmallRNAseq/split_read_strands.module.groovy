SplitReadStrands = {
   doc title: "Splits plus and minus strands",
        desc: "Creates two bam files, fot reads in the plus and minus strands. Basis to create big wig files for these strands",
        constraints: "",
        author: "Antonio Domingues"

   output.dir = MAPPED_STRANDS

   transform(".bam") to (".sense.bam", ".antisense.bam"){
      exec """
         ${TOOL_SAMTOOLS} view -hbF 16 $input | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_THREADS - -o $output1 &&
         ${TOOL_SAMTOOLS} index $output1 &&

         ${TOOL_SAMTOOLS} view -hbf 16 $input | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_THREADS - -o $output2 &&
         ${TOOL_SAMTOOLS} index $output2

      """, "SplitReadStrands"
   }
}
