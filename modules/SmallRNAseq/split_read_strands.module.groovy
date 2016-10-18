SplitReadStrands = {
   doc title: "Splits plus and minus strands",
        desc: "Creates two bam files, for reads in the plus and minus strands. Basis to create big wig files for these strands",
        constraints: "",
        author: "Antonio Domingues"

   output.dir = MAPPED_STRANDS

   transform(".bam") to (".sense.bam", ".antisense.bam"){
      exec """
         ${TOOL_SAMTOOLS}/samtools view -hbF 16 $input | ${TOOL_SAMTOOLS}/samtools sort -@ $BOWTIE_THREADS - -o $output1 &&
         ${TOOL_SAMTOOLS}/samtools index $output1 &&

         ${TOOL_SAMTOOLS}/samtools view -hbf 16 $input | ${TOOL_SAMTOOLS}/samtools sort -@ $BOWTIE_THREADS - -o $output2 &&
         ${TOOL_SAMTOOLS}/samtools index $output2

      """, "SplitReadStrands"
   }
}
