SplitReadStrands = {
   doc title: "Splits plus and minus strands",
        desc: "Creates two bam files, fot reads in the plus and minus strands. Basis to create big wig files for these strands",
        constraints: "",
        author: "Antonio Domingues"

   output.dir = MAPPED_STRANDS

   transform(".bam") to (".sense.bam", ".antisense.bam"){
      exec """
         samtools view -hbF 16 $input | samtools sort -@ $BOWTIE_THREADS - $output1.prefix &&
         samtools index $output1 &&

         samtools view -hbf 16 $input | samtools sort -@ $BOWTIE_THREADS - $output2.prefix &&
         samtools index $output2

      """, "SplitReadStrands"
   }
}
