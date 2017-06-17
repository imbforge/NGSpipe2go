//rule for task SelectUniqMappers from catalog SmallRNAseq, version 0.1
SelectUniqMappers = {
   doc title: "SelectUniqMappers",
      desc:  "Call samtools to Select uniquely mapped reads.",
      author: "Antonio Domingues"

   output.dir = UNIQUEMAP_OUT_DIR

   transform(".bam") to(".unique.bam") {
      exec """
         echo 'VERSION INFO'  1>&2 &&
         ${TOOL_SAMTOOLS}/samtools --version 1>&2 &&
         echo '/VERSION INFO' 1>&2 &&

         ${TOOL_SAMTOOLS}/samtools view -hb -q 255 $input | ${TOOL_SAMTOOLS}/samtools sort -@ $BOWTIE_THREADS - -o $output &&
         ${TOOL_SAMTOOLS}/samtools index $output
      ""","SelectUniqMappers"
   }
}
