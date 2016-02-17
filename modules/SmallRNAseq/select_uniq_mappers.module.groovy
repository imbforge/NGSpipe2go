//rule for task SelectUniqMappers from catalog NGS, version 1
//desc: Select uniquely mapped reads
SelectUniqMappers = {
   doc title: "SelectUniqMappers",
      desc:  "Call samtools to create a new BAM with only the uniquely mapped reads.",

      author: "Sergi Sayols"

   output.dir = UNIQUEMAP_OUT_DIR

   transform(".bam") to(".unique.bam") {
      exec """

         echo 'VERSION INFO'  1>&2 &&
         ${TOOL_SAMTOOLS} --version 1>&2 &&
         echo '/VERSION INFO' 1>&2 &&

         ${TOOL_SAMTOOLS} view -hb -q 255 $input | ${TOOL_SAMTOOLS} sort -@ $BOWTIE_THREADS - -o $output &&
         ${TOOL_SAMTOOLS} index $output
      ""","BAMindexer"
   }
}
