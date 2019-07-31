load MODULE_FOLDER + "SmallRNAseq/select_uniq_mappers.vars.groovy"

SelectUniqMappers = {
   doc title: "SelectUniqMappers",
      desc:  "Call samtools to Select uniquely mapped reads.",
      author: "Antonio Domingues"

   output.dir = UNIQUEMAP_OUT_DIR

   transform(".bam") to(".unique.bam") {
      exec """
         module load samtools/${SAMTOOLS_VERSION} &&

         samtools view -hb -q 255 $input | samtools sort -@ $ESSENTIAL_THREADS - -o $output &&
         samtools index $output
      ""","SelectUniqMappers"
   }
}
