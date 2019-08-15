load MODULE_FOLDER + "SmallRNAseq/select_unmapped.vars.groovy"

SelectUnMapped = {
   doc title: "SelectUnMapped",
      desc:  "Call samtools to create a new BAM with only the reads that failed to map to the genome.",
      author: "Antonio Domingues"

   output.dir = UNIQUEMAP_OUT_DIR

   transform(".bam") to(".unmapped.bam") {
      exec """
         module load samtools/${SAMTOOLS_VERSION} &&

         samtools view -hb -f 4 $input | samtools sort -@ $BOWTIE_THREADS - -o $output &&
         samtools index $output
      ""","SelectUnMapped"
   }
}
