load MODULE_FOLDER + "SmallRNAseq/merge_fastq.vars.groovy"

CatFastQ = {
   doc title: "CatFastQ",
      desc:  "merge single ended fastq files.",
      author: "Antonio Domingues"

   output.dir = FQ_OUT_DIR
   def SAMPLE = new File(input1.prefix.prefix)
   produce(SAMPLE + ".merged.fq.gz") {
      exec """

         zcat $inputs > $output

      ""","CatFastQ"
   }
}
