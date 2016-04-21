//rule for task CatFastQ from catalog NGS, version 1
//desc: Merge SE FastQ files
CatFastQ = {
   doc title: "CatFastQ",
      desc:  "merge single ended fastq files.",

      author: "Antonio Domingues"

   output.dir = FQ_OUT_DIR
   def SAMPLE = new File(input1.prefix.prefix)
   produce(SAMPLE + ".merged.fq.gz") {
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         zcat $inputs > $output

      ""","CatFastQ"
   }
}
