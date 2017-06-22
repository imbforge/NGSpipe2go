//rule for task Bam2FastQ from catalog NGS, version 1
//desc: Convert Bam to FastQ using bedtools
Bam2FastQ = {
   doc title: "Bam2FastQ",
      desc:  "Call bedtools to create a new BAM with only the reads mapped to the genome.",

      author: "Antonio Domingues"

   output.dir = FQ_OUT_DIR

   transform(".bam") to(".fq.gz") {
      exec """
         module load bedtools/${BEDTOOLS_VERSION} &&

         bamToFastq -i $input | gzip > $output

      ""","Bam2FastQ"
   }
}
