load MODULE_FOLDER + "SmallRNAseq/bam2fastq.vars.groovy"
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
