load MODULE_FOLDER + "SmallRNAseq/count_mapped_reads.vars.groovy"

CountMappedReads = {
   doc title: "CountMappedReads",
      desc:  "Count reads mapped in a bam file. Useful for normalization. See: https://www.biostars.org/p/138116/#138118 for explanation of which reads are counted.",
      author: "Antonio Domingues"

   output.dir = READ_COUNTS_OUT

   def SAMPLE_NAME = input.prefix

   transform(".bam") to(".mappedreads.txt") {
      exec """
         module load samtools/${SAMTOOLS_VERSION} &&

         nreads=\$(samtools view -F 0x904 -c $input) &&

         echo ${SAMPLE_NAME} ${nreads} > $output
      ""","CountMappedReads"
   }
}
