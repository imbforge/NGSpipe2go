Cutadapt = {
   doc title: "Cutadapt",
      desc:  "remove adapter from reads",
      constraints: "Only supports compressed FASTQ files",
      author: "Antonio Domingues"

   output.dir = CUTADAPT_OUTDIR

   transform(".fastq.gz") to (".cutadapt.fastq.gz") {
      exec """

          module load cutadapt/${CUTADAPT_VERSION} &&

         cutadapt $ADAPTER_SEQUENCE -O $MINIMUM_OVERLAP -m $MINIMUM_LENGTH_KEEP -M $MAXIMUM_LENGTH_KEEP -o $output $input 2>&1 >> ${FASTQ_QUALITY_FILTER_OUTDIR}/cutadapt.log
      ""","Cutadapt"
   }
}
