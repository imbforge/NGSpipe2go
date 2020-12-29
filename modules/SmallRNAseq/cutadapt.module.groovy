Cutadapt = {
   doc title: "Cutadapt",
      desc:  "remove adapter from reads",
      constraints: "Only supports compressed FASTQ files",
      author: "Antonio Domingues"

   output.dir = CUTADAPT_OUTDIR
   def EXP = input.split("/")[-1].replaceAll(".fastq.gz", "")

    // create the log folder if it doesn't exists
    def CUTADAPT_LOGDIR = new File( CUTADAPT_OUTDIR + "/logs")
    if (!CUTADAPT_LOGDIR.exists()) {
        CUTADAPT_LOGDIR.mkdirs()
    }

   transform(".fastq.gz") to (".cutadapt.fastq.gz", ".cutadapt_discarded.fastq.gz") {
      exec """

         module load cutadapt/${CUTADAPT_VERSION} &&

         cutadapt $ADAPTER_SEQUENCE -O $MINIMUM_OVERLAP -m $MINIMUM_LENGTH_KEEP -M $MAXIMUM_LENGTH_KEEP -o $output1 --too-long-output $output2 $input 2>&1 >> ${CUTADAPT_OUTDIR}/${EXP}.cutadapt.log
         
      ""","Cutadapt"
   }
}
