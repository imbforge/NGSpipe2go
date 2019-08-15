load MODULE_FOLDER + "SmallRNAseq/fastq_quality_filter.vars.groovy"

FastQQualityFilter = {
    doc title: "Remove sequences",
        desc:  "filter reads containing low-quality (Phred score below 20) bases in order to facilitate the PCR duplicates removal.",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues"

    output.dir = FASTQ_QUALITY_FILTER_OUTDIR

    // create the log folder if it doesn't exists
    def FASTQ_QUALITY_FILTER_LOGDIR = new File( FASTQ_QUALITY_FILTER_OUTDIR + "/logs")
    if (!FASTQ_QUALITY_FILTER_LOGDIR.exists()) {
            FASTQ_QUALITY_FILTER_LOGDIR.mkdirs()
    }

    def EXP = input.split("/")[-1].replaceAll(".cutadapt.fastq.gz", "")

    def FASTQ_QUALITY_FILTER_FLAGS = " -q "   + MIN_QUAL  +
                                     " -p "   + MIN_PERCENT  +
                                     " -Q " + QUAL_FORMAT +
                                     FASTQ_QUALITY_FILTER_OTHER

    transform(".fastq.gz") to (".highQ.fastq.gz") {


        exec """
            module load fastx_toolkit/${FASTX_VERSION} &&

            zcat $input | fastq_quality_filter $FASTQ_QUALITY_FILTER_FLAGS -o $output 2>&1 >> ${FASTQ_QUALITY_FILTER_LOGDIR}/${EXP}.fastq_quality_filter.log
      ""","FastQQualityFilter"
   }
}
