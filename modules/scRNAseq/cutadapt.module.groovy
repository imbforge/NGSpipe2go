load MODULE_FOLDER + "scRNAseq/cutadapt.vars.groovy"

Cutadapt = {
        doc title: "Cutadapt",
        desc:  "remove adapter from reads",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch, Nastasja Kreim"

    output.dir = CUTADAPT_OUTDIR

    // create the log folder if it doesn't exists
    def CUTADAPT_LOGDIR = new File( LOGS + "/Cutadapt")
    if (!CUTADAPT_LOGDIR.exists()) {
        CUTADAPT_LOGDIR.mkdirs()
    }
    // create the discarded folder if it doesn't exists
    def CUTADAPT_DISCARDED_DIR = new File( CUTADAPT_OUTDIR + "/discarded")
    if (!CUTADAPT_DISCARDED_DIR.exists()) {
        CUTADAPT_DISCARDED_DIR.mkdirs()
    }
    def CUTADAPT_FLAGS = CUTADAPT_ADAPTER_SEQUENCE + 
                         " " + CUTADAPT_POLYA + 
                         " " + CUTADAPT_MINIMUM_OVERLAP +
                         " " + CUTADAPT_MINIMUM_LENGTH_KEEP +
                         " " + CUTADAPT_ERRORRATE +
                         " " + CUTADAPT_EXTRA

    transform(".fastq.gz") to (".cutadapt.fastq.gz") {
        def SAMPLENAME = input.prefix.prefix
        exec """
            module load cutadapt/${CUTADAPT_VERSION} &&
            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            cutadapt $CUTADAPT_FLAGS --too-short-output=${CUTADAPT_DISCARDED_DIR}/${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz --output=$output $input 2>&1 >> ${CUTADAPT_LOGDIR}/${SAMPLENAME_BASE}.cutadapt.log --info-file=${CUTADAPT_LOGDIR}/${SAMPLENAME_BASE}.cutadapt.info
        ""","Cutadapt"
    }
}
