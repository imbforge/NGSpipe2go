// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/scRNAseq/cutadapt.vars.groovy"

Cutadapt = {
        doc title: "Cutadapt",
        desc:  "remove adapter from reads",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch, Nastasja Kreim, Frank Rühle"

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

    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble("Cutadapt")

    transform(".fastq.gz") to (".cutadapt.fastq.gz") {
        def SAMPLENAME = input.prefix.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&

            cutadapt $CUTADAPT_FLAGS --too-short-output=\${TMP}/\${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz --output=\${TMP}/\${SAMPLENAME_BASE}.cutadapt.fastq.gz $input 2>&1 >> ${CUTADAPT_LOGDIR}/${SAMPLENAME_BASE}.cutadapt.log &&
		
            mv \${TMP}/\${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz ${CUTADAPT_DISCARDED_DIR}/\${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz &&
            mv \${TMP}/\${SAMPLENAME_BASE}.cutadapt.fastq.gz $output

        ""","Cutadapt"
    }
}
