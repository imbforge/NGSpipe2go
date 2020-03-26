// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/scRNAseq/cutadapt.vars.groovy"

Cutadapt = {
        doc title: "Cutadapt",
        desc:  "remove adapter from reads",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch, Nastasja Kreim, Frank Rühle"

    output.dir = Cutadapt_vars.outdir

    // create the log folder if it doesn't exists
    def Cutadapt_LOGDIR = new File( Cutadapt_vars.logdir)
    if (!Cutadapt_LOGDIR.exists()) {
        Cutadapt_LOGDIR.mkdirs()
    }
    // create the discarded folder if it doesn't exists
    def Cutadapt_DISCARDED_DIR = new File( Cutadapt_vars.outdir + "/discarded")
    if (!Cutadapt_DISCARDED_DIR.exists()) {
        Cutadapt_DISCARDED_DIR.mkdirs()
    }

    def Cutadapt_FLAGS =
        (Cutadapt_vars.adapter_sequence    ? " --adapter "        + Cutadapt_vars.adapter_sequence    : "") +
        (Cutadapt_vars.polya               ?                        Cutadapt_vars.polya               : "") +
        (Cutadapt_vars.minimum_overlap     ? " --overlap="        + Cutadapt_vars.minimum_overlap     : "") +
        (Cutadapt_vars.minimum_length_keep ? " --minimum-length " + Cutadapt_vars.minimum_length_keep : "") +
        (Cutadapt_vars.errorrate           ? " --error-rate "     + Cutadapt_vars.errorrate           : "") +
        (Cutadapt_vars.extra               ? " "                  + Cutadapt_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble("Cutadapt")

    transform(".fastq.gz") to (".cutadapt.fastq.gz") {
        def SAMPLENAME = input.prefix.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&

            cutadapt $Cutadapt_FLAGS --too-short-output=\${TMP}/\${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz --output=\${TMP}/\${SAMPLENAME_BASE}.cutadapt.fastq.gz $input 2>&1 >> ${Cutadapt_LOGDIR}/\${SAMPLENAME_BASE}.cutadapt.log &&
		
            mv \${TMP}/\${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz ${Cutadapt_DISCARDED_DIR}/\${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz &&
            mv \${TMP}/\${SAMPLENAME_BASE}.cutadapt.fastq.gz $output

        ""","Cutadapt"
    }
}
