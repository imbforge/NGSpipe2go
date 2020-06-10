Cutadapt = {
        doc title: "Cutadapt",
        desc:  "remove adapter from reads for single end and paired end design",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch, Nastasja Kreim, Frank Rühle"

    output.dir = Cutadapt_vars.outdir

    def CUTADAPT_PAIRED = Cutadapt_vars.paired

    def File SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    if (Cutadapt_vars.paired == true) { // optional read2 in paired end design
       def File SAMPLENAME2 = new File(input2.optional.prefix.prefix)
       def SAMPLENAME_BASE2 = SAMPLENAME2.getName()
    }


    // create the log folder if it doesn't exists
    def CUTADAPT_LOGDIR = new File( Cutadapt_vars.logdir)
    if (!CUTADAPT_LOGDIR.exists()) {
        CUTADAPT_LOGDIR.mkdirs()
    }
    // create the discarded folder if it doesn't exists
    def CUTADAPT_DISCARDED_DIR = new File( Cutadapt_vars.outdir + "/discarded")
    if (!CUTADAPT_DISCARDED_DIR.exists()) {
        CUTADAPT_DISCARDED_DIR.mkdirs()
    }

    def CUTADAPT_FLAGS =
        (Cutadapt_vars.adapter_sequence    ? " --adapter "        + Cutadapt_vars.adapter_sequence    : "") +
        (Cutadapt_vars.polya               ?                        Cutadapt_vars.polya               : "") +
        (Cutadapt_vars.minimum_overlap     ? " --overlap="        + Cutadapt_vars.minimum_overlap     : "") +
        (Cutadapt_vars.minimum_length_keep ? " --minimum-length " + Cutadapt_vars.minimum_length_keep : "") +
        (Cutadapt_vars.errorrate           ? " --error-rate "     + Cutadapt_vars.errorrate           : "") +
        (Cutadapt_vars.extra               ? " "                  + Cutadapt_vars.extra               : "")

    def CUTADAPT_FLAGS_PAIRED = 
        (Cutadapt_vars.paired ?  " --too-short-paired-output \${TMP}/${SAMPLENAME_BASE2}.cutadapt_discarded.fastq.gz" +
                                 " -p \${TMP}/${SAMPLENAME_BASE2}.cutadapt.fastq.gz" : "")

    def CUTADAPT_INPUT_SECOND = (Cutadapt_vars.paired ?  " " + input2.optional : "")

    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble("Cutadapt")

    transform("*.fastq.gz") to (".cutadapt.fastq.gz") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            cutadapt $CUTADAPT_FLAGS $CUTADAPT_FLAGS_PAIRED --too-short-output=\${TMP}/${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz --output=\${TMP}/${SAMPLENAME_BASE}.cutadapt.fastq.gz $input1 $CUTADAPT_INPUT_SECOND 2>&1 >> ${CUTADAPT_LOGDIR}/${SAMPLENAME_BASE}.cutadapt.log &&
		
            mv -t ${CUTADAPT_DISCARDED_DIR}/ \${TMP}/*.cutadapt_discarded.fastq.gz &&
            mv \${TMP}/${SAMPLENAME_BASE}.cutadapt.fastq.gz $output &&
            if [ $CUTADAPT_PAIRED = true ]; then
                mv \${TMP}/${SAMPLENAME_BASE2}.cutadapt.fastq.gz $output2;
            fi 

        ""","Cutadapt"
    }
}


