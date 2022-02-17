Cutadapt = {
        doc title: "Cutadapt",
        desc:  "remove adapter from reads for single end and paired end design",
        constraints: "Only supports compressed FASTQ files. Required naming scheme is .fastq.gz for se and .R[12].fastq.gz for pe design.",
        author: "Antonio Domingues, Anke Busch, Nastasja Kreim, Frank Rühle"

    output.dir = Cutadapt_vars.outdir

    def CUTADAPT_PAIRED = Cutadapt_vars.paired

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_PRUNED = SAMPLENAME_BASE.replace(".R1", "") // delete .R1 in combined log file of pe design
    def SAMPLENAME_BASE_R2 = SAMPLENAME_BASE.replace(".R1", ".R2") // needed for paired end only

    // create the log folder if it doesn't exists
    def CUTADAPT_STATSDIR = new File( Cutadapt_vars.statsdir)
    if (!CUTADAPT_STATSDIR.exists()) {
        CUTADAPT_STATSDIR.mkdirs()
    }
    // create the discarded folder if it doesn't exists
    def CUTADAPT_DISCARDED_DIR = new File( Cutadapt_vars.outdir + "/discarded")
    if (!CUTADAPT_DISCARDED_DIR.exists()) {
        CUTADAPT_DISCARDED_DIR.mkdirs()
    }

    def CUTADAPT_FLAGS =
        (Cutadapt_vars.adapter_sequence    ? " --adapter "        + Cutadapt_vars.adapter_sequence    : "") +
        (Cutadapt_vars.minimum_overlap     ? " --overlap="        + Cutadapt_vars.minimum_overlap     : "") +
        (Cutadapt_vars.minimum_length_keep ? " --minimum-length " + Cutadapt_vars.minimum_length_keep : "") +
        (Cutadapt_vars.maximum_length_keep ? " --maximum-length " + Cutadapt_vars.maximum_length_keep : "") +
        (Cutadapt_vars.qualitycutoff       ? (Cutadapt_vars.nextseqtrim ? " --nextseq-trim=" : " --quality-cutoff=") + Cutadapt_vars.qualitycutoff : "") +
        (Cutadapt_vars.errorrate           ? " --error-rate "     + Cutadapt_vars.errorrate           : "") +
        (Cutadapt_vars.extra               ? " "                  + Cutadapt_vars.extra               : "")

    def CUTADAPT_FLAGS_PAIRED = 
        (Cutadapt_vars.paired ?  " --too-short-paired-output \${TMP}/${SAMPLENAME_BASE_R2}.cutadapt_discarded.fastq.gz" +
                                 " -p \${TMP}/${SAMPLENAME_BASE_R2}.cutadapt.fastq.gz" : "")

    def CUTADAPT_FLAGS_TOOLONG =
        (Cutadapt_vars.maximum_length_keep ? " --too-long-output \${TMP}/${SAMPLENAME_BASE}.cutadapt_discardedTooLong.fastq.gz" : "")

    def CUTADAPT_INPUT_SECOND = (Cutadapt_vars.paired ?  " " + input2.optional : "")

    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform("*.fastq.gz") to (".cutadapt.fastq.gz") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            cutadapt $CUTADAPT_FLAGS $CUTADAPT_FLAGS_PAIRED $CUTADAPT_FLAGS_TOOLONG --too-short-output=\${TMP}/${SAMPLENAME_BASE}.cutadapt_discarded.fastq.gz --output=\${TMP}/${SAMPLENAME_BASE}.cutadapt.fastq.gz $input1 $CUTADAPT_INPUT_SECOND 1> ${CUTADAPT_STATSDIR}/${SAMPLENAME_BASE_PRUNED}.cutadapt.log &&
		
            mv -t ${CUTADAPT_DISCARDED_DIR}/ \${TMP}/${SAMPLENAME_BASE_PRUNED}*cutadapt_discarded*.fastq.gz &&
            mv -t $output.dir \${TMP}/${SAMPLENAME_BASE_PRUNED}*cutadapt.fastq.gz
        ""","Cutadapt"
    }
}


