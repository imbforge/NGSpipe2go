Cutadapt = {
    doc title: "Cutadapt",
        desc:  "remove adapter from reads",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch"

    output.dir = Cutadapt_vars.outdir

    // create the log folder if it doesn't exists
    def Cutadapt_LOGDIR = new File(Cutadapt_vars.logdir)
    if (!Cutadapt_LOGDIR.exists()) {
        Cutadapt_LOGDIR.mkdirs()
    }

    def Cutadapt_FLAGS_1 =
        (Cutadapt_vars.adapter_sequence    ? " -a " + Cutadapt_vars.adapter_sequence    : "") +
        (Cutadapt_vars.minimum_overlap     ? " -O " + Cutadapt_vars.minimum_overlap     : "") +
        (Cutadapt_vars.minimum_length_keep ? " -m " + Cutadapt_vars.minimum_length_keep : "") +
        (Cutadapt_vars.maximum_length_keep ? " -M " + Cutadapt_vars.maximum_length_keep : "")

    def Cutadapt_FLAGS_2 =
        (Cutadapt_vars.adapter_sequence    ? " -a " + Cutadapt_vars.adapter_sequence          : "") +
        (Cutadapt_vars.minimum_overlap     ? " -O " + Cutadapt_vars.minimum_overlap           : "") +
        (Cutadapt_vars.minimum_length_keep ? " -m " + Cutadapt_vars.maximum_length_keep_plus1 : "")

    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".fastq.gz") to (".cutadapt.fastq.gz",".cutadapt_discarded.fastq.gz") {
         def SAMPLENAME = input.prefix.prefix
         exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            cutadapt $Cutadapt_FLAGS_1 -o $output1 $input 2>&1 >> ${Cutadapt_vars.logdir}/\${SAMPLENAME_BASE}.cutadapt.log &&
            cutadapt $Cutadapt_FLAGS_2 -o $output2 $input 2>&1 >> ${Cutadapt_vars.logdir}/\${SAMPLENAME_BASE}.cutadapt_discarded.log
        ""","Cutadapt"
    }
}
