umicount = {
    doc title: "Deduplication and Counting reads per gene",
        desc: "Deduplication and counting of mapped data and splitting accoring to cellbarcode with umi_tools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Nastasja Kreim"

    output.dir = umicount_vars.outdir

    // create the log folder if it doesn't exists
    def umicount_LOGDIR = new File(umicount_vars.logdir)
    if (!umicount_LOGDIR.exists()) {
        umicount_LOGDIR.mkdirs()
    }

    def umicount_FLAGS = 
        (umicount_vars.verbose ? "--verbose=1 " : "") +
        (umicount_vars.paired  ? "--paired "    : "") +
        (umicount_vars.param   ? " " + umicount_vars.param : "") +
        (umicount_vars.extra   ? " " + umicount_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:"", input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to (".umicount.tsv.gz") {
        def SAMPLENAME = input.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            umi_tools count $umicount_FLAGS -I $input -S $output -L ${umicount_LOGDIR}/\${SAMPLENAME_BASE}.umicount.log -E ${umicount_LOGDIR}/\${SAMPLENAME_BASE}.umicount.error 
        ""","umicount"
    }
}
