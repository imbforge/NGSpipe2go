umidedup = {
    doc title: "deduplication based on UMIs",
        desc: "Deduplication of mapped data using UMIs with umi_tools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Nastasja Kreim"

    output.dir = umidedup_vars.outdir
    def umidedup_FLAGS =
        (umidedup_vars.verbose ? "--verbose=1 " : "") +
        (umidedup_vars.paired  ? "--paired "    : "") +
        (umidedup_vars.param   ? " " + umidedup_vars.param : "") +
        (umidedup_vars.extra   ? " " + umidedup_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to (".umidedup.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            umi_tools dedup $umidedup_FLAGS -I $input -S $output --output-stats=${output.prefix}.stats
        ""","umidedup"
    }
}
