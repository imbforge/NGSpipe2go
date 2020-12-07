bamCoverage = {
    doc title: "bamCoverage",
        desc:  "bamCoverage wrapper",
        constraints: "normalised bigwig track for RNA/ChipSeq PE data",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    var subdir : ""
    output.dir = bamCoverage_vars.outdir + "/$subdir"

    def BAMCOVERAGE_FLAGS =
        (bamCoverage_vars.cores     ? " --numberOfProcessors " + bamCoverage_vars.cores : "") +
        (bamCoverage_vars.fragments ? " --extendReads " + (bamCoverage_vars.paired ? "" : bamCoverage_vars.fraglength + " ") : "") +
        (bamCoverage_vars.extra     ? " "                      + bamCoverage_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("deeptools", tools["deeptools"]["version"], tools["deeptools"]["runenv"])
    def PREAMBLE = get_preamble("bamCoverage")

    transform(".bam") to(".bw") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            bamCoverage $BAMCOVERAGE_FLAGS --bam $input -o ${output};
        ""","bamCoverage"
    }
}

