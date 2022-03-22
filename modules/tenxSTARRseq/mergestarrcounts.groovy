MergeSTARRCounts = {
    doc title: "MergeSTARRCounts",
        desc:  "Merge STARR mRNA counts with endogenous mRNA counts from 10X snRNA-seq run.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Martin Oti"

    def SCEFILE = new File(MergeSTARRCounts_vars.outputsce)
    output.dir = SCEFILE.getParentFile()

    def MERGESTARRCOUNTS_FLAGS =
        (MergeSTARRCounts_vars.starrgtf   ? " starrgtf="  + MergeSTARRCounts_vars.starrgtf   : "" ) +
        (MergeSTARRCounts_vars.umicounts  ? " umicounts=" + MergeSTARRCounts_vars.umicounts  : "" ) +
        (MergeSTARRCounts_vars.aggr       ? " aggr="      + MergeSTARRCounts_vars.aggr       : "" ) +
        (MergeSTARRCounts_vars.sce        ? " sce="       + MergeSTARRCounts_vars.sce        : "" )

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("MergeSTARRCounts")

    // run the chunk
    produce(MergeSTARRCounts_vars.outputsce) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [[ ! -e "$output.dir" ]]; then
                mkdir -p "$output.dir";
            fi &&

            Rscript ${PIPELINE_ROOT}/tools/MergeSTARRCounts/MergeSTARRCounts.R $MERGESTARRCOUNTS_FLAGS
        ""","MergeSTARRCounts"
    }
    forward inputs
}
