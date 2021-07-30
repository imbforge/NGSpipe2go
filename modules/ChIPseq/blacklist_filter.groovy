blacklist_filter = {
    doc title: "blacklist_filter",
        desc: "Remove peaks overlapping blacklisted genomic regions",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino, modified by Frank Ruehle"

    var subdir : ""
    output.dir = blacklist_filter_vars.outdir + "/$subdir" 

    def BLACKLIST_FILTER_FLAGS =
        (blacklist_filter_vars.files     ? " peakData="         + blacklist_filter_vars.files   + "/$subdir" : "") +
        (blacklist_filter_vars.blacklist ? " blacklistRegions=" + blacklist_filter_vars.blacklist            : "") +
        (blacklist_filter_vars.outdir    ? " out="              + blacklist_filter_vars.outdir  + "/$subdir" : "") +
        (blacklist_filter_vars.extra     ?                        blacklist_filter_vars.extra                : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("BlackList_Filter.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/BlackList_Filter/BlackList_Filter.R $BLACKLIST_FILTER_FLAGS;
        ""","blacklist_filter"
    }
}
