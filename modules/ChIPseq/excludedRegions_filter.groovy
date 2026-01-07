excludedRegions_filter = {
    doc title: "excludedRegions_filter",
        desc: "Remove peaks overlapping genomic regions specified in excludedRegions",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino, modified by Frank Ruehle"

    var subdir : ""
    var excludedRegionsBED: excludedRegions_filter_vars.excludedRegionsBED
    output.dir = excludedRegions_filter_vars.outdir + "/$subdir" 
    
    println "excludedRegions: " + excludedRegionsBED

    def EXCLUDEDREGIONS_FILTER_FLAGS =
        (excludedRegions_filter_vars.files     ? " peakData="           + excludedRegions_filter_vars.files   + "/$subdir" : "") +
        (excludedRegionsBED                    ? " excludedRegionsBED=" + excludedRegionsBED                               : "") +
        (excludedRegions_filter_vars.outdir    ? " out="                + excludedRegions_filter_vars.outdir  + "/$subdir" : "") +
        (excludedRegions_filter_vars.extra     ?                          excludedRegions_filter_vars.extra                : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("excludedRegions_filter.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/excludedRegions_filter/excludedRegions_filter.R $EXCLUDEDREGIONS_FILTER_FLAGS;
        ""","excludedRegions_filter"
    }
}
