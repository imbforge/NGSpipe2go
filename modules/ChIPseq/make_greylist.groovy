make_greylist = {
    doc title: "make greylist",
        desc: "Create greylist from control bam files",
        constraints: "",
        bpipe_version:"",
        author:"Frank Ruehle"

    var subdir : ""
    output.dir = make_greylist_vars.outdir + "/$subdir" 

    def MAKE_GREYLIST_FLAGS =
        (make_greylist_vars.targets           ? " targets="           + make_greylist_vars.targets             : "") +
        (make_greylist_vars.outdir            ? " out="               + make_greylist_vars.outdir + "/$subdir" : "") +        
        (make_greylist_vars.bams              ? " bams="              + make_greylist_vars.bams                : "") +
        (make_greylist_vars.peaks             ? " peaks="             + make_greylist_vars.peaks + "/$subdir"  : "") +
        (make_greylist_vars.karyoFile         ? " kar="               + make_greylist_vars.karyoFile           : "") +
        (make_greylist_vars.reps              ? " reps="              + make_greylist_vars.reps                : "") +
        (make_greylist_vars.sampleSize        ? " sampleSize="        + make_greylist_vars.sampleSize          : "") +
        (make_greylist_vars.pThreshold        ? " pThreshold="        + make_greylist_vars.pThreshold          : "") +
        (make_greylist_vars.maxgap            ? " maxgap="            + make_greylist_vars.maxgap              : "") +
        (make_greylist_vars.extra             ? " "                   + make_greylist_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("greylist.rds") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/excludedRegions_filter/make_greylist.R $MAKE_GREYLIST_FLAGS;
        ""","make_greylist"
    }
}
