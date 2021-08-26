upsetPlot = {
    doc title: "upset plot",
        desc: "prepare combination matrix and UpSet Plot for peak data",
        constraints: "calculation of combination matrix may take hours for larger projects",
        bpipe_version:"tested with bpipe 0.9.9.8",
        author:"Frank Rühle"

    var subdir : ""
    output.dir = UPSET_vars.outdir + "/$subdir"

    def UPSET_FLAGS =
        (UPSET_vars.files      ? " peakData="       + UPSET_vars.files + "/$subdir"       : "") +
        (UPSET_vars.targets    ? " targets="        + UPSET_vars.targets                  : "") +
        (UPSET_vars.outdir     ? " out="            + UPSET_vars.outdir + "/$subdir"      : "") +
        (UPSET_vars.mode       ? " mode="           + UPSET_vars.mode                     : "") +
        (UPSET_vars.peakOverlapMode  ? " peakOverlapMode=" + UPSET_vars.peakOverlapMode   : "") +
        (UPSET_vars.setsize    ? " setsize="        + UPSET_vars.setsize                  : "") +
        (UPSET_vars.addBarAnnotation ? " addBarAnnotation=" + UPSET_vars.addBarAnnotation : "") +
        (UPSET_vars.extra      ? " "                + UPSET_vars.extra                    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("upsetPlot.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/upsetPlot/upsetPlot.R $UPSET_FLAGS;
        ""","upsetPlot"
    }
}

