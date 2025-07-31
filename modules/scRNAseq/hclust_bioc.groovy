hclust_bioc = {
    doc title: "hclust_bioc",
        desc:  "Apply hierarchical clustering on singlecellexperiment object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = hclust_bioc_vars.outdir
    
    def hclust_bioc_FLAGS =
        (hclust_bioc_vars.outdir             ? " outdir="             + hclust_bioc_vars.outdir             : "") +
        (hclust_bioc_vars.seqtype            ? " seqtype="            + hclust_bioc_vars.seqtype            : "") +
        (hclust_bioc_vars.pipeline_root      ? " pipeline_root="      + hclust_bioc_vars.pipeline_root      : "") +
        (hclust_bioc_vars.res                ? " res="                + hclust_bioc_vars.res                : "") +
        (hclust_bioc_vars.data2clust         ? " data2clust="         + hclust_bioc_vars.data2clust         : "") +
        (hclust_bioc_vars.hclust_method      ? " hclust_method="      + hclust_bioc_vars.hclust_method      : "") +
        (hclust_bioc_vars.deepSplit          ? " deepSplit="          + hclust_bioc_vars.deepSplit          : "") +
        (hclust_bioc_vars.minClusterSize     ? " minClusterSize="     + hclust_bioc_vars.minClusterSize     : "") +
        (hclust_bioc_vars.annocat_plot       ? " annocat_plot="       + hclust_bioc_vars.annocat_plot       : "") +
        (hclust_bioc_vars.annocat_plot2      ? " annocat_plot2="      + hclust_bioc_vars.annocat_plot2      : "") +
        (hclust_bioc_vars.plot_pointsize     ? " plot_pointsize="     + hclust_bioc_vars.plot_pointsize     : "") +
        (hclust_bioc_vars.plot_pointalpha    ? " plot_pointalpha="    + hclust_bioc_vars.plot_pointalpha    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The hclust_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if hclust_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("hclust_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_clust/hclust_bioc.R $hclust_bioc_FLAGS
        ""","hclust_bioc"
    }
}

