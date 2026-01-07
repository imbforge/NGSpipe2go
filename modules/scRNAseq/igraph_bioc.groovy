igraph_bioc = {
    doc title: "igraph_bioc",
        desc:  "Apply graph-based clustering on singlecellexperiment object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = igraph_bioc_vars.outdir
    
    def igraph_bioc_FLAGS =
        (igraph_bioc_vars.outdir             ? " outdir="             + igraph_bioc_vars.outdir             : "") +
        (igraph_bioc_vars.seqtype            ? " seqtype="            + igraph_bioc_vars.seqtype            : "") +
        (igraph_bioc_vars.pipeline_root      ? " pipeline_root="      + igraph_bioc_vars.pipeline_root      : "") +
        (igraph_bioc_vars.res                ? " res="                + igraph_bioc_vars.res                : "") +
        (igraph_bioc_vars.data2clust         ? " data2clust="         + igraph_bioc_vars.data2clust         : "") +
        (igraph_bioc_vars.pre_kmeans         ? " pre_kmeans="         + igraph_bioc_vars.pre_kmeans         : "") +
        (igraph_bioc_vars.n_neighbors        ? " n_neighbors="        + igraph_bioc_vars.n_neighbors        : "") +
        (igraph_bioc_vars.weighting_scheme   ? " weighting_scheme="   + igraph_bioc_vars.weighting_scheme   : "") +
        (igraph_bioc_vars.algorithm          ? " algorithm="          + igraph_bioc_vars.algorithm          : "") +
        (igraph_bioc_vars.annocat_plot       ? " annocat_plot="       + igraph_bioc_vars.annocat_plot       : "") +
        (igraph_bioc_vars.annocat_plot2      ? " annocat_plot2="      + igraph_bioc_vars.annocat_plot2      : "") +
        (igraph_bioc_vars.plot_pointsize     ? " plot_pointsize="     + igraph_bioc_vars.plot_pointsize     : "") +
        (igraph_bioc_vars.plot_pointalpha    ? " plot_pointalpha="    + igraph_bioc_vars.plot_pointalpha    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The igraph_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if igraph_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("igraph_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_clust/igraph_bioc.R $igraph_bioc_FLAGS
        ""","igraph_bioc"
    }
}

