norm_bioc = {
    doc title: "norm_bioc",
        desc:  "Normalize gene expression data in singlecellexperiment object and calculate HVGs and reduced dims",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = norm_bioc_vars.outdir

    def norm_bioc_FLAGS =
        (norm_bioc_vars.outdir             ? " outdir="             + norm_bioc_vars.outdir             : "") +
        (norm_bioc_vars.seqtype            ? " seqtype="            + norm_bioc_vars.seqtype            : "") +
        (norm_bioc_vars.pipeline_root      ? " pipeline_root="      + norm_bioc_vars.pipeline_root      : "") +
        (norm_bioc_vars.res                ? " res="                + norm_bioc_vars.res                : "") +
        (norm_bioc_vars.spikein_norm       ? " spikein_norm="       + norm_bioc_vars.spikein_norm       : "") +
        (norm_bioc_vars.maxgenes2plot      ? " maxgenes2plot="      + norm_bioc_vars.maxgenes2plot      : "") +
        (norm_bioc_vars.org                ? " org="                + norm_bioc_vars.org                : "") +
        (norm_bioc_vars.explanatory_vars   ? " explanatory_vars="   + norm_bioc_vars.explanatory_vars   : "") +
        (norm_bioc_vars.hvg_n              ? " hvg_n="              + norm_bioc_vars.hvg_n              : "") +
        (norm_bioc_vars.block_var          ? " block_var="          + norm_bioc_vars.block_var          : "") +
        (norm_bioc_vars.pca_components     ? " pca_components="     + norm_bioc_vars.pca_components     : "") +
        (norm_bioc_vars.perplexity         ? " perplexity="         + norm_bioc_vars.perplexity         : "") +
        (norm_bioc_vars.n_neighbors        ? " n_neighbors="        + norm_bioc_vars.n_neighbors        : "") +
        (norm_bioc_vars.annocat_plot       ? " annocat_plot="       + norm_bioc_vars.annocat_plot       : "") +
        (norm_bioc_vars.annocat_plot2      ? " annocat_plot2="      + norm_bioc_vars.annocat_plot2      : "") +
        (norm_bioc_vars.plot_pointsize     ? " plot_pointsize="     + norm_bioc_vars.plot_pointsize     : "") +
        (norm_bioc_vars.plot_pointalpha    ? " plot_pointalpha="    + norm_bioc_vars.plot_pointalpha    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The norm_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if norm_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("norm_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_norm/norm_bioc.R $norm_bioc_FLAGS
        ""","norm_bioc"
    }
}

