sc_bioc_norm = {
    doc title: "sc_bioc_norm",
        desc:  "Normalize gene expression data in singlecellexperiment object and calculate HVGs and reduced dims",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_bioc_norm_vars.outdir
    
    def sc_bioc_norm_FLAGS =
        (sc_bioc_norm_vars.outdir             ? " outdir="             + sc_bioc_norm_vars.outdir             : "") +
        (sc_bioc_norm_vars.seqtype            ? " seqtype="            + sc_bioc_norm_vars.seqtype            : "") +
        (sc_bioc_norm_vars.pipeline_root      ? " pipeline_root="      + sc_bioc_norm_vars.pipeline_root      : "") +
        (sc_bioc_norm_vars.res                ? " res="                + sc_bioc_norm_vars.res                : "") +
        (sc_bioc_norm_vars.annocat_plot       ? " annocat_plot="       + sc_bioc_norm_vars.annocat_plot       : "") +
        (sc_bioc_norm_vars.annocat_plot2      ? " annocat_plot2="      + sc_bioc_norm_vars.annocat_plot2      : "") +
        (sc_bioc_norm_vars.maxgenes2plot      ? " maxgenes2plot="      + sc_bioc_norm_vars.maxgenes2plot      : "") +
        (sc_bioc_norm_vars.plot_pointsize     ? " plot_pointsize="     + sc_bioc_norm_vars.plot_pointsize     : "") +
        (sc_bioc_norm_vars.plot_pointalpha    ? " plot_pointalpha="    + sc_bioc_norm_vars.plot_pointalpha    : "") +
        (sc_bioc_norm_vars.org                ? " org="                + sc_bioc_norm_vars.org                : "") +
        (sc_bioc_norm_vars.explanatory_vars   ? " explanatory_vars="   + sc_bioc_norm_vars.explanatory_vars   : "") +
        (sc_bioc_norm_vars.hvg_prop           ? " hvg_prop="           + sc_bioc_norm_vars.hvg_prop           : "") +
        (sc_bioc_norm_vars.block_var          ? " block_var="          + sc_bioc_norm_vars.block_var          : "") +
        (sc_bioc_norm_vars.perplexity         ? " perplexity="         + sc_bioc_norm_vars.perplexity         : "") +
        (sc_bioc_norm_vars.n_neighbors        ? " n_neighbors="        + sc_bioc_norm_vars.n_neighbors        : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_bioc_norm module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_bioc_norm should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_bioc_norm.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_norm/sc_bioc_norm.R $sc_bioc_norm_FLAGS
        ""","sc_bioc_norm"
    }
}

