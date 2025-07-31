collectCl_bioc = {
    doc title: "collectCl_bioc",
        desc:  "Collects all clustering tsv files and merges them to the sce object.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = collectCl_bioc_vars.clustdir
    
    def collectCl_bioc_FLAGS =
        (collectCl_bioc_vars.pipeline_root       ? " pipeline_root="       + collectCl_bioc_vars.pipeline_root       : "") +
        (collectCl_bioc_vars.seqtype             ? " seqtype="             + collectCl_bioc_vars.seqtype             : "") +
        (collectCl_bioc_vars.res                 ? " res="                 + collectCl_bioc_vars.res                 : "") +
        (collectCl_bioc_vars.clustdir            ? " clustdir="            + collectCl_bioc_vars.clustdir            : "") +
        (collectCl_bioc_vars.selected_clustering ? " selected_clustering=" + collectCl_bioc_vars.selected_clustering : "") +
        (collectCl_bioc_vars.exprPlotDir         ? " exprPlotDir="         + collectCl_bioc_vars.exprPlotDir         : "") +
        (collectCl_bioc_vars.selected_genes      ? " selected_genes="      + collectCl_bioc_vars.selected_genes      : "") +
        (collectCl_bioc_vars.annocat_plot        ? " annocat_plot="        + collectCl_bioc_vars.annocat_plot        : "") +
        (collectCl_bioc_vars.annocat_plot2       ? " annocat_plot2="       + collectCl_bioc_vars.annocat_plot2       : "") +
        (collectCl_bioc_vars.plot_pointsize      ? " plot_pointsize="      + collectCl_bioc_vars.plot_pointsize      : "") +
        (collectCl_bioc_vars.plot_pointalpha     ? " plot_pointalpha="     + collectCl_bioc_vars.plot_pointalpha     : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The collectCl_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if collectCl_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("collectCl_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_clust/collectCl_bioc.R $collectCl_bioc_FLAGS
        ""","collectCl_bioc"
    }
}

