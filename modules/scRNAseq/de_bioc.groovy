de_bioc = {
    doc title: "de_bioc",
        desc:  "Differential expression analysis per cluster and run GO enrichment analysis",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = de_bioc_vars.outdir
    
    def de_bioc_FLAGS =
        (de_bioc_vars.outdir              ? " outdir="              + de_bioc_vars.outdir              : "") +
        (de_bioc_vars.seqtype             ? " seqtype="             + de_bioc_vars.seqtype             : "") +
        (de_bioc_vars.pipeline_root       ? " pipeline_root="       + de_bioc_vars.pipeline_root       : "") +
        (de_bioc_vars.res                 ? " res="                 + de_bioc_vars.res                 : "") +
        (de_bioc_vars.contrasts           ? " contrasts="           + de_bioc_vars.contrasts           : "") +
        (de_bioc_vars.selected_clustering ? " selected_clustering=" + de_bioc_vars.selected_clustering : "") +
        (de_bioc_vars.selected_CTanno     ? " selected_CTanno="     + de_bioc_vars.selected_CTanno     : "") +
        (de_bioc_vars.minClusterSize      ? " minClusterSize="      + de_bioc_vars.minClusterSize      : "") +
        (de_bioc_vars.lfc_threshold_DE    ? " lfc_threshold_DE="    + de_bioc_vars.lfc_threshold_DE    : "") +
        (de_bioc_vars.FDR_threshold_DE    ? " FDR_threshold_DE="    + de_bioc_vars.FDR_threshold_DE    : "") +
        (de_bioc_vars.org                 ? " org="                 + de_bioc_vars.org                 : "") +
        (de_bioc_vars.top_genes_for_GO    ? " top_genes_for_GO="    + de_bioc_vars.top_genes_for_GO    : "") +
        (de_bioc_vars.p_threshold_GO      ? " p_threshold_GO="      + de_bioc_vars.p_threshold_GO      : "") +
        (de_bioc_vars.annocat_plot        ? " annocat_plot="        + de_bioc_vars.annocat_plot        : "") +
        (de_bioc_vars.annocat_plot2       ? " annocat_plot2="       + de_bioc_vars.annocat_plot2       : "") +
        (de_bioc_vars.plot_pointsize      ? " plot_pointsize="      + de_bioc_vars.plot_pointsize      : "") +
        (de_bioc_vars.plot_pointalpha     ? " plot_pointalpha="     + de_bioc_vars.plot_pointalpha     : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The de_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if de_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("de_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_de_bioc/de_bioc.R $de_bioc_FLAGS
        ""","de_bioc"
    }
}

