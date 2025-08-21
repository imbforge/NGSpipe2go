collectCT_bioc = {
    doc title: "collectCT_bioc",
        desc:  "Collects all cell type annotation tsv files and merges them to the sce object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = collectCT_bioc_vars.ctanno_dir
    
    def collectCT_bioc_FLAGS =
        (collectCT_bioc_vars.pipeline_root       ? " pipeline_root="       + collectCT_bioc_vars.pipeline_root       : "") +
        (collectCT_bioc_vars.seqtype             ? " seqtype="             + collectCT_bioc_vars.seqtype             : "") +
        (collectCT_bioc_vars.res                 ? " res="                 + collectCT_bioc_vars.res                 : "") +
        (collectCT_bioc_vars.ctanno_dir          ? " ctanno_dir="          + collectCT_bioc_vars.ctanno_dir          : "") +
        (collectCT_bioc_vars.selected_CTanno     ? " selected_CTanno="     + collectCT_bioc_vars.selected_CTanno     : "") +
        (collectCT_bioc_vars.exprPlotDir         ? " exprPlotDir="         + collectCT_bioc_vars.exprPlotDir         : "") +
        (collectCT_bioc_vars.selected_genes      ? " selected_genes="      + collectCT_bioc_vars.selected_genes      : "") +
        (collectCT_bioc_vars.annocat_plot        ? " annocat_plot="        + collectCT_bioc_vars.annocat_plot        : "") +
        (collectCT_bioc_vars.annocat_plot2       ? " annocat_plot2="       + collectCT_bioc_vars.annocat_plot2       : "") +
        (collectCT_bioc_vars.plot_pointsize      ? " plot_pointsize="      + collectCT_bioc_vars.plot_pointsize      : "") +
        (collectCT_bioc_vars.plot_pointalpha     ? " plot_pointalpha="     + collectCT_bioc_vars.plot_pointalpha     : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The collectCT_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if collectCT_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("collectCT_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/CTanno/collectCT_bioc.R $collectCT_bioc_FLAGS
        ""","collectCT_bioc"
    }
}

