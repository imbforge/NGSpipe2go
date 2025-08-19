findmarkers_bioc = {
    doc title: "findmarkers_bioc",
        desc:  "Identify marker genes per cluster and run GO enrichment analysis",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = findmarkers_bioc_vars.outdir
    
    def findmarkers_bioc_FLAGS =
        (findmarkers_bioc_vars.outdir              ? " outdir="              + findmarkers_bioc_vars.outdir              : "") +
        (findmarkers_bioc_vars.seqtype             ? " seqtype="             + findmarkers_bioc_vars.seqtype             : "") +
        (findmarkers_bioc_vars.pipeline_root       ? " pipeline_root="       + findmarkers_bioc_vars.pipeline_root       : "") +
        (findmarkers_bioc_vars.res                 ? " res="                 + findmarkers_bioc_vars.res                 : "") +
        (findmarkers_bioc_vars.selected_clustering ? " selected_clustering=" + findmarkers_bioc_vars.selected_clustering : "") +
        (findmarkers_bioc_vars.rank_effectsize     ? " rank_effectsize="     + findmarkers_bioc_vars.rank_effectsize     : "") +
        (findmarkers_bioc_vars.block_var           ? " block_var="           + findmarkers_bioc_vars.block_var           : "") +
        (findmarkers_bioc_vars.maxgenes2plot       ? " maxgenes2plot="       + findmarkers_bioc_vars.maxgenes2plot       : "") +
        (findmarkers_bioc_vars.maxgenes_hm         ? " maxgenes_hm="         + findmarkers_bioc_vars.maxgenes_hm         : "") +
        (findmarkers_bioc_vars.org                 ? " org="                 + findmarkers_bioc_vars.org                 : "") +
        (findmarkers_bioc_vars.top_genes_for_GO    ? " top_genes_for_GO="    + findmarkers_bioc_vars.top_genes_for_GO    : "") +
        (findmarkers_bioc_vars.p_threshold_GO      ? " p_threshold_GO="      + findmarkers_bioc_vars.p_threshold_GO      : "") +
        (findmarkers_bioc_vars.annocat_plot        ? " annocat_plot="        + findmarkers_bioc_vars.annocat_plot        : "") +
        (findmarkers_bioc_vars.annocat_plot2       ? " annocat_plot2="       + findmarkers_bioc_vars.annocat_plot2       : "") +
        (findmarkers_bioc_vars.plot_pointsize      ? " plot_pointsize="      + findmarkers_bioc_vars.plot_pointsize      : "") +
        (findmarkers_bioc_vars.plot_pointalpha     ? " plot_pointalpha="     + findmarkers_bioc_vars.plot_pointalpha     : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The findmarkers_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if findmarkers_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("findmarkers_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_findmarkers/findmarkers_bioc.R $findmarkers_bioc_FLAGS
        ""","findmarkers_bioc"
    }
}

