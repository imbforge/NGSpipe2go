filtCells_bioc = {
    doc title: "filtCells_bioc",
        desc:  "Apply quality control filter for single cell RNA-Seq experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = filtCells_bioc_vars.outdir
    
    def filtCells_bioc_FLAGS =
        (filtCells_bioc_vars.outdir                       ? " outdir="                       + filtCells_bioc_vars.outdir                       : "") +
        (filtCells_bioc_vars.project                      ? " project="                      + filtCells_bioc_vars.project                      : "") +
        (filtCells_bioc_vars.pipeline_root                ? " pipeline_root="                + filtCells_bioc_vars.pipeline_root                : "") +
        (filtCells_bioc_vars.res                          ? " res="                          + filtCells_bioc_vars.res                          : "") +
        (filtCells_bioc_vars.samples2exclude              ? " samples2exclude="              + filtCells_bioc_vars.samples2exclude              : "") +
        (filtCells_bioc_vars.type_of_threshold            ? " type_of_threshold="            + filtCells_bioc_vars.type_of_threshold            : "") +
        (filtCells_bioc_vars.threshold_total_counts_min   ? " threshold_total_counts_min="   + filtCells_bioc_vars.threshold_total_counts_min   : "") +
        (filtCells_bioc_vars.threshold_total_counts_max   ? " threshold_total_counts_max="   + filtCells_bioc_vars.threshold_total_counts_max   : "") +
        (filtCells_bioc_vars.threshold_total_detected     ? " threshold_total_detected="     + filtCells_bioc_vars.threshold_total_detected     : "") +
        (filtCells_bioc_vars.threshold_pct_counts_Mt      ? " threshold_pct_counts_Mt="      + filtCells_bioc_vars.threshold_pct_counts_Mt      : "") +
        (filtCells_bioc_vars.threshold_pct_counts_spikein ? " threshold_pct_counts_spikein=" + filtCells_bioc_vars.threshold_pct_counts_spikein : "") +
        (filtCells_bioc_vars.NMADS                        ? " NMADS="                        + filtCells_bioc_vars.NMADS                        : "") +
        (filtCells_bioc_vars.category_NMADS               ? " category="                     + filtCells_bioc_vars.category_NMADS               : "") +
        (filtCells_bioc_vars.threshold_doubletscore       ? " threshold_doubletscore="       + filtCells_bioc_vars.threshold_doubletscore       : "") +
        (filtCells_bioc_vars.threshold_low_abundance      ? " threshold_low_abundance="      + filtCells_bioc_vars.threshold_low_abundance      : "") +
        (filtCells_bioc_vars.annocat_plot                 ? " annocat_plot="                 + filtCells_bioc_vars.annocat_plot                 : "") +
        (filtCells_bioc_vars.plot_pointsize               ? " plot_pointsize="               + filtCells_bioc_vars.plot_pointsize               : "") +
        (filtCells_bioc_vars.plot_pointalpha              ? " plot_pointalpha="              + filtCells_bioc_vars.plot_pointalpha              : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The filtCells_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if filtCells_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("filtCells_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_qc/filtCells_bioc.R $filtCells_bioc_FLAGS
        ""","filtCells_bioc"
    }
}

