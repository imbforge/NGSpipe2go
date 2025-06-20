sc_bioc_filter = {
    doc title: "sc_bioc_filter",
        desc:  "Apply quality control filter for single cell RNA-Seq experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_bioc_filter_vars.outdir
    
    def sc_bioc_filter_FLAGS =
        (sc_bioc_filter_vars.outdir                     ? " outdir="                     + sc_bioc_filter_vars.outdir                     : "") +
        (sc_bioc_filter_vars.project                    ? " project="                    + sc_bioc_filter_vars.project                    : "") +
        (sc_bioc_filter_vars.pipeline_root              ? " pipeline_root="              + sc_bioc_filter_vars.pipeline_root              : "") +
        (sc_bioc_filter_vars.res                        ? " res="                        + sc_bioc_filter_vars.res                        : "") +
        (sc_bioc_filter_vars.samples2exclude            ? " samples2exclude="            + sc_bioc_filter_vars.samples2exclude            : "") +
        (sc_bioc_filter_vars.type_of_threshold          ? " type_of_threshold="          + sc_bioc_filter_vars.type_of_threshold          : "") +
        (sc_bioc_filter_vars.threshold_total_counts_min ? " threshold_total_counts_min=" + sc_bioc_filter_vars.threshold_total_counts_min : "") +
        (sc_bioc_filter_vars.threshold_total_counts_max ? " threshold_total_counts_max=" + sc_bioc_filter_vars.threshold_total_counts_max : "") +
        (sc_bioc_filter_vars.threshold_total_detected   ? " threshold_total_detected="   + sc_bioc_filter_vars.threshold_total_detected   : "") +
        (sc_bioc_filter_vars.threshold_pct_counts_Mt    ? " threshold_pct_counts_Mt="    + sc_bioc_filter_vars.threshold_pct_counts_Mt    : "") +
        (sc_bioc_filter_vars.NMADS                      ? " NMADS="                      + sc_bioc_filter_vars.NMADS                      : "") +
        (sc_bioc_filter_vars.category_NMADS             ? " category_NMADS="             + sc_bioc_filter_vars.category_NMADS             : "") +
        (sc_bioc_filter_vars.threshold_low_abundance    ? " threshold_low_abundance="    + sc_bioc_filter_vars.threshold_low_abundance    : "") +
        (sc_bioc_filter_vars.annocat_plot               ? " annocat_plot="               + sc_bioc_filter_vars.annocat_plot               : "") +
        (sc_bioc_norm_vars.plot_pointsize               ? " plot_pointsize="             + sc_bioc_norm_vars.plot_pointsize               : "") +
        (sc_bioc_norm_vars.plot_pointalpha              ? " plot_pointalpha="            + sc_bioc_norm_vars.plot_pointalpha              : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_bioc_filter module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_bioc_filter should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_bioc_filter.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_qc/sc_bioc_filter.R $sc_bioc_filter_FLAGS
        ""","sc_bioc_filter"
    }
}

