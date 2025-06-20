sc_bioc_qc = {
    doc title: "sc_bioc_qc",
        desc:  "Quality control for single cell RNA-Seq experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_bioc_qc_vars.outdir
    
    def sc_bioc_qc_FLAGS =
        (sc_bioc_qc_vars.outdir             ? " outdir="             + sc_bioc_qc_vars.outdir             : "") +
        (sc_bioc_qc_vars.seqtype            ? " seqtype="            + sc_bioc_qc_vars.seqtype            : "") +
        (sc_bioc_qc_vars.pipeline_root      ? " pipeline_root="      + sc_bioc_qc_vars.pipeline_root      : "") +
        (sc_bioc_qc_vars.res                ? " res="                + sc_bioc_qc_vars.res                : "") +
        (sc_bioc_qc_vars.mito_genes         ? " mito_genes="         + sc_bioc_qc_vars.mito_genes         : "") +
        (sc_bioc_qc_vars.annocat_plot       ? " annocat_plot="       + sc_bioc_qc_vars.annocat_plot       : "") +
        (sc_bioc_qc_vars.annocat_plot2      ? " annocat_plot2="      + sc_bioc_qc_vars.annocat_plot2      : "") +
        (sc_bioc_norm_vars.plot_pointsize   ? " plot_pointsize="     + sc_bioc_norm_vars.plot_pointsize   : "") +
        (sc_bioc_norm_vars.plot_pointalpha  ? " plot_pointalpha="    + sc_bioc_norm_vars.plot_pointalpha  : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_bioc_qc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_bioc_qc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_bioc_qc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_qc/sc_bioc_qc.R $sc_bioc_qc_FLAGS
        ""","sc_bioc_qc"
    }
}

