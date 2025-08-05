qc_bioc = {
    doc title: "qc_bioc",
        desc:  "Quality control for single cell RNA-Seq experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = qc_bioc_vars.outdir
    
    def qc_bioc_FLAGS =
        (qc_bioc_vars.outdir             ? " outdir="             + qc_bioc_vars.outdir             : "") +
        (qc_bioc_vars.seqtype            ? " seqtype="            + qc_bioc_vars.seqtype            : "") +
        (qc_bioc_vars.pipeline_root      ? " pipeline_root="      + qc_bioc_vars.pipeline_root      : "") +
        (qc_bioc_vars.res                ? " res="                + qc_bioc_vars.res                : "") +
        (qc_bioc_vars.mito_genes         ? " mito_genes="         + qc_bioc_vars.mito_genes         : "") +
        (qc_bioc_vars.spikein_genes      ? " spikein_genes="      + qc_bioc_vars.spikein_genes      : "") +
        (qc_bioc_vars.doubletscore_group ? " doubletscore_group=" + qc_bioc_vars.doubletscore_group : "") +
        (qc_bioc_vars.annocat_plot       ? " annocat_plot="       + qc_bioc_vars.annocat_plot       : "") +
        (qc_bioc_vars.annocat_plot2      ? " annocat_plot2="      + qc_bioc_vars.annocat_plot2      : "") +
        (qc_bioc_vars.plot_pointsize     ? " plot_pointsize="     + qc_bioc_vars.plot_pointsize     : "") +
        (qc_bioc_vars.plot_pointalpha    ? " plot_pointalpha="    + qc_bioc_vars.plot_pointalpha    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The qc_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if qc_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("qc_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_qc/qc_bioc.R $qc_bioc_FLAGS
        ""","qc_bioc"
    }
}

