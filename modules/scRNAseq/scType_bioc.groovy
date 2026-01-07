scType_bioc = {
    doc title: "scType_bioc",
        desc:  "Cell type annotation with scType",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = scType_bioc_vars.outdir
    
    def scType_bioc_FLAGS =
        (scType_bioc_vars.outdir              ? " outdir="              + scType_bioc_vars.outdir              : "") +
        (scType_bioc_vars.seqtype             ? " seqtype="             + scType_bioc_vars.seqtype             : "") +
        (scType_bioc_vars.pipeline_root       ? " pipeline_root="       + scType_bioc_vars.pipeline_root       : "") +
        (scType_bioc_vars.res                 ? " res="                 + scType_bioc_vars.res                 : "") +
        (scType_bioc_vars.selected_clustering ? " selected_clustering=" + scType_bioc_vars.selected_clustering : "") +
        (scType_bioc_vars.dbfile              ? " dbfile="              + scType_bioc_vars.dbfile              : "") +
        (scType_bioc_vars.tissueType          ? " tissueType="          + scType_bioc_vars.tissueType          : "") +
        (scType_bioc_vars.ctcolumn            ? " ctcolumn="            + scType_bioc_vars.ctcolumn            : "") +
        (scType_bioc_vars.annocat_plot        ? " annocat_plot="        + scType_bioc_vars.annocat_plot        : "") +
        (scType_bioc_vars.annocat_plot2       ? " annocat_plot2="       + scType_bioc_vars.annocat_plot2       : "") +
        (scType_bioc_vars.plot_pointsize      ? " plot_pointsize="      + scType_bioc_vars.plot_pointsize      : "") +
        (scType_bioc_vars.plot_pointalpha     ? " plot_pointalpha="     + scType_bioc_vars.plot_pointalpha     : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The scType_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if scType_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("scType_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/CTanno/scType_bioc.R $scType_bioc_FLAGS
        ""","scType_bioc"
    }
}

