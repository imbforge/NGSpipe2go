sc_filter = {
    doc title: "sc_filter",
        desc:  "Apply quality control filter for single cell multiome experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_filter_vars.outdir
    
    def sc_filter_FLAGS =

        (sc_filter_vars.outdir             ? " outdir="             + sc_filter_vars.outdir             : "") +
        (sc_filter_vars.project            ? " project="            + sc_filter_vars.project            : "") +
        (sc_filter_vars.res                ? " res="                + sc_filter_vars.res                : "") +
        (sc_filter_vars.nCount_ATAC_min    ? " nCount_ATAC_min="    + sc_filter_vars.nCount_ATAC_min    : "") +
        (sc_filter_vars.nCount_ATAC_max    ? " nCount_ATAC_max="    + sc_filter_vars.nCount_ATAC_max    : "") +
        (sc_filter_vars.nCount_RNA_min     ? " nCount_RNA_min="     + sc_filter_vars.nCount_RNA_min     : "") +
        (sc_filter_vars.nCount_RNA_max     ? " nCount_RNA_max="     + sc_filter_vars.nCount_RNA_max     : "") +
        (sc_filter_vars.FRiPmin            ? " FRiPmin="            + sc_filter_vars.FRiPmin            : "") +
        (sc_filter_vars.FRiBLmax           ? " FRiBLmax="           + sc_filter_vars.FRiBLmax           : "") +
        (sc_filter_vars.nucleosome_sig_max ? " nucleosome_sig_max=" + sc_filter_vars.nucleosome_sig_max : "") +
        (sc_filter_vars.TSS_enrich_min     ? " TSS_enrich_min="     + sc_filter_vars.TSS_enrich_min     : "") +
        (sc_filter_vars.extra              ? " "                    + sc_filter_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_filter module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_filter should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_filter.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_qc/sc_filter_multiome.R $sc_filter_FLAGS
        ""","sc_filter"
    }
}

