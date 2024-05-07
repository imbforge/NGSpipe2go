sc_integrateATAC = {
    doc title: "sc_integrateATAC",
        desc:  "Batch correct ATAC data in the multiome Seurat object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sivarajan Karunanithi"

    output.dir = sc_integrateATAC_vars.outdir
    
    def sc_integrateATAC_FLAGS =
        (sc_integrateATAC_vars.outdir             ? " outdir="             + sc_integrateATAC_vars.outdir             : "") +
        (sc_integrateATAC_vars.project            ? " project="            + sc_integrateATAC_vars.project            : "") +
        (sc_integrateATAC_vars.res                ? " res="                + sc_integrateATAC_vars.res                : "") +
        (sc_integrateATAC_vars.featureCutoff      ? " featureCutoff="      + sc_integrateATAC_vars.featureCutoff      : "") +
        (sc_integrateATAC_vars.skipFirstLSIcomp   ? " skipFirstLSIcomp="   + sc_integrateATAC_vars.skipFirstLSIcomp   : "") +
        (sc_integrateATAC_vars.extra              ? " "                    + sc_integrateATAC_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_integrateATAC module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_integrateATAC should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_integrateATAC.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_integrate/sc_integrateATAC.R $sc_integrateATAC_FLAGS
        ""","sc_integrateATAC"
    }
}

