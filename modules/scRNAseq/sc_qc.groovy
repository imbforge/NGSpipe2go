sc_qc = {
    doc title: "sc_qc",
        desc:  "Quality control for single cell multiome experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_qc_vars.outdir
    
    def sc_qc_FLAGS =
        (sc_qc_vars.outdir             ? " outdir="             + sc_qc_vars.outdir             : "") +
        (sc_qc_vars.project            ? " project="            + sc_qc_vars.project            : "") +
        (sc_qc_vars.res                ? " res="                + sc_qc_vars.res                : "") +
        (sc_qc_vars.db                 ? " db="                 + sc_qc_vars.db                 : "") +
        (sc_qc_vars.extra              ? " "                    + sc_qc_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_qc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_qc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_qc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_qc/sc_qc_multiome.R $sc_qc_FLAGS
        ""","sc_qc"
    }
}

