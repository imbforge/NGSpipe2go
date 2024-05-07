sc_integrateRNA = {
    doc title: "sc_integrateRNA",
        desc:  "Batch correct gene expression data in Seurat object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sivarajan Karunanithi"

    output.dir = sc_integrateRNA_vars.outdir
    def sc_integrateRNA_FLAGS =
        (sc_integrateRNA_vars.outdir             ? " outdir="             + sc_integrateRNA_vars.outdir             : "") +
        (sc_integrateRNA_vars.project            ? " project="            + sc_integrateRNA_vars.project            : "") +
        (sc_integrateRNA_vars.res                ? " result="             + sc_integrateRNA_vars.res                : "") +
        (sc_integrateRNA_vars.batch              ? " batch="              + sc_integrateRNA_vars.batch              : "") +
        (sc_integrateRNA_vars.n_features         ? " n_features="         + sc_integrateRNA_vars.n_features         : "") +
        (sc_integrateRNA_vars.rdtype             ? " rdtype="             + sc_integrateRNA_vars.rdtype             : "") +
        (sc_integrateRNA_vars.extra              ? " "                    + sc_integrateRNA_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_integrateRNA module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_integrateRNA should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_integrateRNA.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_integrate/sc_integrateRNA.R $sc_integrateRNA_FLAGS
        ""","sc_integrateRNA"
    }
}

