SCTransform = {
    doc title: "SCTransform",
        desc:  "normalize gene expression data in Seurat object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = SCTransform_vars.outdir
    
    def SCTransform_FLAGS =
        (SCTransform_vars.outdir             ? " outdir="             + SCTransform_vars.outdir             : "") +
        (SCTransform_vars.project            ? " project="            + SCTransform_vars.project            : "") +
        (SCTransform_vars.res                ? " res="                + SCTransform_vars.res                : "") +
        (SCTransform_vars.extra              ? " "                    + SCTransform_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The SCTransform module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if SCTransform should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("SCTransform.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_norm/SCTransform.R $SCTransform_FLAGS
        ""","SCTransform"
    }
}

