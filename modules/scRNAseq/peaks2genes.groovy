peaks2genes = {
    doc title: "peaks2genes",
        desc:  "Computing the correlation between gene expression and accessibility at nearby peaks",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = peaks2genes_vars.outdir
    
    def peaks2genes_FLAGS =
        (peaks2genes_vars.outdir             ? " outdir="             + peaks2genes_vars.outdir             : "") +
        (peaks2genes_vars.project            ? " project="            + peaks2genes_vars.project            : "") +
        (peaks2genes_vars.res                ? " res="                + peaks2genes_vars.res                : "") +
        (peaks2genes_vars.db                 ? " db="                 + peaks2genes_vars.db                 : "") +
        (peaks2genes_vars.genes2use          ? " genes2use="          + peaks2genes_vars.genes2use          : "") +
        (peaks2genes_vars.genes2plot         ? " genes2plot="         + peaks2genes_vars.genes2plot         : "") +
        (peaks2genes_vars.plotUpstream       ? " plotUpstream="       + peaks2genes_vars.plotUpstream       : "") +
        (peaks2genes_vars.plotDownstream     ? " plotDownstream="     + peaks2genes_vars.plotDownstream     : "") +
        (peaks2genes_vars.extra              ? " "                    + peaks2genes_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The peaks2genes module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if peaks2genes should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("peaks2genes.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_peaks2genes/peaks2genes.R $peaks2genes_FLAGS
        ""","peaks2genes"
    }
}

