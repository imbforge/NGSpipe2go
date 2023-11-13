diffExprSeurat = {
    doc title: "diffExprSeurat",
        desc:  "Differential expression testing with Seurat between all pairs of groups for each cluster and celltype",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = diffExprSeurat_vars.outdir
    
    def diffExprSeurat_FLAGS =
        (diffExprSeurat_vars.outdir             ? " outdir="             + diffExprSeurat_vars.outdir             : "") +
        (diffExprSeurat_vars.project            ? " project="            + diffExprSeurat_vars.project            : "") +
        (diffExprSeurat_vars.res                ? " res="                + diffExprSeurat_vars.res                : "") +
        (diffExprSeurat_vars.assay              ? " assay="              + diffExprSeurat_vars.assay              : "") +
        (diffExprSeurat_vars.minCells           ? " minCells="           + diffExprSeurat_vars.minCells           : "") +
        (diffExprSeurat_vars.clusterVar         ? " clusterVar="         + diffExprSeurat_vars.clusterVar         : "") +
        (diffExprSeurat_vars.CTannoSelected     ? " CTannoSelected="     + diffExprSeurat_vars.CTannoSelected     : "") +
        (diffExprSeurat_vars.test               ? " test="               + diffExprSeurat_vars.test               : "") +
        (diffExprSeurat_vars.latentVars         ? " latentVars="         + diffExprSeurat_vars.latentVars         : "") +
        (diffExprSeurat_vars.batchCorrection    ? " batchCorrection="    + diffExprSeurat_vars.batchCorrection    : "") +
        (diffExprSeurat_vars.extra              ? " "                    + diffExprSeurat_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The diffExprSeurat module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if diffExprSeurat should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("diffExprSeurat.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_diffExprSeurat/diffExprSeurat.R $diffExprSeurat_FLAGS
        ""","diffExprSeurat"
    }
}

