diffPeaks = {
    doc title: "diffPeaks",
        desc:  "Identify differentially accessible peaks between all pairs of groups for each cluster and celltype",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = diffPeaks_vars.outdir
    
    def diffPeaks_FLAGS =
        (diffPeaks_vars.outdir             ? " outdir="             + diffPeaks_vars.outdir             : "") +
        (diffPeaks_vars.project            ? " project="            + diffPeaks_vars.project            : "") +
        (diffPeaks_vars.res                ? " res="                + diffPeaks_vars.res                : "") +
        (diffPeaks_vars.assay              ? " assay="              + diffPeaks_vars.assay              : "") +
        (diffPeaks_vars.minCells           ? " minCells="           + diffPeaks_vars.minCells           : "") +
        (diffPeaks_vars.clusterVar         ? " clusterVar="         + diffPeaks_vars.clusterVar         : "") +
        (diffPeaks_vars.CTannoSelected     ? " CTannoSelected="     + diffPeaks_vars.CTannoSelected     : "") +
        (diffPeaks_vars.test               ? " test="               + diffPeaks_vars.test               : "") +
        (diffPeaks_vars.latentVars         ? " latentVars="         + diffPeaks_vars.latentVars         : "") +
        (diffPeaks_vars.extra              ? " "                    + diffPeaks_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The diffPeaks module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if diffPeaks should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("diffPeaks.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_DNAaccess/diffPeaks.R $diffPeaks_FLAGS
        ""","diffPeaks"
    }
}

