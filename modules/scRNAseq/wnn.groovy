wnn = {
    doc title: "wnn",
        desc:  "Compute a joint neighbor graph that represent both the gene expression and DNA accessibility measurements",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = wnn_vars.outdir
    
    def wnn_FLAGS =
        (wnn_vars.outdir             ? " outdir="             + wnn_vars.outdir             : "") +
        (wnn_vars.project            ? " project="            + wnn_vars.project            : "") +
        (wnn_vars.res                ? " res="                + wnn_vars.res                : "") +
        (wnn_vars.knn                ? " knn="                + wnn_vars.knn                : "") +
        (wnn_vars.knnRange           ? " knnRange="           + wnn_vars.knnRange           : "") +        
        (wnn_vars.clusterAlg         ? " clusterAlg="         + wnn_vars.clusterAlg         : "") +        
        (wnn_vars.clusterRes         ? " clusterRes="         + wnn_vars.clusterRes         : "") +        
        (wnn_vars.skipFirstLSIcomp   ? " skipFirstLSIcomp="   + wnn_vars.skipFirstLSIcomp   : "") +        
        (wnn_vars.extra              ? " "                    + wnn_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The wnn module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if wnn should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("wnn.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_wnn/wnn.R $wnn_FLAGS
        ""","wnn"
    }
}

