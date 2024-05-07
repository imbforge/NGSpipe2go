CRmotifCounts = {
    doc title: "CRmotifCounts",
        desc:  "load motif counts and feature linkage data from Cellranger",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = CRmotifCounts_vars.outdir
    
    def CRmotifCounts_FLAGS =

        (CRmotifCounts_vars.outdir             ? " outdir="             + CRmotifCounts_vars.outdir             : "") +
        (CRmotifCounts_vars.project            ? " project="            + CRmotifCounts_vars.project            : "") +
        (CRmotifCounts_vars.res                ? " res="                + CRmotifCounts_vars.res                : "") +
        (CRmotifCounts_vars.cellranger_aggr_id ? " cellranger_aggr_id=" + CRmotifCounts_vars.cellranger_aggr_id : "") +
        (CRmotifCounts_vars.extra              ? " "                    + CRmotifCounts_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())


    // Only make sense to have this module when using cellranger_aggr. However, we cannot move this module within the 
    // ESSENTIAL_USE_AGGR_DATA as we need the quality filtered seurat object in this step.

    // The CRmotifCounts module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if CRmotifCounts should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("CRmotifCounts.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_motifs/CRmotifCounts.R $CRmotifCounts_FLAGS
        ""","CRmotifCounts"
    }
}

