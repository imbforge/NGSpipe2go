// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

MPSprofiling = {
    doc title: "Multiplexed protein stability (MPS) profiling",
        desc:  "Processing barcode count data and calculation of protein stability indices (PSIs)",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = MPSprofiling_vars.outdir
    def MPSprofiling_FLAGS =
        (MPSprofiling_vars.targets   ? " targets="   + MPSprofiling_vars.targets   : "") +
        (MPSprofiling_vars.prefix    ? " prefix="    + MPSprofiling_vars.prefix    : "") +
        (MPSprofiling_vars.suffix    ? " suffix="    + MPSprofiling_vars.suffix    : "") +
        (MPSprofiling_vars.inputdir  ? " inputdir="  + MPSprofiling_vars.inputdir  : "") +
        (MPSprofiling_vars.outdir    ? " out="       + MPSprofiling_vars.outdir    : "") +
        (MPSprofiling_vars.logdir    ? " log="       + MPSprofiling_vars.logdir    : "") +       
        (MPSprofiling_vars.expdesign ? " expdesign=" + MPSprofiling_vars.expdesign : "") +
        (MPSprofiling_vars.threshold_rel_countssum      ? " threshold_rel_countssum="      + MPSprofiling_vars.threshold_rel_countssum      : "") +
        (MPSprofiling_vars.removeLowCountsRaw           ? " removeLowCountsRaw="           + MPSprofiling_vars.removeLowCountsRaw           : "") +
        (MPSprofiling_vars.minCountThreshold            ? " minCountThreshold="            + MPSprofiling_vars.minCountThreshold            : "") +
        (MPSprofiling_vars.excludeSeqsNotInAllFractions ? " excludeSeqsNotInAllFractions=" + MPSprofiling_vars.excludeSeqsNotInAllFractions : "") +
        (MPSprofiling_vars.remove_NA_bynuc_PSI          ? " remove_NA_bynuc_PSI="          + MPSprofiling_vars.remove_NA_bynuc_PSI          : "") +
        (MPSprofiling_vars.remove_NA_byaa_pooledPSI     ? " remove_NA_byaa_pooledPSI="     + MPSprofiling_vars.remove_NA_byaa_pooledPSI     : "") +
        (MPSprofiling_vars.invertedDesign               ? " invertedDesign="               + MPSprofiling_vars.invertedDesign               : "") +
        (MPSprofiling_vars.extra                        ? " "                              + MPSprofiling_vars.extra                        : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    produce("MPSprofiling.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/MPSprofiling/MPSprofiling.R $MPSprofiling_FLAGS
        ""","MPSprofiling"
    }
}

