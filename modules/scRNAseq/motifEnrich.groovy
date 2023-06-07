motifEnrich = {
    doc title: "motifEnrich",
        desc:  "Identify enriched motifs in sets of differentially accessible peaks between all pairs of groups for each cluster and celltype",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = motifEnrich_vars.outdir
    
    def motifEnrich_FLAGS =
        (motifEnrich_vars.outdir             ? " outdir="             + motifEnrich_vars.outdir             : "") +
        (motifEnrich_vars.project            ? " project="            + motifEnrich_vars.project            : "") +
        (motifEnrich_vars.res                ? " res="                + motifEnrich_vars.res                : "") +
        (motifEnrich_vars.db                 ? " db="                 + motifEnrich_vars.db                 : "") +
        (motifEnrich_vars.diffPeaks_dir      ? " diffPeaks_dir="      + motifEnrich_vars.diffPeaks_dir      : "") +
        (motifEnrich_vars.pval_thresh        ? " pval_thresh="        + motifEnrich_vars.pval_thresh        : "") +
        (motifEnrich_vars.min_peaks          ? " min_peaks="          + motifEnrich_vars.min_peaks          : "") +
        (motifEnrich_vars.extra              ? " "                    + motifEnrich_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The motifEnrich module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if motifEnrich should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("motifEnrich.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_motifs/motifEnrich.R $motifEnrich_FLAGS
        ""","motifEnrich"
    }
}

