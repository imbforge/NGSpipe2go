motifFootprinting = {
    doc title: "motifFootprinting",
        desc:  "Perform motif footprinting with enriched TF motifs",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = motifFootprinting_vars.outdir
    
    def motifFootprinting_FLAGS =
        (motifFootprinting_vars.outdir             ? " outdir="             + motifFootprinting_vars.outdir             : "") +
        (motifFootprinting_vars.project            ? " project="            + motifFootprinting_vars.project            : "") +
        (motifFootprinting_vars.res                ? " res="                + motifFootprinting_vars.res                : "") +
        (motifFootprinting_vars.db                 ? " db="                 + motifFootprinting_vars.db                 : "") +
        (motifFootprinting_vars.motifEnrich_dir    ? " motifEnrich_dir="    + motifFootprinting_vars.motifEnrich_dir    : "") +
        (motifFootprinting_vars.assay              ? " assay="              + motifFootprinting_vars.assay              : "") +
        (motifFootprinting_vars.motifsPerContrast  ? " motifsPerContrast="  + motifFootprinting_vars.motifsPerContrast  : "") +
        (motifFootprinting_vars.motifsByName       ? " motifsByName="       + motifFootprinting_vars.motifsByName       : "") +
        (motifFootprinting_vars.inPeaks            ? " inPeaks="            + motifFootprinting_vars.inPeaks            : "") +
        (motifFootprinting_vars.upstream           ? " upstream="           + motifFootprinting_vars.upstream           : "") +
        (motifFootprinting_vars.downstream         ? " downstream="         + motifFootprinting_vars.downstream         : "") +
        (motifFootprinting_vars.groupPlotBy        ? " groupPlotBy="        + motifFootprinting_vars.groupPlotBy        : "") +
        (motifFootprinting_vars.splitPlotBy        ? " splitPlotBy="        + motifFootprinting_vars.splitPlotBy        : "") +
        (motifFootprinting_vars.extra              ? " "                    + motifFootprinting_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The motifFootprinting module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if motifFootprinting should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("motifFootprinting.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_motifs/motifFootprinting.R $motifFootprinting_FLAGS
        ""","motifFootprinting"
    }
}

