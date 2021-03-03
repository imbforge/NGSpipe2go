diffbind2 = {
    doc title: "diffbind",
        desc:  "Differential binding analysis using Diffbind version 2",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Frank Rühle"

    var subdir : ""
    output.dir = diffbind2_vars.outdir + "/$subdir"

    def DIFFBIND_FLAGS =
        (diffbind2_vars.targets           ? " targets="           + diffbind2_vars.targets             : "") +
        (diffbind2_vars.contrasts         ? " contrasts="         + diffbind2_vars.contrasts           : "") +
        (diffbind2_vars.cwd               ? " cwd="               + diffbind2_vars.cwd                 : "") +
        (diffbind2_vars.outdir            ? " out="               + diffbind2_vars.outdir + "/$subdir" : "") +        
        (diffbind2_vars.bams              ? " bams="              + diffbind2_vars.bams                : "") +
        (diffbind2_vars.peaks             ? " peaks="             + diffbind2_vars.peaks + "/$subdir"  : "") +
        (diffbind2_vars.fragsize          ? " fragsize="          + diffbind2_vars.fragsize            : "") +
        (diffbind2_vars.summits           ? " summits="           + diffbind2_vars.summits             : "") +
        (diffbind2_vars.filter            ? " filter="            + diffbind2_vars.filter              : "") +
        (diffbind2_vars.substractControl  ? " substractControl="  + diffbind2_vars.substractControl    : "") +
        (diffbind2_vars.analysisMethod    ? " analysisMethod="    + diffbind2_vars.analysisMethod      : "") +
        (diffbind2_vars.librarySize       ? " librarySize="       + diffbind2_vars.librarySize         : "") +
        (diffbind2_vars.tagwisedispersion ? " tagwiseDispersion=" + diffbind2_vars.tagwisedispersion   : "") + 
        (diffbind2_vars.fdr_threshold     ? " fdr_threshold="     + diffbind2_vars.fdr_threshold       : "") +
        (diffbind2_vars.fold              ? " fold="              + diffbind2_vars.fold                : "") +
        (diffbind2_vars.annotate          ? " annotate="          + diffbind2_vars.annotate            : "") +
        (diffbind2_vars.tss               ? " tss="               + diffbind2_vars.tss                 : "") +
        (diffbind2_vars.txdb              ? " txdb="              + diffbind2_vars.txdb                : "") +
        (diffbind2_vars.annodb            ? " annodb="            + diffbind2_vars.annodb              : "") +
        (diffbind2_vars.genomedb          ? " genomedb="          + diffbind2_vars.genomedb            : "") +
        (diffbind2_vars.paired            ? " pe="                + diffbind2_vars.paired              : "") +
        (diffbind2_vars.extra             ? " "                   + diffbind2_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("R", "3.6.0", "lmod") 
    def PREAMBLE = get_preamble("diffbind")

    // run the chunk
    produce("diffbind.RData", "diffbind.xlsx", "diffbind.rds") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/diffbind/diffbind2.R $DIFFBIND_FLAGS
        ""","diffbind2"
    }

    forward input
}

