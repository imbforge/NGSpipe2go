diffbind3 = {
    doc title: "diffbind",
        desc:  "Differential binding analysis using Diffbind version 3",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Frank Rühle"

    var subdir : ""
    output.dir = diffbind3_vars.outdir + "/$subdir"

    def DIFFBIND_FLAGS =
        (diffbind3_vars.diffbindversion   ? " diffbindversion="   + diffbind3_vars.diffbindversion     : "") +
        (diffbind3_vars.targets           ? " targets="           + diffbind3_vars.targets             : "") +
        (diffbind3_vars.contrasts         ? " contrasts="         + diffbind3_vars.contrasts           : "") +
        (diffbind3_vars.cwd               ? " cwd="               + diffbind3_vars.cwd                 : "") +
        (diffbind3_vars.outdir            ? " out="               + diffbind3_vars.outdir + "/$subdir" : "") +        
        (diffbind3_vars.bams              ? " bams="              + diffbind3_vars.bams                : "") +
        (diffbind3_vars.peaks             ? " peaks="             + diffbind3_vars.peaks + "/$subdir"  : "") +
        (diffbind3_vars.fragsize          ? " fragsize="          + diffbind3_vars.fragsize            : "") +
        (diffbind3_vars.blacklist         ? " blacklist="         + diffbind3_vars.blacklist           : "") +
        (diffbind3_vars.greylist          ? " greylist="          + diffbind3_vars.greylist            : "") +
        (diffbind3_vars.summits           ? " summits="           + diffbind3_vars.summits             : "") +
        (diffbind3_vars.filter            ? " filter="            + diffbind3_vars.filter              : "") +
        (diffbind3_vars.analysisMethod    ? " analysisMethod="    + diffbind3_vars.analysisMethod      : "") +
        (diffbind3_vars.librarySize       ? " librarySize="       + diffbind3_vars.librarySize         : "") +
        (diffbind3_vars.normalization     ? " normalization="     + diffbind3_vars.normalization       : "") +
        (diffbind3_vars.substractControl  ? " substractControl="  + diffbind3_vars.substractControl    : "") +
        (diffbind3_vars.conditionColumn   ? " conditionColumn="   + diffbind3_vars.conditionColumn     : "") +
        (diffbind3_vars.design            ? " design="            + diffbind3_vars.design              : "") +
        (diffbind3_vars.fdr_threshold     ? " fdr_threshold="     + diffbind3_vars.fdr_threshold       : "") +
        (diffbind3_vars.fold              ? " fold="              + diffbind3_vars.fold                : "") +
        (diffbind3_vars.annotate          ? " annotate="          + diffbind3_vars.annotate            : "") +
        (diffbind3_vars.tss               ? " tss="               + diffbind3_vars.tss                 : "") +
        (diffbind3_vars.txdb              ? " txdb="              + diffbind3_vars.txdb                : "") +
        (diffbind3_vars.annodb            ? " annodb="            + diffbind3_vars.annodb              : "") +
        (diffbind3_vars.genomedb          ? " genomedb="          + diffbind3_vars.genomedb            : "") +
        (diffbind3_vars.paired            ? " pe="                + diffbind3_vars.paired              : "") +
        (diffbind3_vars.extra             ? " "                   + diffbind3_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("R", "4.0.3", "lmod")
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    produce("diffbind.RData", "diffbind.xlsx", "diffbind.rds") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/diffbind/diffbind3.R $DIFFBIND_FLAGS
        ""","diffbind3"
    }

    forward input
}

