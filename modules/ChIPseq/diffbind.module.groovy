// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/diffbind.vars.groovy"

diffbind = {
    doc title: "diffbind",
        desc:  "Differential binding analysis using Diffbind",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = diffbind_vars.outdir
    def DIFFBIND_FLAGS =
        (diffbind_vars.targets           ? " targets="           + diffbind_vars.targets           : "") +
        (diffbind_vars.contrasts         ? " contrasts="         + diffbind_vars.contrasts         : "") +
        (diffbind_vars.cwd               ? " cwd="               + diffbind_vars.cwd               : "") +
        (diffbind_vars.outdir            ? " out="               + diffbind_vars.outdir            : "") +
        (diffbind_vars.bams              ? " bams="              + diffbind_vars.bams              : "") +
        (diffbind_vars.fragsize          ? " fragsize="          + diffbind_vars.fragsize          : "") +
        (diffbind_vars.substractcontrol  ? " substractControl="  + diffbind_vars.substractcontrol  : "") +
        (diffbind_vars.fulllibrarysize   ? " fullLibrarySize="   + diffbind_vars.fulllibrarysize   : "") +
        (diffbind_vars.tagwisedispersion ? " tagwiseDispersion=" + diffbind_vars.tagwisedispersion : "") +
        (diffbind_vars.annotate          ? " annotate="          + diffbind_vars.annotate          : "") +
        (diffbind_vars.tss               ? " tss="               + diffbind_vars.tss               : "") +
        (diffbind_vars.txdb              ? " txdb="              + diffbind_vars.txdb              : "") +
        (diffbind_vars.annodb            ? " annodb="            + diffbind_vars.annodb            : "") +
        (diffbind_vars.paired            ? " pe="                + diffbind_vars.paired            : "") +
        (diffbind_vars.extra             ? " "                   + diffbind_vars.extra             : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("diffbind")

    // run the chunk
    produce("diffbind.pdf", "diffbind.xlsx", "diffbind.rds") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/diffbind/diffbind.R $DIFFBIND_FLAGS
        ""","diffbind"
    }

    forward input
}

