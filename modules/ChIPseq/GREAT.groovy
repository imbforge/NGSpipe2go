GREAT = {
    doc title: "GREAT",
        desc: "Genomic Regions Enrichment Analysis",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino"

    var subdir : ""
    output.dir = GREAT_vars.outdir + "/$subdir"

    def GREAT_FLAGS =
        (GREAT_vars.files      ? " peakData="       + GREAT_vars.files + "/$subdir"  : "") +
        (GREAT_vars.targets    ? " targets="        + GREAT_vars.targets             : "") +
        (GREAT_vars.outdir     ? " out="            + GREAT_vars.outdir + "/$subdir" : "") +
        (GREAT_vars.padj       ? " padj="           + GREAT_vars.padj                : "") +
        (GREAT_vars.nterms     ? " nterms="         + GREAT_vars.nterms              : "") +
        (GREAT_vars.db         ? " db="             + GREAT_vars.db                  : "") +
        (GREAT_vars.upstream   ? " adv_upstream="   + GREAT_vars.upstream            : "") +
        (GREAT_vars.downstream ? " adv_downstream=" + GREAT_vars.downstream          : "") +
        (GREAT_vars.extra      ? " "                + GREAT_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("GREAT")

    produce("GREAT.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/GO_Enrichment/GREAT.R $GREAT_FLAGS
        ""","GREAT"
    }
}

