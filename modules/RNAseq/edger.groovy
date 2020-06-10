DE_edgeR = {
    doc title: "DE_edgeR",
        desc:  "Differential expression analysis using linears models and edgeR",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = DE_edgeR_vars.outdir
    def DE_edgeR_FLAGS =
        (DE_edgeR_vars.targets   ? " targets="   + DE_edgeR_vars.targets   : "" ) +
        (DE_edgeR_vars.contrasts ? " contrasts=" + DE_edgeR_vars.contrasts : "" ) +
        (DE_edgeR_vars.mmatrix   ? " mmatrix="   + DE_edgeR_vars.mmatrix   : "" ) +
        (DE_edgeR_vars.filter    ? " filter="    + DE_edgeR_vars.filter    : "" ) +
        (DE_edgeR_vars.prefix    ? " prefix="    + DE_edgeR_vars.prefix    : "" ) +
        (DE_edgeR_vars.suffix    ? " suffix="    + DE_edgeR_vars.suffix    : "" ) +
        (DE_edgeR_vars.cwd       ? " cwd="       + DE_edgeR_vars.cwd       : "" ) +
        (DE_edgeR_vars.robust    ? " robust="    + DE_edgeR_vars.robust    : "" ) +
        (DE_edgeR_vars.gtf       ? " gtf="       + DE_edgeR_vars.gtf       : "" ) +
        (DE_edgeR_vars.outdir    ? " out="       + DE_edgeR_vars.outdir    : "" ) +
        (DE_edgeR_vars.extra     ? " "           + DE_edgeR_vars.extra     : "" )

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("DE_edgeR")

    // run the chunk
    produce("DE_edgeR.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/DE_edgeR/DE_edgeR.R $DE_edgeR_FLAGS
        ""","DE_edgeR"
    }

    forward input
}

