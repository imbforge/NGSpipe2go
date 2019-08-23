// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/deseq2.vars.groovy"

DE_DESeq2 = {
    doc title: "DE_DESeq2",
        desc:  "Differential expression analysis using linear models and DESeq2",
        constraints: "Only simple contrasts in 1-factor design. Include always the intercept. Always gene filtering",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = DE_DESeq2_vars.outdir
    def DE_DESeq2_FLAGS =
        (DE_DESeq2_vars.targets   ? " targets="   + DE_DESeq2_vars.targets   : "") +
        (DE_DESeq2_vars.contrasts ? " contrasts=" + DE_DESeq2_vars.contrasts : "") +
        (DE_DESeq2_vars.mmatrix   ? " mmatrix="   + DE_DESeq2_vars.mmatrix   : "") +
        (DE_DESeq2_vars.filter    ? " filter="    + DE_DESeq2_vars.filter    : "") +
        (DE_DESeq2_vars.prefix    ? " prefix="    + DE_DESeq2_vars.prefix    : "") +
        (DE_DESeq2_vars.suffix    ? " suffix="    + DE_DESeq2_vars.suffix    : "") +
        (DE_DESeq2_vars.cwd       ? " cwd="       + DE_DESeq2_vars.cwd       : "") +
        (DE_DESeq2_vars.outdir    ? " out="       + DE_DESeq2_vars.outdir    : "") +
        (DE_DESeq2_vars.genes     ? " gtf="       + DE_DESeq2_vars.genes     : "") +
        (DE_DESeq2_vars.extra     ? " "           + DE_DESeq2_vars.extra     : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("DE_DESeq2")

    // run the chunk
    produce("DE_DESeq2.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/DE_DESeq2/DE_DESeq2.R $DE_DESeq2_FLAGS
        ""","DE_DESeq2"
    }
}

