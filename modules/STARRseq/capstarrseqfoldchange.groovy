CapSTARRseq_FoldChange = {
    doc title: "CapSTARRseq_FoldChange",
        desc:  "Differential expression analysis using RNA / input DNA log2 Fold Change",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = CapSTARRseq_FoldChange_vars.outdir
    def CapSTARRseq_FoldChange_FLAGS =
        (CapSTARRseq_FoldChange_vars.targets   ? " targets="   + CapSTARRseq_FoldChange_vars.targets   : "") +
        (CapSTARRseq_FoldChange_vars.contrasts ? " contrasts=" + CapSTARRseq_FoldChange_vars.contrasts : "") +
        (CapSTARRseq_FoldChange_vars.filter    ? " filter="    + CapSTARRseq_FoldChange_vars.filter    : "") +
        (CapSTARRseq_FoldChange_vars.prefix    ? " prefix="    + CapSTARRseq_FoldChange_vars.prefix    : "") +
        (CapSTARRseq_FoldChange_vars.suffix    ? " suffix="    + CapSTARRseq_FoldChange_vars.suffix    : "") +
        (CapSTARRseq_FoldChange_vars.cwd       ? " cwd="       + CapSTARRseq_FoldChange_vars.cwd       : "") +
        (CapSTARRseq_FoldChange_vars.outdir    ? " out="       + CapSTARRseq_FoldChange_vars.outdir    : "") +
        (CapSTARRseq_FoldChange_vars.genes     ? " gtf="       + CapSTARRseq_FoldChange_vars.genes     : "") +
        (CapSTARRseq_FoldChange_vars.pattern   ? " pattern="   + CapSTARRseq_FoldChange_vars.pattern   : "") +
        (CapSTARRseq_FoldChange_vars.FDR       ? " FDR="       + CapSTARRseq_FoldChange_vars.FDR       : "") +
        (CapSTARRseq_FoldChange_vars.FC        ? " FC="        + CapSTARRseq_FoldChange_vars.FC        : "") +
        (CapSTARRseq_FoldChange_vars.extra     ? " "           + CapSTARRseq_FoldChange_vars.extra     : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    produce("CapSTARRseq_FoldChange.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/CapSTARRseq_FC/CapSTARRseq_FC.R $CapSTARRseq_FoldChange_FLAGS
        ""","CapSTARRseq_FoldChange"
    }
}

