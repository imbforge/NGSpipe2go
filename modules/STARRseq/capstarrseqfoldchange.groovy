CapSTARRseq_FoldChange = {
    doc title: "CapSTARRseq_FoldChange",
        desc:  "Differential expression analysis using RNA / input DNA log2 Fold Change",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Martin Oti, Sergi Sayols"

    output.dir = CapSTARRseq_FoldChange_vars.outdir + "/$subdir"
    def CapSTARRseq_FoldChange_FLAGS =
        (CapSTARRseq_FoldChange_vars.targets   ? " targets="   + CapSTARRseq_FoldChange_vars.targets             : "") +
        (CapSTARRseq_FoldChange_vars.cwd       ? " cwd="       + CapSTARRseq_FoldChange_vars.cwd + "/$subdir"    : "") +
        (CapSTARRseq_FoldChange_vars.outdir    ? " out="       + CapSTARRseq_FoldChange_vars.outdir + "/$subdir" : "") +
        (CapSTARRseq_FoldChange_vars.gtf       ? " gtf="       + CapSTARRseq_FoldChange_vars.gtf                 : "") +
        (CapSTARRseq_FoldChange_vars.pattern   ? " pattern="   + CapSTARRseq_FoldChange_vars.pattern             : "") +
        (CapSTARRseq_FoldChange_vars.extra     ? " "           + CapSTARRseq_FoldChange_vars.extra               : "") 

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

