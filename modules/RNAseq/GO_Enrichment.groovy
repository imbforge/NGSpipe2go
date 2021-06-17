GO_Enrichment = {
    doc title: "GO_Enrichment",
        desc: "Gene Ontology enrichment analysis",
        constraints: "",
        bpipe_version: "",
        author: ""

    output.dir = GO_Enrichment_vars.outdir
    def GO_Enrichment_FLAGS =
        (GO_Enrichment_vars.log2fold ? " log2Fold="     + GO_Enrichment_vars.log2fold : "" ) +
        (GO_Enrichment_vars.padj     ? " padj="         + GO_Enrichment_vars.padj     : "" ) +
        (GO_Enrichment_vars.org      ? " organism="     + GO_Enrichment_vars.org      : "" ) +
        (GO_Enrichment_vars.univ     ? " univ="         + GO_Enrichment_vars.univ     : "" ) +
        (GO_Enrichment_vars.type     ? " type="         + GO_Enrichment_vars.type     : "" ) +
        (GO_Enrichment_vars.category ? " plotCategory=" + GO_Enrichment_vars.category : "" ) +
        (GO_Enrichment_vars.outdir   ? " out="          + GO_Enrichment_vars.outdir   : "" ) +
        (GO_Enrichment_vars.cores    ? " cores="        + GO_Enrichment_vars.cores    : "" ) +
        (GO_Enrichment_vars.extra    ? " "              + GO_Enrichment_vars.extra    : "" ) 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(module:"GO_Enrichment", branch:branch, branch_outdir:"")

    transform(".RData") to("_GO.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            touch $output &&
            Rscript ${PIPELINE_ROOT}/tools/GO_Enrichment/GO_Enrichment.R rData=$input $GO_Enrichment_FLAGS &&
            if [ \$? -ne 0 ]; then rm $output; fi;
        ""","GO_Enrichment"
    }
}
