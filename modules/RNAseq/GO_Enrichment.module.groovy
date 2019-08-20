// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/GO_Enrichment.vars.groovy"

GO_Enrichment = {
    doc title: "GO_Enrichment",
        desc: "Gene Ontology enrichment analysis",
        constraints: "",
        bpipe_version: "",
        author: ""

    output.dir = GO_Enrichment_OUTDIR.replaceFirst("out=", "")
    def GO_Enrichment_FLAGS = GO_Enrichment_LOG2FOLD + " " + 
                              GO_Enrichment_PADJ     + " " +
                              GO_Enrichment_ORG      + " " +
                              GO_Enrichment_UNIV     + " " +
                              GO_Enrichment_TYPE     + " " +
                              GO_Enrichment_CATEGORY + " " +
                              GO_Enrichment_OUTDIR   + " " +
                              GO_Enrichment_CORES    + " " +
                              GO_Enrichment_EXTRA

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("GO_Enrichment")

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
