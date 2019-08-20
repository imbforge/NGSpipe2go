// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/edger.vars.groovy"

DE_edgeR = {
    doc title: "DE_edgeR",
        desc:  "Differential expression analysis using linears models and edgeR",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = DE_edgeR_OUTDIR.replcaeFirst("out=", "")
    def DE_edgeR_FLAGS = DE_edgeR_TARGETS   + " " + 
                         DE_edgeR_CONTRASTS + " " +
                         DE_edgeR_MMATRIX   + " " +
                         DE_edgeR_FILTER    + " " +
                         DE_edgeR_PREFIX    + " " +
                         DE_edgeR_SUFFIX    + " " +
                         DE_edgeR_CWD       + " " +
                         DE_edgeR_ROBUST    + " " +
                         DE_edgeR_GTF       + " " +
                         DE_edgeR_OUTDIR    + "/DE_edgeR " +
                         DE_edgeR_EXTRA

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

