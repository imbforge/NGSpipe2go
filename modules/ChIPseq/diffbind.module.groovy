// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/diffbind.vars.groovy"

diffbind = {
    doc title: "diffbind",
        desc:  "Differential binding analysis using Diffbind",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = RESULTS + "/diffbind"
    def DIFFBIND_FLAGS = DIFFBIND_TARGETS   + " " + 
                         DIFFBIND_CONTRASTS + " " +
                         DIFFBIND_CWD       + " " +
                         DIFFBIND_OUTDIR    + " " +
                         DIFFBIND_BAMS      + " " +
                         DIFFBIND_FRAGSIZE  + " " +
                         DIFFBIND_SUBSTRACTCONTROL  + " " +
                         DIFFBIND_FULLLIBRARYSIZE   + " " +
                         DIFFBIND_TAGWISEDISPERSION + " " +
                         DIFFBIND_ANNOTATE  + " " +
                         DIFFBIND_TSS       + " " +
                         DIFFBIND_TXDB      + " " +
                         DIFFBIND_ANNODB    + " " +
                         DIFFBIND_PAIRED    + " " +
                         DIFFBIND_EXTRA

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    // run the chunk
    produce("diffbind.pdf", "diffbind.xlsx", "diffbind.rds") {
        exec """
            ${TOOL_ENV} &&

            Rscript ${PIPELINE_ROOT}/tools/diffbind/diffbind.R $DIFFBIND_FLAGS
        ""","diffbind"
    }

    forward input
}

