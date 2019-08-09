// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/GREAT.vars.groovy"

GREAT = {
    doc title: "GREAT",
        desc: "Genomic Regions Enrichment Analysis",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino"

    output.dir=GREAT_OUTDIR.replaceFirst("out=","")
    def GREAT_FLAGS = GREAT_FILES + " " +
        GREAT_TARGETS + " " + 
        GREAT_OUTDIR + " " +
        GREAT_PADJ + " " +
        GREAT_NTERMS + " " +
        GREAT_DB + " " +
        GREAT_UPSTREAM + " " +
        GREAT_DOWNSTREAM + " " +
        GREAT_EXTRA

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    produce("GREAT.RData") {
        exec """
           ${TOOL_ENV} &&

           Rscript ${PIPELINE_ROOT}/tools/GO_Enrichment/GREAT.R $GREAT_FLAGS
        ""","GREAT"
    }
}

