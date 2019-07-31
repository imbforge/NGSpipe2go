load MODULE_FOLDER + "ChIPseq/GREAT.vars.groovy"

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

    produce("GREAT.RData") {
        exec """
           module load R/${R_VERSION} &&

           Rscript ${TOOL_GO}/GREAT.R $GREAT_FLAGS
        ""","GREAT"
    }
}

