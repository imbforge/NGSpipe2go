// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup_stats.vars.groovy"

DedupStats = {
    doc title: "Statistics of unique reads",
        desc:  "Counts the number of reads in the original reads file, and after PCR duplicate removal, and plots results",
        author: "Antonio Domingues, Anke Busch"

    output.dir = REMOVE_DUP_PLOTDIR

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("DedupStats")

    produce(REMOVE_DUP_PLOTDIR + "/dedupReads.pdf", REMOVE_DUP_PLOTDIR + "/dedupReads.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/dedup_stats.R ${REMOVE_DUP_DATADIR} ${REMOVE_DUP_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}
        ""","DedupStats"
    }
}
