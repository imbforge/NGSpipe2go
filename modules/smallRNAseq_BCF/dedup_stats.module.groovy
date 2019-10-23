// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/dedup_stats.vars.groovy"

DedupStats = {
    doc title: "Statistics of unique reads",
        desc:  "Counts the number of reads in the original reads file, and after PCR duplicate removal, and plots results",
        author: "Antonio Domingues, Anke Busch"

    output.dir = DedupStats_vars.plotdir

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("DedupStats")

    produce("dedupReads.pdf", "dedupReads.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/dedup_stats.R ${DedupStats_vars.datadir} ${DedupStats_vars.plotdir} ${DedupStats_vars.sample_prefix}
        ""","DedupStats"
    }
}
