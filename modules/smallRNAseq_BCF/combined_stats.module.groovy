// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/combined_stats.vars.groovy"

CombinedStats = {
    doc title: "Summary of all trimming and filtering statistics of raw reads",
        desc:  "Counts the number of reads in the original reads file, after adapter trimming, after quality filtering and after duplicate removal and plots results",
        author: "Anke Busch"

    output.dir = COMBINED_STATS_PLOTDIR

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("CombinedStats")

    produce(COMBINED_STATS_PLOTDIR + "/allTrimmingStats.pdf", COMBINED_STATS_PLOTDIR + "/allTrimmingStats.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/summary_stats.R ${COMBINED_STATS_DATADIR} ${COMBINED_STATS_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}
        ""","CombinedStats"
    }
}
