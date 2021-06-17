CombinedStats = {
    doc title: "Summary of all trimming and filtering statistics of raw reads",
        desc:  "Counts the number of reads in the original reads file, after adapter trimming, after quality filtering and after duplicate removal and plots results",
        author: "Anke Busch"

    output.dir = CombinedStats_vars.plotdir

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(module:"CombinedStats", branch:branch, branch_outdir:"")

    produce("allTrimmingStats.pdf", "allTrimmingStats.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/summary_stats.R $CombinedStats_vars.datadir $CombinedStats_vars.plotdir $CombinedStats_vars.sample_prefix
        ""","CombinedStats"
    }
}
