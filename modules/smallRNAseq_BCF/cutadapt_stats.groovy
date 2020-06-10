CutadaptStats = {
    doc title: "Adapter trimming statistics",
        desc:  "Creates a plot summarizing the amount of reads left after trimming off the adapter",
        author: "Anke Busch"

    output.dir = CutadaptStats_vars.plotdir

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("CutadaptStats")

    produce("trimmedReads.pdf", "trimmedReads.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/cutadapt_stats.R ${CutadaptStats_vars.datadir} ${CutadaptStats_vars.plotdir} ${CutadaptStats_vars.sample_prefix}
        ""","CutadaptStats"
    }
}
