// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/cutadapt_stats.vars.groovy"

CutadaptStats = {
    doc title: "Adapter trimming statistics",
        desc:  "Creates a plot summarizing the amount of reads left after trimming off the adapter",
        author: "Anke Busch"

    output.dir = REMOVE_ADAPTER_PLOTDIR

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    produce(REMOVE_ADAPTER_PLOTDIR + "/trimmedReads.pdf", REMOVE_ADAPTER_PLOTDIR + "/trimmedReads.png") {
        exec """
            ${TOOL_ENV} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/cutadapt_stats.R ${REMOVE_ADAPTER_DATADIR} ${REMOVE_ADAPTER_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}
        ""","CutadaptStats"
    }
}
