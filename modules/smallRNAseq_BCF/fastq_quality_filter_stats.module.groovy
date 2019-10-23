// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter_stats.vars.groovy"

FastQQualityFilterStats = {
    doc title: "Statistics of quality filtered reads",
        desc:  "Creates a plot summarizing the numbers of removed and kept reads due to low and high qualities, respectively",
        author: "Anke Busch"

    output.dir = FastQQualityFilterStats_vars.plotdir

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("FastQQualityFilterStats")

    produce("qualityFilteredReads.pdf", "qualityFilteredReads.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/fastq_stats/fastq_quality_filter_stats.R ${FastQQualityFilterStats_vars.datadir} ${FastQQualityFilterStats_vars.plotdir} ${FastQQualityFilterStats_vars.sample_prefix} 
        ""","FastQQualityFilterStats"
    }
}
