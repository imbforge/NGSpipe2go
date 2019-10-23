// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mapping_stats.vars.groovy"

MappingStats = {
    doc title: "Statistics of mapping efficiency",
        desc:  "Counts the number of reads in the mapped bam, including total, unique, and mapped. Returns a plot of the results.",
        constraints: "Bam files produced by Bowtie 1.x. Might not work for other mappers/versions.",
        author: "Antonio Domingues, Anke Busch"

    output.dir = MappingStats_vars.plotdir

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("MappingStats")

    produce("totalReads.pdf", "totalReads.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/mapping_stats/mapping_stats_bowtie1_BCF.R ${MappingStats_vars.datadir} ${MappingStats_vars.plotdir} ${MappingStats_vars.sample_prefix}
        ""","MappingStats"
    }
}
