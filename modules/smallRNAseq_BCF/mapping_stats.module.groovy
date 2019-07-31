load MODULE_FOLDER + "smallRNAseq_BCF/mapping_stats.vars.groovy"

MappingStats = {
    doc title: "Statistics of mapping efficiency",
        desc:  "Counts the number of reads in the mapped bam, including total, unique, and mapped. Returns a plot of the results.",
        constraints: "Bam files produced by Bowtie 1.x. Might not work for other mappers/versions.",
        author: "Antonio Domingues, Anke Busch"

    output.dir = MAPPING_STATS_PLOTDIR

    produce(MAPPING_STATS_PLOTDIR + "/totalReads.pdf", MAPPING_STATS_PLOTDIR + "/totalReads.png") {
        exec """
            module load R/${R_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&

            Rscript ${MAPPING_STATS_TOOL} ${MAPPING_STATS_DATADIR} ${MAPPING_STATS_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}
        ""","MappingStats"
    }
}
