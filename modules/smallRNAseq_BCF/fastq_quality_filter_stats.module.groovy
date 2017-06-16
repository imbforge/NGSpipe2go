FastQQualityFilterStats = {
		doc title: "Statistics of quality filtered reads",
		desc:  "Creates a plot summarizing the numbers of removed and kept reads due to low and high qualities, respectively",
		author: "Anke Busch"

		output.dir = REMOVE_LOWQUAL_PLOTDIR

		produce(REMOVE_LOWQUAL_PLOTDIR + "/qualityFilteredReads.pdf", REMOVE_LOWQUAL_PLOTDIR + "/qualityFilteredReads.png") {

			exec """

            module load R/${R_VERSION} && 

			        Rscript ${QUALITY_FILTER_PLOT_TOOL} ${REMOVE_LOWQUAL_DATADIR} ${REMOVE_LOWQUAL_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX} 

			""","FastQQualityFilterStats"
		}
}
