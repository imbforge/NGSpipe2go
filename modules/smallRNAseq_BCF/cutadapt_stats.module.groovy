CutadaptStats = {
		doc title: "Adapter trimming statistics",
		desc:  "Creates a plot summarizing the amount of reads left after trimming off the adapter",
		author: "Anke Busch"

		output.dir = REMOVE_ADAPTER_PLOTDIR

		produce(REMOVE_ADAPTER_PLOTDIR + "/trimmedReads.pdf", REMOVE_ADAPTER_PLOTDIR + "/trimmedReads.png") {

			exec """
				export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
				source ${TOOL_R}/env.sh &&

				echo 'VERSION INFO'  1>&2 &&
				echo \$(Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 &&
				echo '/VERSION INFO'  1>&2 &&

        	 		Rscript ${CUTADAPT_PLOT_TOOL} ${REMOVE_ADAPTER_DATADIR} ${REMOVE_ADAPTER_PLOTDIR} ${ESSENTIAL_SAMPLE_PREFIX}

			""","CutadaptStats"
		}
}
