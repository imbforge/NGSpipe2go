CollectPlots = {
	doc title: "Compile results",
		desc:  "Merges all plots in a single PDF for ease of visualization",
		author: "Antonio Domingues"

   output.dir = COLLECT_OUTDIR

	produce("all_plots.pdf") {

      exec """
         pdfjoin --paper a4paper --rotateoversize false --landscape \
         $COLLECT_OUTDIR/processed_reads/figure/PCRDuplicates.pdf \
         $COLLECT_OUTDIR/mapped/multimapped/figure/totalReads.pdf \
         $COLLECT_OUTDIR/piRNA_quantification/figure/*.pdf \
         $COLLECT_OUTDIR/nucleotide_signature/*/figure/*.NucleotideDistributionOnPiRNA.pdf \
         $COLLECT_OUTDIR/ping-pong/*/figure/*.ppPlot.pdf \
          -o $output

		""","CollectPlots"
	}
}
