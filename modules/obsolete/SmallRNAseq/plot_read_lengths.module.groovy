load MODULE_FOLDER + "SmallRNAseq/plot_read_lengths.vars.groovy"

PlotReadLengths = {
    doc title: "PlotReadLengths",
    desc: "Plots the read length distribution of libraries (fastq).",
    constraints: "none",
    author: "António Domingues"

    output.dir = PLOT_READ_LENGTH_OUTDIR
    def OUT_PLOT_DIR = PLOT_READ_LENGTH_OUTDIR + "/" + "plots"
    def OUT_PLOT = OUT_PLOT_DIR + "/PercentageReadsLengthDistribution.pdf"

    produce(OUT_PLOT){
        exec """
            module load R/${R_VERSION} &&
            Rscript $PLOT_READ_LENGTH_TOOL_PATH $PLOT_READ_LENGTH_OUTDIR

      ""","PlotReadLengths"
    }
}
