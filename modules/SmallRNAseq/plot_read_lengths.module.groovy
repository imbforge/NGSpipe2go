PlotReadLengths = {
    doc title: "PlotReadLengths",
    desc: "Plots the read length distribution of libraries (fastq).",
    constraints: "none",
    author: "António Domingues"

    output.dir = PLOT_READ_LENGTH_OUTDIR
    OUT_PLOT_DIR = PLOT_READ_LENGTH_OUTDIR + "/" + "plots"
    OUT_PLOT = OUT_PLOT_DIR + "/PercentageReadsLengthDistribution.pdf"

    produce(OUT_PLOT){
        
        exec """
            module load R/${R_VERSION} &&
            Rscript $PLOT_READ_LENGTH_TOOL_PATH $PLOT_READ_LENGTH_OUTDIR

      ""","PlotReadLengths"
    }
}
