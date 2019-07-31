load MODULE_FOLDER + "SmallRNAseq/plot_sequence_bias.vars.groovy"

PlotSequenceBias = {
    doc title: "PlotSequenceBias",
    desc: "Plots the read length distribution and first nucleotide frequency for mapped reads.",
    constraints: "none",
    author: "António Domingues"

    output.dir = PLOT_BIAS_OUTDIR
    def OUT_PLOT_DIR = PLOT_BIAS_OUTDIR + "/" + "figure"
    def OUT_PLOT = OUT_PLOT_DIR + "/nucleotide_bias_read_length.normalized.pdf"

    from("*.nuc_bias.txt") produce(OUT_PLOT){
        exec """
            module load R/${R_VERSION} &&
            Rscript $PLOT_SMALL_RNA_TOOL_PATH --inputs $inputs --libsizes ${PLOT_BIAS_MAPPED} --out $output.dir 

      ""","PlotSequenceBias"
    }
}
