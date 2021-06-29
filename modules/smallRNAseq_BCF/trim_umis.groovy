TrimUMIs = {
    doc title: "Trim UMIs",
        desc:  """Trims random barcodes that help in the identification of PCR duplicates and are in adapter-removed reads: NNNN-insert-NNNN.""",
        constraints: "Requires seqtk.",
        author: "Antonio Domingues, Anke Busch"

    output.dir = TrimUMIs_vars.outdir

    def TRIMFQ_FLAGS =
        (TrimUMIs_vars.left_trim  ? " -b " + TrimUMIs_vars.left_trim  : "") +
        (TrimUMIs_vars.right_trim ? " -e " + TrimUMIs_vars.right_trim : "")

    def TOOL_ENV = prepare_tool_env("seqtk", tools["seqtk"]["version"], tools["seqtk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".fastq.gz") to (".trimmed.fastq.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            seqtk trimfq $TRIMFQ_FLAGS $input | gzip > $output
        ""","TrimUMIs"
    }
}
