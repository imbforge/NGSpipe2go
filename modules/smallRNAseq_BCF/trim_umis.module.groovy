// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/trim_umis.vars.groovy"

TrimUMIs = {
    doc title: "Trim UMIs",
        desc:  """Trims random barcodes that help in the identification of PCR duplicates and are in adapter-removed reads: NNNN-insert-NNNN.""",
        constraints: "Requires seqtk.",
        author: "Antonio Domingues, Anke Busch"

    output.dir = TRIM_OUTDIR
    def SEQTK_FLAGS = " -l " + LEFT_TRIM +
                      " -b " + RIGHT_TRIM

    def TOOL_ENV = prepare_tool_env("seqtk", tools["seqtk"]["version"], tools["seqtk"]["runenv"])
    def PREAMBLE = get_preamble("TrimUMIs")

    transform(".fastq.gz") to (".trimmed.fastq.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            seqtk trimfq -b ${LEFT_TRIM} -e ${RIGHT_TRIM} $input | gzip > $output
        ""","TrimUMIs"
    }
}
