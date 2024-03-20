filter2htseq = {
    doc title: "filter2htseq",
        desc: "filter featureCount output to fit HTSeq format, extract column 1 and 7 as well as skipping the header",
        constraints: "none.",
        author: "Oliver Drechsel, Antonio Domingues, Anke Busch"

    var subdir : ""
    output.dir = filter2htseq_vars.outdir + "/$subdir"

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".raw_readcounts.tsv") to (".readcounts.tsv") {
        exec """
            ${PREAMBLE} &&

            tail -n +3 $input | awk '{print \$1\"\\t\"\$7}' > $output
        ""","filter2htseq"
    }
}
