load MODULE_FOLDER + "smallRNAseq_BCF/filter2htseq.vars.groovy"

Filter2HTSeq = {
    doc title: "Filter2HTSeq",
        desc: "filter featureCount output to fit HTSeq format, extract column 1 and 7 as well as skipping the header",
        constraints: "none.",
        author: "Oliver Drechsel, Antonio Domingues, Anke Busch"

    output.dir = SUBREAD_OUTDIR

    transform(".raw_readcounts.tsv") to (".readcounts.tsv") {
        exec """
            tail -n +3 $input | awk '{print \$1\"\\t\"\$7}' > $output
        ""","Filter2HTSeq"
    }
}
