filter_miRNA_counts = {
    doc title: "filter_miRNA_counts",
        desc:  "Extract count of miRNAs to separate count table files",
        constraints: "based on subread (featurecounts) run",
        author: "Anke Busch"

    var subdir : ""
    output.dir = filter_miRNA_counts_vars.outdir + "/$subdir"

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".readcounts.tsv") to (".miRNA.readcounts.tsv") {

        exec """

            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/miRNA/extract_miRNA.R gtf=$filter_miRNA_counts_vars.genesgtf input=$input outdir=$output.dir

        ""","filter_miRNA_counts"
    }
}

