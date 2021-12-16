filter_smallRNA_counts = {
    doc title: "filter_smallRNA_counts",
        desc:  "Extract count of a selected type of smallRNAs to separate count table files",
        constraints: "based on subread (featurecounts) run",
        author: "Anke Busch"

    var subdir : ""
    output.dir = filter_smallRNA_counts_vars.outdir + "/$subdir"

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".readcounts.tsv") to ("." + filter_smallRNA_counts_vars.smallrna + ".readcounts.tsv") {

        exec """

            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/smallRNA_BCF/extract_smallRNA.R gtf=$filter_smallRNA_counts_vars.genesgtf input=$input outdir=$output.dir type=$filter_smallRNA_counts_vars.type smallrna=$filter_smallRNA_counts_vars.smallrna

        ""","filter_smallRNA_counts"
    }
}

