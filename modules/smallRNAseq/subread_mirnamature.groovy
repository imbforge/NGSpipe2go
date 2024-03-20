subread_miRNAmature_count = {
    doc title: "subread_miRNAmature_count",
        desc:  "Counting reads on mature miRNAs with featurecount of the subread package",
        constraints: "miRNA gff (from miRBase) needs to be available.",
        author: "Anke Busch"

    var subdir : ""
    output.dir = subread_miRNAmature_count_vars.outdir + "/$subdir"

    def SUBREAD_MIRNAMATURE_FLAGS =
        "--donotsort " +
        (subread_miRNAmature_count_vars.threads  ? " -T " + subread_miRNAmature_count_vars.threads  : "") +
        (subread_miRNAmature_count_vars.genesgff ? " -a " + subread_miRNAmature_count_vars.genesgff : "") +
        (subread_miRNAmature_count_vars.count_multimapping ? " -M "                     : "") +
        (subread_miRNAmature_count_vars.feature  ? " -t " + subread_miRNAmature_count_vars.feature  : "") +
        (subread_miRNAmature_count_vars.attribute? " -g " + subread_miRNAmature_count_vars.attribute: "") +
        (subread_miRNAmature_count_vars.extra    ? " "    + subread_miRNAmature_count_vars.extra    : "") +
        (subread_miRNAmature_count_vars.stranded == "no" ? " -s0 " : (subread_miRNAmature_count_vars.stranded == "yes" ? " -s1 " : " -s2 "))
    
    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to (".miRNAmature.raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            featureCounts $SUBREAD_MIRNAMATURE_FLAGS -o $output $input 2> ${output.prefix}_subreadlog.stderr
        ""","subread_miRNAmature_count"
    }
}
