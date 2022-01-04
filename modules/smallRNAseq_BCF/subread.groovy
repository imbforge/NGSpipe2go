subread_count = {
    doc title: "subread_count",
        desc:  "Counting reads in features with feature-count out of the subread package",
        constraints: "Default: strand specific counting.",
        author: "Oliver Drechsel, Antonio Domingues, Anke Busch"

    var subdir : ""
    output.dir = subread_count_vars.outdir + "/$subdir"

    def SUBREAD_FLAGS =
        "--donotsort " +
        (subread_count_vars.threads  ? " -T " + subread_count_vars.threads  : "") +
        (subread_count_vars.genesgtf ? " -a " + subread_count_vars.genesgtf : "") +
        (subread_count_vars.count_multimapping ? " -M "                     : "") +
        (subread_count_vars.count_ambiguous    ? " -O "                     : "") +
        (subread_count_vars.feature  ? " -t " + subread_count_vars.feature  : "") +
        (subread_count_vars.attribute? " -g " + subread_count_vars.attribute: "") +
        (subread_count_vars.extra    ? " "    + subread_count_vars.extra    : "") +
        (subread_count_vars.stranded == "no" ? " -s0 " : (subread_count_vars.stranded == "yes" ? " -s1 " : " -s2 "))
    
    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to (".raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            featureCounts $SUBREAD_FLAGS -o $output $input 2> ${output.prefix}_subreadlog.stderr
        ""","subread_count"
    }
}
