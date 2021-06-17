SubreadCount = {
    doc title: "SubreadCount",
        desc:  "Counting reads in features with feature-count out of the subread package",
        constraints: "Default: strand specific counting.",
        author: "Oliver Drechsel, Antonio Domingues, Anke Busch"

    output.dir  = SubreadCount_vars.outdir
    def SUBREAD_FLAGS =
        "--donotsort " +
        (SubreadCount_vars.threads  ? " -T " + SubreadCount_vars.threads  : "") +
        (SubreadCount_vars.genesgtf ? " -a " + SubreadCount_vars.genesgtf : "") +
        (SubreadCount_vars.ignore_duplicates  ? " --ignoreDup "            : "") +
        (SubreadCount_vars.count_multimapping ? " -M "                     : "") +
        (SubreadCount_vars.extra    ? " "    + SubreadCount_vars.extra    : "") +
        (SubreadCount_vars.stranded == "no" ? " -s0 " : (subread2rnatypes_vars.stranded == "yes" ? " -s1 " : " -s2 "))
    
    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble(module:"SubreadCount", branch:branch, branch_outdir:"")

    // run the chunk
    transform(".bam") to (".raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            featureCounts $SUBREAD_FLAGS -o $output $input 2> ${output.prefix}_subreadlog.stderr
        ""","SubreadCount"
    }
}
