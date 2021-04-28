subread_count = {
    doc title: "subread_count_se",
        desc:  "Counting reads in features with feature-count out of the subread package",
        constraints: """Default: strand specific counting.""",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir  = subread_count_vars.outdir
    def SUBREAD_FLAGS =
        "--donotsort " +
        (subread_count_vars.threads  ? " -T " + subread_count_vars.threads  : "") +
        (subread_count_vars.genesgtf ? " -a " + subread_count_vars.genesgtf : "") +
        (subread_count_vars.feature  ? " -t " + subread_count_vars.feature  : "") +
        (subread_count_vars.paired   ? " -p "                               : "") +
        (subread_count_vars.extra    ? " "    + subread_count_vars.extra    : "") +
        (subread_count_vars.stranded == "no" ? " -s0 " : (subread_count_vars.stranded == "yes" ? " -s1 " : " -s2 "))

    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble("subread_count")

    // run the chunk
    transform(".bam") to (".raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            base=`basename $input` &&
            if [[ "${subread_count_vars.paired}" == "true" ]]; 
            then
                echo "We are resorting and doing the repair\n" &&
                repair -i $input -T ${subread_count_vars.threads} -o \${TMP}/\${base} &&
                featureCounts $SUBREAD_FLAGS -o $output \${TMP}/\${base} 2> ${output.prefix}_subreadlog.stderr;
            else
                featureCounts $SUBREAD_FLAGS -o $output $input 2> ${output.prefix}_subreadlog.stderr;
            fi
        ""","subread_count"
    }
}
