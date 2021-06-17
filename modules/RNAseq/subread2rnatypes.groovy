subread2rnatypes = {
    doc title: "subread2rnatypes",
        desc:  "Counting gene biotypes in features with featureCounts of the subread package",
        constraints: "Default: strand specific counting.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = subread2rnatypes_vars.outdir

    def RNATYPES_FLAGS =
        "-F GTF --donotsort " +
        (subread2rnatypes_vars.genesgtf   ? " -a " + subread2rnatypes_vars.genesgtf   : "") +
        (subread2rnatypes_vars.feature    ? " -t " + subread2rnatypes_vars.feature    : "") +
        (subread2rnatypes_vars.accumulate ? " -g " + subread2rnatypes_vars.accumulate : "") +
        (subread2rnatypes_vars.threads    ? " -T " + subread2rnatypes_vars.threads    : "") +
        (subread2rnatypes_vars.paired     ? " -p "                                    : "") +
        (subread2rnatypes_vars.extra      ? " "    + subread2rnatypes_vars.extra      : "") +
        (subread2rnatypes_vars.stranded == "no" ? " -s0 " : (subread2rnatypes_vars.stranded == "yes" ? " -s1 " : " -s2 "))

    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])
    def PREAMBLE = get_preamble(module:"subread2rnatypes", branch:branch, branch_outdir:"")

    // run the chunk
    transform(".bam") to ("_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            base=`basename $input` &&
            if [[ "${subread2rnatypes_vars.paired}" == "true" ]]; then            
                echo "We are resorting and doing the repair\n" &&
                repair -i $input -T ${subread2rnatypes_vars.threads} -o \${TMP}/\${base} &&
                featureCounts $RNATYPES_FLAGS -o \${TMP}/\${base} \${TMP}/\${base} 2> ${output.prefix}_rnatypeslog.stderr;
            else
                featureCounts $RNATYPES_FLAGS -o \${TMP}/\${base} $input 2> ${output.prefix}_rnatypeslog.stderr;
            fi &&
            cut -f1,6,7 \${TMP}/\${base} > $output &&
            rm \${TMP}/\${base}
        ""","subread2rnatypes"
    }
}

