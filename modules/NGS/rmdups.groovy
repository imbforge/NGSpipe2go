RmDups = {
    doc title: "RmDups",
        desc:  "Call picard tools to mark with/without removing duplicated reads from a bam file",
        constraints: "Picard tools version >= 1.141. Expects an env var called `picard` with the path to picard's jar",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = RmDups_vars.outdir
    def RMDUPS_FLAGS =
        " REMOVE_DUPLICATES=" + (RmDups_vars.remove_dups   ? "TRUE" : "FALSE") +
        " ASSUME_SORTED="     + (RmDups_vars.assume_sorted ? "TRUE" : "FALSE") +
        (RmDups_vars.extra ? " " + RmDups_vars.extra : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble("RmDups")

    transform(".bam") to (".duprm.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${RmDups_vars.java_flags} -jar \${picard} MarkDuplicates $RMDUPS_FLAGS INPUT=$input OUTPUT=$output METRICS_FILE=${input.prefix}_duprm_dupmetrics.tsv TMP_DIR=\${TMP}
        ""","RmDups"
    }
}

