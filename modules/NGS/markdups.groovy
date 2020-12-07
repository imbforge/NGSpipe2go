MarkDups = {
    doc title: "MarkDups",
        desc:  "Call picard tools to mark with/without removing duplicated reads from a bam file",
        constraints: "Picard tools version >= 1.141. Expects an env var called `picard` with the path to picard's jar",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = MarkDups_vars.outdir
    def MARKDUPS_FLAGS =
        " REMOVE_DUPLICATES=" + (MarkDups_vars.remove_dups   ? "TRUE" : "FALSE") +
        " ASSUME_SORTED="     + (MarkDups_vars.assume_sorted ? "TRUE" : "FALSE") +
        (MarkDups_vars.extra ? " " + MarkDups_vars.extra : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble("MarkDups")

    transform(".bam") to (".dupmarked.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${MarkDups_vars.java_flags} -jar \${picard} MarkDuplicates $MARKDUPS_FLAGS INPUT=$input OUTPUT=$output METRICS_FILE=${input.prefix}_dupmarked_dupmetrics.tsv
        ""","MarkDups"
    }
}

