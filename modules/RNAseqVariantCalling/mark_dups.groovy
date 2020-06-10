MarkDups = {
    doc title: "MarkDups",
        desc: "Call picard tools to mark with/without removing duplicated reads from a bam file",
        constraints: "Picard tools version >= 1.141"
        author: "Sergi Sayols, modified by Antonio Domingues"

    output.dir = MarkDups_vars.outdir
    def MarkDups_FLAGS =
        " REMOVE_DUPLICATES=" + (MarkDups_vars.remove_dups   ? "TRUE" : "FALSE") +
        " CREATE_INDEX="      + (MarkDups_vars.index         ? "TRUE" : "FALSE") +
        " ASSUME_SORTED="     + (MarkDups_vars.assume_sorted ? "TRUE" : "FALSE") +
        (MarkDups_vars.validation ? " VALIDATION_STRINGENCY=" + MarkDups_vars.validation : "") +
        (MarkDups_vars.extra      ? " "                       + MarkDups_vars.extra      : "")

    def TOOL_ENV = prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble("MarkDups")

    transform(".rg.bam") to (".rg.duprm.bam"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${MarkDups_vars.java_flags} -jar \${picard} MarkDuplicates $MarkDups_FLAGS INPUT=$input OUTPUT=$output METRICS_FILE=${input.prefix}_dupmetrics.tsv
        ""","MarkDups"
    }
}
