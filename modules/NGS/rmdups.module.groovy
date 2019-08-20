// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/rmdups.vars.groovy"

RmDups = {
    doc title: "MarkDups",
        desc:  "Call picard tools to mark with/without removing duplicated reads from a bam file",
        constraints: "Picard tools version >= 1.141",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir=MAPPED
    def JAVA_FLAGS  = "-Xmx" + MARKDUPS_MAXMEM + "m"
    def MARKDUPS_FLAGS  = "REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE"

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble("RmDups")

    transform(".bam") to (".duprm.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java $JAVA_FLAGS -jar \${picard} MarkDuplicates $MARKDUPS_FLAGS INPUT=$input OUTPUT=$output METRICS_FILE=${input.prefix}_dupmetrics.tsv
        ""","RmDups"
    }
}

