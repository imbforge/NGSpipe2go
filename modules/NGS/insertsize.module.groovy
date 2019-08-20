// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/insertsize.vars.groovy"

InsertSize = {
    doc title: "InsertSize",
        desc:  "Call picard tools create insert size values",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir=INSERTSIZE_OUTDIR
    def INSERTSIZE_FLAGS = INSERTSIZE_OTHER
    def JAVA_FLAGS  = INSERTSIZE_MAXMEM 

    def TOOL_ENV = prepare_tool_env("R", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && "+
                   prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble("InsertSize")

    transform(".bam") to ("_insertsizemetrics.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java $JAVA_FLAGS -jar \${picard} CollectInsertSizeMetrics $INSERTSIZE_FLAGS INPUT=$input OUTPUT=$output HISTOGRAM_FILE=${output.prefix}_hist.pdf
        ""","InsertSize"
    }
}

