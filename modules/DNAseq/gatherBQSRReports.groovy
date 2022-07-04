GatherBQSRReports = {
    doc title: "GATK GatherBQSRReports",
        desc:  "This tool gathers scattered BQSR recalibration reports into a single file.",
        constraints: "",
        bpipe_version: "tested with 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = GatherBQSRReports_vars.outdir

    def GatherBQSRReports_vars_FLAGS = 
            (GatherBQSRReports_vars.extra          ? " "         + GatherBQSRReports_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    produce("gatheredBQSR.report") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${GatherBQSRReports_vars.java_flags}" GatherBQSRReports $GatherBQSRReports_vars_FLAGS --tmp-dir \${TMP} -O $output \$(for f in \$(ls ${output.dir}/*.table);do echo " -I " "\$f"; done)
        ""","GatherBQSRReports"
    }
    forward input
}


