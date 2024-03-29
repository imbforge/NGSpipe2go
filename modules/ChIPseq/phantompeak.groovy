phantompeak = {
    doc title: "Phantompeak QC  plot",
        desc:  "Phantompeak",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    var subdir : ""
    output.dir = phantompeak_vars.outdir + "/$subdir" 

    def PHANTOMPEAK_FLAGS =
        (phantompeak_vars.minshift ? " " + phantompeak_vars.minshift : "") +
        (phantompeak_vars.maxshift ? " " + phantompeak_vars.maxshift : "") +
        (phantompeak_vars.binsize  ? " " + phantompeak_vars.binsize  : "") +
        (phantompeak_vars.readlen  ? " " + phantompeak_vars.readlen  : "") +
        (phantompeak_vars.threads  ? " " + phantompeak_vars.threads  : "") +
        (phantompeak_vars.extra    ? " " + phantompeak_vars.extra    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bam") to("_phantompeak.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/ENCODEqc/phantompeak.R $input \$(basename $input.prefix) $PHANTOMPEAK_FLAGS &&

            mv \$(basename $input.prefix)_phantompeak.* $output.dir
        ""","phantompeak"
    }

    forward input
}

