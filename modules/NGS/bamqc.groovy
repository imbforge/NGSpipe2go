BamQC = {
    doc title: "BamQC",
        desc:  "Quality control of bam file",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Giuseppe Petrosino"

    output.dir = BamQC_vars.outdir
    def BAMQC_FLAGS = 
        (BamQC_vars.extra ? " " + BamQC_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("bamqc", tools["bamqc"]["version"], tools["bamqc"]["runenv"])
    def PREAMBLE = get_preamble("BamQC")

    transform(".bam") to ("_bamqc.zip") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bamqc $BAMQC_FLAGS -o $output.dir $input
        ""","BamQC"
    }

    forward input
}
