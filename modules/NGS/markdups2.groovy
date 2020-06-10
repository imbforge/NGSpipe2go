MarkDups2 = {
    doc title: "MarkDups2",
        desc:  "Call bamUtil dedup tool to mark with/without removing duplicated reads from a bam file",
        constraints: "bamUtil tool version >= 1.0.13",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Giuseppe Petrosino"

    output.dir=MarkDups2_vars.outdir

    def TOOL_ENV = prepare_tool_env("bamutil", tools["bamutil"]["version"], tools["bamutil"]["runenv"])
    def PREAMBLE = get_preamble("MarkDups2")

    transform(".bam") to (".dupmarked.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bam dedup --in $input --out $output --log ${input.prefix}_dupmetrics.log --noPhoneHome
        ""","MarkDups2"
    }
}
