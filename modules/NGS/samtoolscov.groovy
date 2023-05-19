samtoolscov = {
    doc title: "samtoolscov",
        desc:  "Call samtools to generate coverage statistics",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = samtoolscov_vars.outdir
    def SAMTOOLSCOV_FLAGS = 
        (samtoolscov_vars.extra ? " " + samtoolscov_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bam\$") to(".coverage.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            samtools coverage $SAMTOOLSCOV_FLAGS -o $output $input
        ""","samtoolscov"
    }

    forward input
}

