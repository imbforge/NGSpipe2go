BAMindexer = {
    doc title: "BAMindexer",
        desc:  "Call samtools to index a bam file",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Nastasja Kreim"

    def File f = new File(input1)
    output.dir = f.getParent()

    def TOOL_ENV = prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("BAMindexer")

    transform(".bam\$") to(".bam.bai") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            samtools index $input
        ""","BAMindexer"
    }

    forward input
}

