// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/bamindexer.vars.groovy"

BAMindexer = {
    doc title: "BAMindexer",
        desc:  "Call samtools to index a bam file",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Nastasja Kreim"

    def OUTPUTDIR = input1
    int path_index = OUTPUTDIR.lastIndexOf("/")
    OUTPUTDIR = OUTPUTDIR.substring(0,path_index)
    output.dir = OUTPUTDIR

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

