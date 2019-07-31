load MODULE_FOLDER + "NGS/bamindexer.vars.groovy"

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

    transform(".bam\$") to(".bam.bai") {
        exec """
            module load samtools/${SAMTOOLS_VERSION} &&

            samtools index $input
        ""","BAMindexer"
    }

    forward input
}

