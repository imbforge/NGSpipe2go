FastQC = {
    doc title: "FastQC",
        desc:  "Quality control of input file",
        constraints: "Only supports compressed FASTQ files",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Frank Rühle"

    var subdir : ""
    output.dir = FastQC_vars.outdir + "/$subdir"

    def FASTQC_FLAGS =
        (FastQC_vars.extra ? " " + FastQC_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("fastqc", tools["fastqc"]["version"], tools["fastqc"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform("*.fastq.gz") to ("_fastqc.zip") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            fastqc --extract $FASTQC_FLAGS -o $output.dir $inputs
        ""","FastQC"
    }

    forward inputs
}
