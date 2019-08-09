// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/fastqc.vars.groovy"

FastQC = {
    doc title: "FastQC",
        desc:  "Quality control of input file",
        constraints: "Only supports compressed FASTQ files",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = FASTQC_OUTDIR
    def FASTQC_FLAGS = "--extract --quiet"

    def TOOL_ENV = prepare_tool_env("fastqc", tools["fastqc"]["version"], tools["fastqc"]["runenv"])

    transform(".fastq.gz") to ("_fastqc.zip") {
        exec """
            ${TOOL_ENV} &&

            fastqc $FASTQC_FLAGS -o $output.dir $input
        ""","FastQC"
    }

    forward input
}
