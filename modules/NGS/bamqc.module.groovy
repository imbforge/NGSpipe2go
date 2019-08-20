// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/bamqc.vars.groovy"

BamQC = {
    doc title: "BamQC",
        desc:  "Quality control of bam file",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Giuseppe Petrosino"

    output.dir   = BAMQC_OUTDIR
    def BAMQC_FLAGS = "--extract --quiet"

    def TOOL_ENV = prepare_tool_env("bamqc", tools["bamqc"]["version"], tools["bamqc"]["runenv"])
    def PREAMBLE = get_preamble("BamQC")

    transform(".bam") to ("_bamqc.zip") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bamqc $BamQC_FLAGS -o $output.dir $input
        ""","BamQC"
    }

    forward input
}
